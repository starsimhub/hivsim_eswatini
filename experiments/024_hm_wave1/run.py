"""Exp 024 — History matching, wave 1.

1000 Latin-hypercube points, 1 replicate each, N = 10,000, bayes_linear
emulator, implausibility threshold 4.0. Emulates one macro feature:
prevalence 15-49 averaged over the three PHIA years.

See README.md for why each sigma is what it is -- that is the decision most
likely to make or break the wave, and it is argued there rather than here.

Usage
  python run.py                 # run wave 1 (~4.5 h on 8 workers)
  python run.py --n_samples 20  # smoke test
  python run.py --resume        # continue an interrupted run
  python run.py --analyse_only  # re-read the checkpoint and re-report

Outputs
  outputs/sims/point_{hash}.parquet  per-point output, keyed by parameter
                                   values so the cache is design-independent
  outputs/hm/wave1/                package diagnostics: pairplot, convergence,
                                   zscores_vs_targets, constrained_dims,
                                   nroy_samples.csv
  outputs/observations.csv         every registered target with its sigma
  outputs/nroy_summary.txt         the surviving region
"""

import os
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                  NUMEXPR_NUM_THREADS="1", MKL_NUM_THREADS="1")

import argparse
import hashlib
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import sciris as sc

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
OUT_DIR = EXP_DIR / "outputs"
SIM_DIR, HM_DIR = OUT_DIR / "sims", OUT_DIR / "hm"
for d in (OUT_DIR, SIM_DIR, HM_DIR):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
from hm_shim import hm, DEFAULT_EMULATOR  # noqa: E402  (stubs gpflow first)
from run_sims import make_sim  # noqa: E402
from analyzers import PopByAgeSex, Cascade  # noqa: E402
from standard_figures import load_targets, targets_15_49  # noqa: E402

N_AGENTS = 10_000
STOP = 2026
SIM_SEED = 1
N_SAMPLES = 1000
THRESHOLD = 4.0          # raised from 3.0: known structural misspecification
RANDOM_SEED = 20260903
PHIA_YEARS = (2007, 2011, 2016)

# --- Parameter box. Bounds only -- HM has no prior density. -------------------
# Log-bounded where 023 sampled on a log scale; the simulator un-transforms.
BOUNDS = {
    "log_beta_m2f":        (np.log(0.0096), np.log(0.025)),
    "log_rel_beta_f2m":    (np.log(0.15),   np.log(0.60)),
    "log_s_f_young":       (np.log(0.8),    np.log(3.0)),
    "age_gap_shift":       (-2.0,           3.0),
    "log_age_gap_sd_mult": (np.log(0.6),    np.log(1.8)),
    "prop_f0":             (0.45,           0.85),
    "prop_m0":             (0.40,           0.80),
}
FIXED = {"rel_init_prev": 0.2, "conc_mult": 1.0}   # dropped in 023

AGE_DIFF_BASE = {
    "teens": [(7, 3), (6, 3), (5, 1)],
    "young": [(8, 3), (7, 3), (5, 2)],
    "adult": [(8, 3), (7, 3), (5, 2)],
}
CONC_BASE = {"f1_conc": 0.15, "f2_conc": 0.25, "m1_conc": 0.15, "m2_conc": 0.5}

# --- Observation uncertainty ---------------------------------------------------
# Model discrepancy on prevalence. Anchored so the best configuration measured
# to date (022 arm A, 15-49 aggregate bias -0.043) sits at ~2 sigma. A modelling
# judgement, not a measurement -- see README. If wave 1 collapses, revisit this
# before revisiting the model.
SIGMA_DISC_PREV = 0.02
# The deaths down-weighting. HM has no weights, only sigmas. The model reaches
# 7,069 against UNAIDS' 11,000, so 2000 puts that gap at ~2 sigma; a survey-like
# sigma would put it at 8 sigma and rule out the entire box.
UNAIDS_PEAK, SIGMA_DEATHS = 11_000.0, 2000.0
SIGMA_INC_FLOOR = 0.25       # incidence CIs are wide already; floor the narrow ones
SIGMA_RATIO_REL = 0.15       # prevalence ratios, relative
SIGMA_INC_RATIO = 0.35       # F:M incidence ratio, on a value near 2.0

WAVE1_FEATURE = "prev_15_49_all_mean"


def build_observations():
    """Every registered target with its sigma. Only some are emulated per wave."""
    tg = load_targets()
    agg = targets_15_49(tg)
    obs, prov = {}, []

    def add(name, mean, sd, tier, note):
        obs[name] = (float(mean), float(sd))
        prov.append(dict(feature=name, mean=mean, sigma=sd, tier=tier, note=note))

    # --- Tier A: macro ---
    for _, r in agg.iterrows():
        b = tg[(tg.year == r.year) & (tg.sex == r.sex) & tg.age_low.between(15, 45)]
        w = b.Count / b.Count.sum()
        sd_ci = float(np.sqrt((w ** 2 * ((b.ub - b.lb) / 3.92) ** 2).sum()))
        add(f"prev_15_49_{r.sex}_{int(r.year)}", r.phia,
            np.hypot(sd_ci, SIGMA_DISC_PREV), "A",
            f"PHIA CI sigma {sd_ci:.4f} + discrepancy {SIGMA_DISC_PREV}")
    both = agg.groupby("year").phia.mean()
    add(WAVE1_FEATURE, both.mean(), np.hypot(0.005, SIGMA_DISC_PREV), "A",
        "mean of sex-specific prevalence 15-49 over 2007/2011/2016 -- wave 1 feature")
    add("peak_deaths", UNAIDS_PEAK, SIGMA_DEATHS, "A",
        "UNAIDS peak 2004. sigma is the deliberate down-weighting, not a measurement")

    # --- Tier B: first-order shape ---
    for year in PHIA_YEARS:
        f = float(agg[(agg.year == year) & (agg.sex == "f")].phia.iloc[0])
        m = float(agg[(agg.year == year) & (agg.sex == "m")].phia.iloc[0])
        add(f"prev_fm_ratio_{year}", f / m, SIGMA_RATIO_REL * (f / m), "B",
            "female:male prevalence 15-49")
        b = tg[(tg.year == year) & (tg.sex == "f")]
        yo = (b[b.age_low.between(15, 20)].phia.mean()
              / b[b.age_low.between(35, 40)].phia.mean())
        add(f"prev_young_old_f_{year}", yo, SIGMA_RATIO_REL * yo, "B",
            "female 15-24 : 35-44 prevalence")

    inc = pd.read_csv(REPO_ROOT / "calibration_data" / "incidence_by_age_sex.csv")
    for _, r in inc.iterrows():
        sd = max(float((r.ub - r.lb) / 3.92), SIGMA_INC_FLOOR)
        add(f"inc_{r.sex}_{int(r.age_low)}_{int(r.age_high)}_{int(r.year)}",
            r.incidence_pct, sd, "B",
            f"{r.source}{' (CI reaches 0)' if r.uninformative else ''}")
    for year, ratio in ((2011, 1.903), (2016, 2.035)):
        add(f"inc_fm_ratio_{year}", ratio, SIGMA_INC_RATIO, "B",
            "female:male incidence 15-49, robust to the recency assay's MDRI")

    # --- Tier C: age-stratified, for later waves ---
    band_edges = [(15, 25), (25, 35), (35, 45), (45, 65)]
    for year in PHIA_YEARS:
        for sex in ("f", "m"):
            for lo, hi in band_edges:
                b = tg[(tg.year == year) & (tg.sex == sex)
                       & tg.age_low.between(lo, hi - 5)]
                if not len(b):
                    continue
                w = b.Count / b.Count.sum()
                mean = float((w * b.phia).sum())
                sd_ci = float(np.sqrt((w ** 2 * ((b.ub - b.lb) / 3.92) ** 2).sum()))
                add(f"prev_{sex}_{lo}_{hi}_{year}", mean,
                    np.hypot(sd_ci, SIGMA_DISC_PREV), "C",
                    "age-stratified prevalence")
    return obs, pd.DataFrame(prov)


# --- Simulator -----------------------------------------------------------------

def point_key(row):
    """Cache key from the parameter VALUES, not the row index.

    Keying on index would be wrong and silently so: the Latin hypercube design
    changes with n_samples, so point 3 of a 1000-point design is a different
    parameter vector from point 3 of a 12-point smoke test. A value hash makes
    the cache correct across designs, and makes resume work even if HM proposes
    the same point again in a later wave.
    """
    vec = ",".join(f"{float(row[k]):.10g}" for k in sorted(BOUNDS))
    return hashlib.sha1(vec.encode()).hexdigest()[:16]


def _one_point(idx, row):
    """Run one parameter point; cache so an interrupted wave costs nothing."""
    path = SIM_DIR / f"point_{point_key(row)}.parquet"
    if path.exists():
        return

    beta = float(np.exp(row["log_beta_m2f"]))
    rel_f2m = float(np.exp(row["log_rel_beta_f2m"]))
    s_f_young = float(np.exp(row["log_s_f_young"]))
    sd_mult = float(np.exp(row["log_age_gap_sd_mult"]))
    shift = float(row["age_gap_shift"])

    hiv_pars = dict(beta_m2f=beta, rel_beta_f2m=rel_f2m,
                    rel_init_prev=FIXED["rel_init_prev"],
                    rel_sus_age=[(15, 25, 'f', s_f_young),
                                 (25, 50, 'f', 1.0), (15, 50, 'm', 1.0)])
    network_pars = dict(
        prop_f0=float(row["prop_f0"]), prop_m0=float(row["prop_m0"]),
        age_diff_pars={g: [(max(m + shift, 1.0), max(s * sd_mult, 0.2))
                           for m, s in v] for g, v in AGE_DIFF_BASE.items()},
        **{k: v * FIXED["conc_mult"] for k, v in CONC_BASE.items()})

    t0 = time.perf_counter()
    sim = make_sim(seed=SIM_SEED, stop=STOP, verbose=-1, hiv_pars=hiv_pars,
                   network_pars=network_pars,
                   analyzers=[PopByAgeSex(), Cascade()])
    sim.pars.n_agents = N_AGENTS
    sim.run()
    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns
            if c == "timevec" or c.startswith(("popagesex.", "hiv.new_deaths"))]
    out = df[keep].copy()
    out["point"], out["runtime_s"] = idx, time.perf_counter() - t0
    out["point_key"] = point_key(row)
    out.to_parquet(path, index=False)


def _prev(g, sex, lo, hi, year):
    sub = g[g.timevec == year]
    if not len(sub):
        return np.nan
    inf = alv = 0.0
    for b in range(lo, hi, 5):
        ic, ac = (f"popagesex.n_infected_{sex}_{b}_{b+5}",
                  f"popagesex.n_alive_{sex}_{b}_{b+5}")
        if ic not in g.columns:
            return np.nan
        inf += float(sub[ic].iloc[0]); alv += float(sub[ac].iloc[0])
    return inf / alv if alv > 0 else np.nan


def _inc(g, sex, lo, hi, year):
    sub = g[g.timevec == year]
    if not len(sub):
        return np.nan
    lo, hi = (lo // 5) * 5, -(-hi // 5) * 5
    new = alv = inf = 0.0
    for b in range(lo, hi, 5):
        nc = f"popagesex.new_infections_{sex}_{b}_{b+5}"
        if nc not in g.columns:
            return np.nan
        new += float(sub[nc].iloc[0])
        alv += float(sub[f"popagesex.n_alive_{sex}_{b}_{b+5}"].iloc[0])
        inf += float(sub[f"popagesex.n_infected_{sex}_{b}_{b+5}"].iloc[0])
    susc = alv - inf
    return 100.0 * new / susc if susc > 0 else np.nan


def summarise_point(g):
    """One simulation -> one row of named outputs matching the observation keys."""
    r = {}
    for year in PHIA_YEARS:
        for sex in ("f", "m"):
            r[f"prev_15_49_{sex}_{year}"] = _prev(g, sex, 15, 50, year)
        f, m = r[f"prev_15_49_f_{year}"], r[f"prev_15_49_m_{year}"]
        r[f"prev_fm_ratio_{year}"] = f / m if m else np.nan
        yo_num, yo_den = _prev(g, "f", 15, 25, year), _prev(g, "f", 35, 45, year)
        r[f"prev_young_old_f_{year}"] = yo_num / yo_den if yo_den else np.nan
        for sex in ("f", "m"):
            for lo, hi in ((15, 25), (25, 35), (35, 45), (45, 65)):
                r[f"prev_{sex}_{lo}_{hi}_{year}"] = _prev(g, sex, lo, hi, year)
    r[WAVE1_FEATURE] = float(np.nanmean(
        [r[f"prev_15_49_{s}_{y}"] for y in PHIA_YEARS for s in ("f", "m")]))

    inc = pd.read_csv(REPO_ROOT / "calibration_data" / "incidence_by_age_sex.csv")
    for _, t in inc.iterrows():
        r[f"inc_{t.sex}_{int(t.age_low)}_{int(t.age_high)}_{int(t.year)}"] = _inc(
            g, t.sex, int(t.age_low), int(t.age_high), int(t.year))
    for year in (2011, 2016):
        keys = [k for k in r if k.startswith("inc_") and k.endswith(f"_{year}")]
        fv = np.nanmean([r[k] for k in keys if k.startswith("inc_f_")] or [np.nan])
        mv = np.nanmean([r[k] for k in keys if k.startswith("inc_m_")] or [np.nan])
        r[f"inc_fm_ratio_{year}"] = fv / mv if mv and np.isfinite(mv) else np.nan

    d = g.groupby("timevec")["hiv.new_deaths"].mean()
    r["peak_deaths"] = float(d.max()) if len(d) else np.nan
    return r


def make_simulator(n_workers):
    """The function HM calls: DataFrame of samples in, DataFrame of outputs out."""
    def simulator(samples: pd.DataFrame) -> pd.DataFrame:
        samples = samples.reset_index(drop=True)
        keys = [point_key(row) for _, row in samples.iterrows()]
        todo = [(i, dict(row)) for i, row in samples.iterrows()
                if not (SIM_DIR / f"point_{keys[i]}.parquet").exists()]
        if todo:
            print(f"  running {len(todo)}/{len(samples)} points "
                  f"({len(samples) - len(todo)} cached)", flush=True)
            sc.parallelize(_one_point, iterarg=todo, ncpus=n_workers)
        rows = [summarise_point(pd.read_parquet(SIM_DIR / f"point_{k}.parquet"))
                for k in keys]
        return pd.DataFrame(rows)
    return simulator


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_samples", type=int, default=N_SAMPLES)
    p.add_argument("--n_workers", type=int, default=8)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--analyse_only", action="store_true")
    args = p.parse_args()

    obs, prov = build_observations()
    prov.to_csv(OUT_DIR / "observations.csv", index=False)
    print(f"{len(obs)} observations registered "
          f"(tier A {sum(prov.tier == 'A')}, B {sum(prov.tier == 'B')}, "
          f"C {sum(prov.tier == 'C')})")
    print(f"wave 1 emulates: {WAVE1_FEATURE} = "
          f"{obs[WAVE1_FEATURE][0]:.4f} +/- {obs[WAVE1_FEATURE][1]:.4f}")
    print(f"threshold {THRESHOLD}, emulator {DEFAULT_EMULATOR}, "
          f"{args.n_samples} points, N={N_AGENTS:,}\n")

    engine = hm.HistoryMatching(
        function=make_simulator(args.n_workers),
        bounds=BOUNDS,
        observations=obs,
        emulator_type=DEFAULT_EMULATOR,
        feature_selection=hm.ManualFeatureSelection([WAVE1_FEATURE]),
        n_samples=args.n_samples,
        max_iterations=1,
        implausibility_threshold=THRESHOLD,
        random_seed=RANDOM_SEED,
        output_dir=str(HM_DIR),
        run_name="wave1",
    )

    if not args.analyse_only:
        t0 = sc.tic()
        engine.run(resume=args.resume)
        sc.toc(t0, label="wave 1")

    print("\n=== Emulator quality ===")
    try:
        engine.print_emulator_quality_metrics()
    except Exception as e:
        print(f"  (unavailable: {e})")

    summary = engine.nroy_summary()
    print(f"\n{summary}")
    (OUT_DIR / "nroy_summary.txt").write_text(str(summary), encoding="utf-8")
    print(f"\ndiagnostics -> {HM_DIR / 'wave1'}")


if __name__ == "__main__":
    main()
