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
import starsim as ss  # noqa: E402
import stisim as sti  # noqa: E402
from run_sims import make_sim  # noqa: E402
from analyzers import PopByAgeSex, Cascade  # noqa: E402
from standard_figures import (load_targets, targets_15_49,  # noqa: E402
                              load_inc_targets, _inc_series)

MODEL_TAG = "model-v1.4"  # v1.3 + hiv_cd4_floor.HIVCD4Floor; see 025's SUMMARY
N_AGENTS = 20_000       # up from 024: wave-2 features are age-stratified (020)
STOP = 2026
SIM_SEED = 1
N_SAMPLES = 1000
THRESHOLD = 4.0          # raised from 3.0: known structural misspecification
RANDOM_SEED = 20260907
PHIA_YEARS = (2007, 2011, 2016)

# Loaded once. summarise_point() runs in a worker process per point, so reading
# the target file there would re-read it 1000 times.
_TG = load_targets()
_INC = load_inc_targets()

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
SIGMA_DISC_PREV = 0.01
# The deaths down-weighting. HM has no weights, only sigmas. The model reaches
# 7,069 against UNAIDS' 11,000, so 2000 puts that gap at ~2 sigma; a survey-like
# sigma would put it at 8 sigma and rule out the entire box.
UNAIDS_PEAK, SIGMA_DEATHS = 11_000.0, 2000.0
SIGMA_INC_DISC_REL = 0.10    # relative model discrepancy on incidence, researcher
                             # decision 2026-09-07: a modelled incidence and a
                             # recency-assay/cohort estimate are not the same
                             # quantity, and prevalence already carries an
                             # allowance -- without one, incidence would bite
                             # harder than prevalence.
SIGMA_RATIO_REL = 0.15       # prevalence ratios, relative
SIGMA_INC_RATIO = 0.35       # F:M incidence ratio, on a value near 2.0

WAVE1_FEATURE = "prev_15_49_all_mean"
# Wave 2 emulates two features, each with a specific job -- see README.
WAVE2_FEATURES = ["inc_all_mean", "prev_young_old_f_mean"]


def tier_c_bands(year, tg):
    """Tier C band edges for one survey year, capped at that survey's own coverage.

    The top band must not extend past the data. SHIMS1 (2011) publishes strata
    only to [45:50), while 2007 and 2016 reach [60:65). A fixed (45, 65) band
    therefore compared a target built from PHIA 45-49 against a model average
    over 45-65 -- and male prevalence falls steeply after 50, so the model looked
    5 sigma low on a band it actually fits (z -4.99 -> +0.30 once matched).

    That was the only real defect wave 1 reported, and it was in this function,
    not in the model. Deriving the cap from the target file makes it impossible
    to reintroduce when a survey with different coverage is added.
    """
    top = int(tg[tg.year == year].age_low.max()) + 5
    edges = [(15, 25), (25, 35), (35, 45), (45, top)]
    return [(lo, hi) for lo, hi in edges if hi > lo]


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
    yos = []
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
        yos.append(yo)

    # CORRECTION (025): read `sigma` from the target file rather than
    # recomputing it. incidence_construction.py owns that number now, and 024's
    # version silently discarded it -- along with the model_age_basis column,
    # which is what stops SHIMS1's 18-49 estimate being compared against a model
    # 15-49 one. That mismatch flipped the sign of the 2011 male residual.
    add("prev_young_old_f_mean", float(np.mean(yos)),
        SIGMA_RATIO_REL * float(np.mean(yos)), "B",
        "female 15-24 : 35-44 prevalence, mean over PHIA years -- wave 2 feature")

    inc = load_inc_targets()
    for _, r in inc.iterrows():
        add(f"inc_{r.sex}_{int(r.age_low)}_{int(r.age_high)}_{int(r.year)}",
            r.incidence_pct,
            float(np.hypot(r.sigma, SIGMA_INC_DISC_REL * r.incidence_pct)), "B",
            f"{r.source} | basis {r.model_age_basis} | "
            f"CI sigma {r.sigma:.3f} + {SIGMA_INC_DISC_REL:.0%} discrepancy")
    add("inc_all_mean", float(inc.incidence_pct.mean()),
        float(np.hypot(inc.sigma.mean(),
                       SIGMA_INC_DISC_REL * inc.incidence_pct.mean())), "B",
        "mean of the four fitted sex-specific incidence aggregates -- wave 2 feature")
    for year, ratio in ((2011, 1.903), (2016, 2.035)):
        add(f"inc_fm_ratio_{year}", ratio, SIGMA_INC_RATIO, "B",
            "female:male incidence 15-49, robust to the recency assay's MDRI")

    # --- Tier C: age-stratified, for later waves ---
    for year in PHIA_YEARS:
        for sex in ("f", "m"):
            for lo, hi in tier_c_bands(year, tg):
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
    # CORRECTION (025): the simulation config must be in the key too. 024 hashed
    # only the parameter values, so raising N from 10,000 to 20,000 would have
    # returned stale 10,000-agent results for any repeated point -- silently, and
    # it would have corrupted the wave. Found while planning the N change.
    vec += f"|N={N_AGENTS}|stop={STOP}|seed={SIM_SEED}"
    # CORRECTION (027): the MODEL VERSION must be in the key too. 025 obs 6 --
    # the key hashed parameters and sim config but nothing about the model, so
    # after the model-v1.4 CD4 floor landed, a --resume would have silently
    # served v1.3 results. Same class of bug as the N_AGENTS omission above, one
    # level further out. Folder separation alone is not a defence: it protects
    # across experiments, not across a model change inside one.
    vec += f"|model={MODEL_TAG}|stisim={sti.__version__}|starsim={ss.__version__}"
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


def _basis_weight(basis):
    """First-band person-time weight implied by a target's model_age_basis.

    "18_50_prorate_first_band_0.4" -> 0.4, because 18-19 is 2/5 of the 15-19
    band. "15_50" -> 1.0. Parsed rather than hardcoded so the target file stays
    the single source of truth.
    """
    if basis is None or "prorate" not in str(basis):
        return 1.0
    return float(str(basis).rsplit("_", 1)[-1])


def _inc(g, sex, lo, hi, year, prorate_first=1.0):
    """Annual incidence % for one sex/age range/year, over SUSCEPTIBLES.

    CORRECTION (025): `prorate_first` down-weights the first 5-year band's
    person-time, which is how an 18-49 estimate is expressed on 5-year bands.
    024 floored 18 to 15 and compared a model 15-49 estimate against SHIMS1's
    18-49 target. Male incidence is near zero at 15-19, so including it diluted
    the average and flipped the sign of the 2011 male residual (-0.67 -> +0.26).
    """
    sub = g[g.timevec == year]
    if not len(sub):
        return np.nan
    lo, hi = (lo // 5) * 5, -(-hi // 5) * 5
    new = alv = inf = 0.0
    for i, b in enumerate(range(lo, hi, 5)):
        w = prorate_first if i == 0 else 1.0
        nc = f"popagesex.new_infections_{sex}_{b}_{b+5}"
        if nc not in g.columns:
            return np.nan
        new += w * float(sub[nc].iloc[0])
        alv += w * float(sub[f"popagesex.n_alive_{sex}_{b}_{b+5}"].iloc[0])
        inf += w * float(sub[f"popagesex.n_infected_{sex}_{b}_{b+5}"].iloc[0])
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
            # Same per-year cap as the observation side, or the two disagree
            # about what "the top band" means -- see tier_c_bands().
            for lo, hi in tier_c_bands(year, _TG):
                r[f"prev_{sex}_{lo}_{hi}_{year}"] = _prev(g, sex, lo, hi, year)
    r[WAVE1_FEATURE] = float(np.nanmean(
        [r[f"prev_15_49_{s}_{y}"] for y in PHIA_YEARS for s in ("f", "m")]))

    # Each row is computed on ITS OWN published age basis, read from the target
    # file's model_age_basis column -- not on a single convention for all rows.
    for _, t in _INC.iterrows():
        r[f"inc_{t.sex}_{int(t.age_low)}_{int(t.age_high)}_{int(t.year)}"] = _inc(
            g, t.sex, int(t.age_low), int(t.age_high), int(t.year),
            prorate_first=_basis_weight(t.model_age_basis))
    r["inc_all_mean"] = float(np.nanmean(
        [r[f"inc_{t.sex}_{int(t.age_low)}_{int(t.age_high)}_{int(t.year)}"]
         for _, t in _INC.iterrows()]))
    r["prev_young_old_f_mean"] = float(np.nanmean(
        [r[f"prev_young_old_f_{y}"] for y in PHIA_YEARS]))
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


# --- Truth point for the re-identification pre-flight -------------------------
# 024's design row 868: the draw that fit 48/48 targets, so it is both plausible
# and already characterised. Stored in natural units; the box is log-space.
TRUTH_NATURAL = dict(beta_m2f=0.014161835611057863,
                     rel_beta_f2m=0.3097261030208422,
                     s_f_young=2.224000314771824,
                     age_gap_shift=0.0699801603950227,
                     age_gap_sd_mult=1.5001144712607821,
                     prop_f0=0.608751696061895,
                     prop_m0=0.4742150061435999)
TRUTH_SEED = 7          # NOT SIM_SEED: the synthetic target must carry one
                        # replicate's noise, not be the model's exact own output


def truth_row():
    """TRUTH_NATURAL in the box's coordinates."""
    return {"log_beta_m2f": np.log(TRUTH_NATURAL["beta_m2f"]),
            "log_rel_beta_f2m": np.log(TRUTH_NATURAL["rel_beta_f2m"]),
            "log_s_f_young": np.log(TRUTH_NATURAL["s_f_young"]),
            "age_gap_shift": TRUTH_NATURAL["age_gap_shift"],
            "log_age_gap_sd_mult": np.log(TRUTH_NATURAL["age_gap_sd_mult"]),
            "prop_f0": TRUTH_NATURAL["prop_f0"],
            "prop_m0": TRUTH_NATURAL["prop_m0"]}


def synthetic_observations(obs, n_workers):
    """Replace every target MEAN with the truth point's own output, keeping sigma.

    Sigmas are kept from the real observations deliberately. The question the
    pre-flight asks is "can this pipeline, with these uncertainties, recover a
    known point" -- shrinking sigma to make recovery easy would answer a
    different and useless question.

    Structural error cancels by construction here, since the data comes from the
    model. This tests the PIPELINE (prior, features, emulator, implausibility),
    not whether the model is right.
    """
    global SIM_SEED
    row = truth_row()
    saved, SIM_SEED = SIM_SEED, TRUTH_SEED
    try:
        path = SIM_DIR / f"point_{point_key(row)}.parquet"
        if not path.exists():
            _one_point(0, row)
        truth_out = summarise_point(pd.read_parquet(path))
    finally:
        SIM_SEED = saved

    syn, missing = {}, []
    for k, (mean, sd) in obs.items():
        v = truth_out.get(k)
        if v is None or not np.isfinite(v):
            missing.append(k)
            continue
        syn[k] = (float(v), float(sd))
    if missing:
        tail = "..." if len(missing) > 4 else ""
        print(f"  ({len(missing)} targets dropped, not finite at truth: "
              f"{missing[:4]}{tail})")
    return syn, truth_out


def run_waves(obs, features, n_workers, n_samples, run_name, resume=False,
              max_iterations=3):
    """One HistoryMatching job over `max_iterations` NROY-refinement waves.

    NOTE on what a "wave" is here. ManualFeatureSelection(list) emulates EVERY
    listed feature in EVERY wave -- it is not a per-wave schedule. So all three
    features are applied simultaneously, and successive waves re-draw the design
    inside the surviving NROY rather than layering a new check.

    That is a deliberate choice, not a workaround. The question 024 left open is
    whether the model can satisfy prevalence AND incidence AND the female
    young:old ratio *jointly* -- 024's headline was a joint 48/48 fit -- and
    simultaneous emulation answers that directly. Sequential layering would
    answer "does each check cut further", which is a weaker question here.

    It does deviate from the history-matching guidance of 1-2 features per wave,
    on the grounds that more features dilute the signal. Mitigation: the package
    writes per-feature emulator metrics and z-scores, so attribution of which
    feature did the cutting survives. If an emulator comes back with R^2 < 0.8 it
    should be dropped rather than carried.
    """
    engine = hm.HistoryMatching(
        function=make_simulator(n_workers),
        bounds=BOUNDS,
        observations=obs,
        emulator_type=DEFAULT_EMULATOR,
        feature_selection=hm.ManualFeatureSelection(features),
        n_samples=n_samples,
        max_iterations=max_iterations,
        implausibility_threshold=THRESHOLD,
        random_seed=RANDOM_SEED,
        output_dir=str(HM_DIR),
        run_name=run_name,
    )
    t0 = sc.tic()
    engine.run(resume=resume)
    sc.toc(t0, label=run_name)
    return engine


def report(engine, label, outstem):
    print(f"\n=== {label}: emulator quality ===")
    try:
        engine.print_emulator_quality_metrics()
    except Exception as e:
        print(f"  (unavailable: {e})")
    summary = engine.nroy_summary()
    print(f"\n{summary}")
    (OUT_DIR / f"{outstem}_nroy_summary.txt").write_text(str(summary),
                                                         encoding="utf-8")
    return summary


def check_recovery(outstem):
    """Is the truth inside the NROY, and did the beta/s_f ridge narrow?"""
    truth = truth_row()
    nroy = None
    for cand in sorted((HM_DIR / outstem).rglob("nroy_samples.csv")):
        nroy = pd.read_csv(cand)          # last = final wave
    if nroy is None or not len(nroy):
        print("  no NROY samples found -- cannot assess recovery")
        return None

    rows = []
    for k in sorted(BOUNDS):
        lo, hi = np.percentile(nroy[k], [2.5, 97.5])
        width_prior = BOUNDS[k][1] - BOUNDS[k][0]
        rows.append(dict(parameter=k, truth=truth[k], nroy_lo=lo, nroy_hi=hi,
                         inside=bool(lo <= truth[k] <= hi),
                         frac_of_prior_width=(hi - lo) / width_prior))
    rec = pd.DataFrame(rows)

    # The degeneracy 024's PC1 identified: beta_m2f and s_f_young enter as a
    # product. Recovery of the PRODUCT with the RATIO still free is the
    # pre-registered failure mode, so measure both explicitly.
    prod = nroy["log_beta_m2f"] + nroy["log_s_f_young"]
    ratio = nroy["log_beta_m2f"] - nroy["log_s_f_young"]
    t_prod = truth["log_beta_m2f"] + truth["log_s_f_young"]
    t_ratio = truth["log_beta_m2f"] - truth["log_s_f_young"]
    deg = pd.DataFrame([
        dict(quantity="log(beta_m2f * s_f_young)", truth=t_prod,
             lo=np.percentile(prod, 2.5), hi=np.percentile(prod, 97.5),
             sd=float(prod.std())),
        dict(quantity="log(beta_m2f / s_f_young)", truth=t_ratio,
             lo=np.percentile(ratio, 2.5), hi=np.percentile(ratio, 97.5),
             sd=float(ratio.std())),
    ])
    deg["inside"] = (deg.lo <= deg.truth) & (deg.truth <= deg.hi)
    deg["corr_beta_sfy"] = float(nroy["log_beta_m2f"].corr(nroy["log_s_f_young"]))

    rec.to_csv(OUT_DIR / f"{outstem}_recovery.csv", index=False)
    deg.to_csv(OUT_DIR / f"{outstem}_degeneracy.csv", index=False)
    print(f"\n=== {outstem}: parameter recovery ===")
    print(rec.round(4).to_string(index=False))
    print(f"\n=== {outstem}: the beta_m2f / s_f_young degeneracy ===")
    print(deg.round(4).to_string(index=False))
    print(f"\ntruth inside NROY on {int(rec.inside.sum())}/{len(rec)} parameters; "
          f"product recovered: {bool(deg.iloc[0].inside)}; "
          f"ratio recovered: {bool(deg.iloc[1].inside)}")
    return rec, deg


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_samples", type=int, default=N_SAMPLES)
    p.add_argument("--n_workers", type=int, default=8)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--max_waves", type=int, default=3,
                   help="NROY-refinement waves; the engine stops early if the "
                        "NROY falls below 20 samples, which IS the collapse signal")
    p.add_argument("--part", choices=["preflight", "wave2", "both"],
                   default="both",
                   help="preflight GATES wave2; 'both' runs preflight then wave2")
    args = p.parse_args()

    obs, prov = build_observations()
    prov.to_csv(OUT_DIR / "observations.csv", index=False)
    print(f"{len(obs)} observations registered "
          f"(tier A {sum(prov.tier == 'A')}, B {sum(prov.tier == 'B')}, "
          f"C {sum(prov.tier == 'C')})")
    # CORRECTION (027): 025 printed "wave 1 emulates: X / wave 2 emulates: Y, Z"
    # while passing all three as one list to ManualFeatureSelection, so the log
    # claimed a per-wave schedule the code does not implement and contradicted
    # its own README. Print what actually happens.
    feats = [WAVE1_FEATURE] + WAVE2_FEATURES
    print(f"emulated every wave (simultaneously): {', '.join(feats)}")
    for f in feats:
        print(f"    {f:24} = {obs[f][0]:.4f} +/- {obs[f][1]:.4f}")
    print(f"threshold {THRESHOLD}, emulator {DEFAULT_EMULATOR}, "
          f"{args.n_samples} points/wave, N={N_AGENTS:,}, "
          f"sigma_disc_prev={SIGMA_DISC_PREV}, model={MODEL_TAG}, "
          f"starsim={ss.__version__}, stisim={sti.__version__}\n")

    features = [WAVE1_FEATURE] + WAVE2_FEATURES

    if args.part in ("preflight", "both"):
        print("=" * 72)
        print("PART 1 -- re-identification pre-flight (synthetic data)")
        print("=" * 72)
        syn, truth_out = synthetic_observations(obs, args.n_workers)
        pd.Series(truth_out).to_csv(OUT_DIR / "preflight_truth_outputs.csv")
        eng = run_waves(syn, features, args.n_workers, args.n_samples,
                        "preflight", resume=args.resume,
                        max_iterations=args.max_waves)
        report(eng, "preflight", "preflight")
        check_recovery("preflight")
        print("\n>>> PART 1 complete. Inspect recovery before trusting wave 2's")
        print(">>> marginals in the beta_m2f / s_f_young directions.\n")

    if args.part in ("wave2", "both"):
        print("=" * 72)
        print("PART 2 -- wave 2 on real data")
        print("=" * 72)
        eng = run_waves(obs, features, args.n_workers, args.n_samples,
                        "wave2", resume=args.resume,
                        max_iterations=args.max_waves)
        report(eng, "wave 2", "wave2")

    print(f"\ndiagnostics -> {HM_DIR}")


if __name__ == "__main__":
    main()
