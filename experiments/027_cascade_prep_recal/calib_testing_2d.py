"""Exp 027 step 2d: age-structured testing and treatment, testing recalibrated.

Researcher's instructions, DTC_check/context_dtc.md, 2026-11-19. Fixed inputs
(all opt-in in make_sim; none is calibrated here):

    ART coverage    forced ART among PLHIV on the testing cohorts 0-24, 25-34,
                    35-49, 50+, from the SHIMS 90-90-90 file; 0-14 take the 15-24
                    values (art_coverage='cohorts';
                    art_coverage_cohort_construction.py)
    ANC testing     logistic per pregnancy: starts 1992.5, midpoint 2003, 95%
                    of the 90% maximum by 2010 (anc_logistic=True)
    PMTCT           as 2c, 0 in 2003 -> 90% 2010 (pmtct.py)
    VLS among ART   as 2c, logistic through SHIMS2 and SHIMS3 ('phia_logistic')
    VMMC            SHIMS2 Table 15.5.A totals in the SHIMS2 row, moved to 2017
                    (vmmc_shims2=True)

plus PrEP, shims2_year=2017, row 868's transmission parameters. Free (fixed
testing rates, goal 5), 10 parameters:

    gp_start_m/f        year general-population testing begins, men / women
    rel_<s>_<cohort>    testing level per sex x cohort (x 0.5/yr in 2020);
                        children are in the 0-24 cohort

Loss, as 2c but by age:
  - SHIMS2 first and second 90 at 2017.0, by sex, for 15-24, 25-34, 35-49 and
    15+ (SHIMS2 has no 50+ band; 15+ carries it). sigma from the Wilson CI on
    the table's unweighted n.
  - ART-delivery penalty (2b): share of the forced cohort ART not delivered,
    worst of 2007.5/2011.5/2014.5, penalised above 5%.
  - Kasturi guide (2c), 15-49 at 2005.5/2008.5/2011.5/2014.5: men:women ratio of
    % diagnosed (sigma 0.10) and each sex's level (sigma 0.20). A guide, not a
    target (researcher, 2026-10-28 and 2026-11-19).
SHIMS3 (2021) is held out and only reported.

Search: optuna TPE (multivariate), BATCH trials at a time evaluated in
parallel on SEEDS_SEARCH seeds, the 2c chosen point enqueued as a start; then
the best TOP_FINAL trials re-run on SEEDS_FINAL fresh seeds and the best of
those chosen. Outputs: outputs/calib2d_*.csv, calib2d_chosen.json.
"""

import os
import sys
from pathlib import Path

import numpy as np
import optuna
import pandas as pd
import sciris as sc
import starsim as ss

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
os.chdir(REPO)

from run_sims import make_sim  # noqa: E402
from interventions import AGE_COHORTS  # noqa: E402
from tune_testing import BEST, build_pars, SHIMS  # noqa: E402

OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

COHORT_KEYS = [f"{lo}" for lo, _ in AGE_COHORTS]          # "0", "25", "35", "50"
BOUNDS = {"gp_start_m": (1990, 2012), "gp_start_f": (1990, 2012),
          **{f"rel_{s}_{c}": (0.05, 3.0) for s in "mf" for c in COHORT_KEYS}}
NAMES = list(BOUNDS)

SHIMS_BANDS = {"15-24": (15, 25), "25-34": (25, 35), "35-49": (35, 50), "15+": (15, 100)}
SHORT_TIMES = (2007.5, 2011.5, 2014.5)
KAS_TIMES = (2005.5, 2008.5, 2011.5, 2014.5)
SNAPS = sorted(set(SHORT_TIMES + KAS_TIMES + (2017.0, 2021.0)))
ART_BINS = [(0, 25), (25, 35), (35, 50), (50, 100)]
KAS_BINS = [(a, a + 5) for a in range(15, 50, 5)]
SHORT_TOL, SHORT_SCALE = 0.05, 0.02
KAS_RATIO_SIGMA, KAS_LEVEL_SIGMA = 0.10, 0.20
N_TRIALS, BATCH, SEEDS_SEARCH = 360, 9, 2
TOP_FINAL, SEEDS_FINAL = 10, 8
NCPUS = 9

MODEL_KW = dict(prep=True, shims2_year=2017, pmtct=True, art_vls_coverage="phia_logistic",
                art_coverage="cohorts", vmmc_shims2=True)

# --- SHIMS2 targets: first and second 90 by band ------------------------------
STAGE = {"Diagnosed": "dx", "Diagnosed_OnART": "art_dx", "OnART_VLS": "vls"}
_s = SHIMS[(SHIMS.Year == 2017) & SHIMS["Age group"].isin(SHIMS_BANDS)].copy()
_s["stage"] = _s.Metric.map(STAGE)
_p, _n, _z = _s.Percentage / 100, _s.Number, 1.96
_s["p"] = _p
_s["sigma"] = _z * np.sqrt(_p * (1 - _p) / _n + _z**2 / (4 * _n**2)) / (1 + _z**2 / _n) / _z
TARGETS = {(r.stage, r.sex, r["Age group"]): (r.p, r.sigma) for _, r in _s.iterrows()}
SCORED = ("dx", "art_dx")

# --- Kasturi % aware by 5-year band (guide) -----------------------------------
_k = pd.read_csv(REPO / "DTC_check" / "SHIMS_Datasets" /
                 "aware_and_art_no_interactions_modeled_output_1995_2021-From_Kasturi.csv",
                 index_col=0)
_k["sex"] = _k.gender.str.lower().map({"female": "f", "male": "m"})
_k["lo"] = _k.agecat.str.split("-").str[0].astype(int)
KAS = {(int(r.year), r.sex, r.lo): r.aware_extrapolated for _, r in _k.iterrows()}


def art_targets(t):
    """Forced cohort ART among PLHIV at time t, per (sex, bin), SHIMS2 at 2017."""
    a = pd.read_csv(REPO / "data" / "art_coverage_cohorts.csv")
    a.loc[a.Year == 2016, "Year"] = 2017
    a["sex"] = a.Gender.map({0: "f", 1: "m"})
    out = {}
    for (sex, ab), s in a.groupby(["sex", "AgeBin"]):
        lo, hi = (int(x) for x in ab.strip("[)").split(","))
        s = s.sort_values("Year")
        out[(sex, (lo, hi))] = float(np.interp(t, s.Year, s.p_art))
    return out


ART_TARGETS = {t: art_targets(t) for t in SHORT_TIMES}


class Snap(ss.Analyzer):
    """PLHIV / diagnosed / on-ART / suppressed counts at exact times."""

    def init_pre(self, sim):
        super().init_pre(sim)
        self.rows = []

    def step(self):
        y = self.t.now("year")
        hit = [s for s in SNAPS if abs(y - s) < 1e-6]
        if not hit:
            return
        ppl, hiv = self.sim.people, self.sim.diseases.hiv
        groups = {**SHIMS_BANDS, "15-49": (15, 50),
                  **{f"art_{lo}_{hi}": (lo, hi) for lo, hi in ART_BINS},
                  **{f"k_{lo}_{hi}": (lo, hi) for lo, hi in KAS_BINS}}
        for sex, sb in (("f", ppl.female), ("m", ppl.male)):
            for name, (lo, hi) in groups.items():
                inb = ppl.alive & sb & (ppl.age >= lo) & (ppl.age < hi)
                self.rows.append(dict(year=hit[0], sex=sex, age=name,
                                      plhiv=(inb & hiv.infected).count(),
                                      dx=(inb & hiv.diagnosed).count(),
                                      art=(inb & hiv.on_art).count(),
                                      vls=(inb & hiv.on_effective_art).count()))


def testing_pars(x):
    """Parameter dict (names as BOUNDS) -> make_sim testing_pars."""
    return dict(gp_start_m=x["gp_start_m"], gp_start_f=x["gp_start_f"],
                rel_test_age={s: [x[f"rel_{s}_{c}"] for c in COHORT_KEYS] for s in "mf"},
                anc=True, anc_logistic=True)


def run_point(pid, x, seed):
    hiv_pars, network_pars = build_pars(BEST)
    sim = make_sim(seed=seed, start=1985, stop=2022, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars,
                   testing_pars=testing_pars(x), analyzers=[Snap(name="snap")],
                   **MODEL_KW)
    sim.run()
    d = pd.DataFrame(sim.analyzers.snap.rows)
    d["pid"], d["seed"] = pid, seed
    return d


def score(raw):
    """Pool seeds per point -> loss, its parts, ART shortfall, Kasturi comparison."""
    g = raw.groupby(["pid", "year", "sex", "age"])[["plhiv", "dx", "art", "vls"]].sum().reset_index()
    g["p_dx"], g["p_art_dx"], g["p_vls"] = g.dx / g.plhiv, g.art / g.dx, g.vls / g.art
    rows = []
    for pid, gp in g.groupby("pid"):
        r = dict(pid=pid, loss_shims=0.0, loss_kas=0.0)
        at = lambda y, sex, age: gp[(gp.year == y) & (gp.sex == sex) & (gp.age == age)].iloc[0]  # noqa: E731
        for sex in ("f", "m"):
            for band in list(SHIMS_BANDS) + ["15-49"]:
                a17, a21 = at(2017.0, sex, band), at(2021.0, sex, band)
                for stage, col in (("dx", "p_dx"), ("art_dx", "p_art_dx"), ("vls", "p_vls")):
                    r[f"{stage}17_{sex}_{band}"] = a17[col]
                    r[f"{stage}21_{sex}_{band}"] = a21[col]
                    if stage in SCORED and band in SHIMS_BANDS:
                        p, sig = TARGETS[(stage, sex, band)]
                        r["loss_shims"] += ((a17[col] - p) / sig) ** 2
            for t in KAS_TIMES:
                bands = [at(t, sex, f"k_{lo}_{hi}") for lo, hi in KAS_BINS]
                w = np.array([b.plhiv for b in bands], dtype=float)
                kas = float(sum(KAS[(int(t), sex, lo)] * wi
                                for (lo, _), wi in zip(KAS_BINS, w)) / w.sum())
                a = at(t, sex, "15-49")
                r[f"dx{int(t)}_{sex}"], r[f"kas{int(t)}_{sex}"] = a.p_dx, kas
                r[f"artdx{int(t)}_{sex}"] = a.p_art_dx
                r["loss_kas"] += ((a.p_dx - kas) / KAS_LEVEL_SIGMA) ** 2
        for t in KAS_TIMES:
            y = int(t)
            ratio = r[f"dx{y}_m"] / r[f"dx{y}_f"]
            r[f"ratio{y}"], r[f"kasratio{y}"] = ratio, r[f"kas{y}_m"] / r[f"kas{y}_f"]
            r["loss_kas"] += ((ratio - r[f"kasratio{y}"]) / KAS_RATIO_SIGMA) ** 2
        short = []
        for t, tg in ART_TARGETS.items():
            want = miss = 0.0
            for (sex, (lo, hi)), target in tg.items():
                row = at(t, sex, f"art_{lo}_{hi}")
                want += target * row.plhiv
                miss += max(0.0, target * row.plhiv - row.art)
            short.append(miss / want if want else 0.0)
        r["art_short"] = max(short)
        r["loss_short"] = (max(0.0, r["art_short"] - SHORT_TOL) / SHORT_SCALE) ** 2
        r["loss"] = r["loss_shims"] + r["loss_short"] + r["loss_kas"]
        rows.append(r)
    return pd.DataFrame(rows)


def evaluate(points, seeds):
    """points: list of (pid, x). Runs every point x seed in parallel, scores."""
    args = [(pid, x, s) for pid, x in points for s in seeds]
    raw = pd.concat(sc.parallelize(run_point, iterarg=args, ncpus=NCPUS))
    return score(raw)


def start_point():
    """2c's chosen point, spread evenly over the cohorts."""
    c = sc.loadjson(OUT / "calib2c_chosen.json")
    x = {"gp_start_m": c["gp_start_m"], "gp_start_f": c["gp_start_f"]}
    x.update({f"rel_{s}_{k}": c[f"rel_test_{s}"] for s in "mf" for k in COHORT_KEYS})
    return x


if __name__ == "__main__":
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=20261119, multivariate=True,
                                           constant_liar=True, n_startup_trials=45))
    study.enqueue_trial(start_point())
    records = []
    while len(study.trials) < N_TRIALS:
        trials = [study.ask() for _ in range(BATCH)]
        xs = [{k: t.suggest_float(k, *BOUNDS[k]) for k in NAMES} for t in trials]
        s = evaluate([(t.number, x) for t, x in zip(trials, xs)], range(1, SEEDS_SEARCH + 1))
        for t, x in zip(trials, xs):
            r = s[s.pid == t.number].iloc[0]
            study.tell(t, float(r.loss))
            records.append({**x, **r.to_dict()})
        best = min(records, key=lambda d: d["loss"])
        print(f"trials {len(study.trials):4d}  best loss {best['loss']:.1f} "
              f"(shims {best['loss_shims']:.1f}, kas {best['loss_kas']:.1f}, "
              f"short {best['art_short']:.3f})", flush=True)
    search = pd.DataFrame(records)
    search.to_csv(OUT / "calib2d_search.csv", index=False)

    # Final: best trials again on fresh seeds
    top = search.nsmallest(TOP_FINAL, "loss")
    pts = [(20_000 + i, {k: r[k] for k in NAMES}) for i, (_, r) in enumerate(top.iterrows())]
    fin = evaluate(pts, range(101, 101 + SEEDS_FINAL))
    fin = pd.DataFrame([{**x, "pid": pid} for pid, x in pts]).merge(fin, on="pid")
    fin.to_csv(OUT / "calib2d_final.csv", index=False)
    chosen = fin.nsmallest(1, "loss").iloc[0]
    sc.savejson(OUT / "calib2d_chosen.json",
                {**{k: float(chosen[k]) for k in NAMES},
                 "loss": float(chosen.loss), "loss_shims": float(chosen.loss_shims),
                 "loss_kas": float(chosen.loss_kas), "art_short": float(chosen.art_short)})

    with pd.option_context("display.width", 240, "display.float_format", "{:.3f}".format):
        print("\nFinal (8 fresh seeds), best first:")
        print(fin.sort_values("loss")[NAMES + ["loss", "loss_shims", "loss_kas", "art_short"]]
              .to_string(index=False))
        print("\nChosen vs SHIMS2 (2017.0) [target]; 2021 held out:")
        for sex in ("f", "m"):
            for band in SHIMS_BANDS:
                line = f"  {sex} {band:5s}"
                for stage in ("dx", "art_dx", "vls"):
                    tgt = TARGETS.get((stage, sex, band), (np.nan,))[0]
                    line += f"  {stage} {chosen[f'{stage}17_{sex}_{band}']:.3f} [{tgt:.3f}]"
                line += f"   | 2021 dx {chosen[f'dx21_{sex}_{band}']:.3f}"
                print(line)
        print("\nChosen, 15-49 % diagnosed vs Kasturi (guide):")
        for sex in ("f", "m"):
            print(f"  {sex}  " + "  ".join(
                f"{int(t)}: {chosen[f'dx{int(t)}_{sex}']:.2f} (K {chosen[f'kas{int(t)}_{sex}']:.2f})"
                for t in KAS_TIMES))
        print("  m:f ratio  " + "  ".join(
            f"{int(t)}: {chosen[f'ratio{int(t)}']:.2f} (K {chosen[f'kasratio{int(t)}']:.2f})"
            for t in KAS_TIMES))
