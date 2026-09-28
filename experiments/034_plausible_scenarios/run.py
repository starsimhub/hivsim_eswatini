"""Exp 033 -- how far can the cascade actually go? 032's ladder, extended.

031 showed that an arm called "ART coverage to 95%" delivered 0.818 in men
25-34, because awareness is a hard ceiling and nothing warns when a target sits
above it. The fix is not a better label: it is to stop treating the cascade as
a control variable at all. Here every cell is plotted at the viral suppression
it ACHIEVED, so a cell that falls short simply lands further left.

Grid: 5 cascade rungs x 5 PrEP rungs x 10 seeds = 250 sims.

Usage (from the repo ROOT -- run_sims reads data/ by relative path)
  python experiments/034_plausible_scenarios/run.py --smoke
  python experiments/034_plausible_scenarios/run.py --n_workers 110

Structure is 026/031's, deliberately: per-cell parquet so an interrupted job
resumes free, and a fingerprint in the cache key so changing a ramp invalidates
only the cells it affects.
"""

import argparse
import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import sciris as sc

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from run_sims import make_sim                                 # noqa: E402
from analyzers import PopByAgeSex, Cascade, CascadeByAge      # noqa: E402
from interventions import TEST_RATE_M, TEST_RATE_F            # noqa: E402
import scenarios as S                                         # noqa: E402

OUT = HERE / "outputs"
SIM_DIR = OUT / "sims"
for d in (OUT, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

START, STOP = 1985, 2041
SCEN_START, SCEN_REACH = 2026, 2030
N_SEEDS = 10
KEEP_PREFIX = ("popagesex.", "cascade.", "cascadeage.", "hiv.")

# 024 design row 868 -- identical to 026, 029 and 031. Copied, not re-derived.
BEST = dict(beta_m2f=0.014161835611057863, rel_beta_f2m=0.3097261030208422,
            s_f_young=2.224000314771824, age_gap_shift=0.0699801603950227,
            age_gap_sd_mult=1.5001144712607821, prop_f0=0.608751696061895,
            prop_m0=0.4742150061435999)

AGE_DIFF_BASE = {"teens": [(7, 3), (6, 3), (5, 1)],
                 "young": [(8, 3), (7, 3), (5, 2)],
                 "adult": [(8, 3), (7, 3), (5, 2)]}
CONC_BASE = {"f1_conc": 0.15, "f2_conc": 0.25, "m1_conc": 0.15, "m2_conc": 0.5}
REL_INIT_PREV, CONC_MULT = 0.2, 1.0


def build_pars(p):
    """Identical to 025/026/031. Copied so the cells stay comparable to the
    calibration that justified their parameters."""
    hiv_pars = dict(
        beta_m2f=p["beta_m2f"], rel_beta_f2m=p["rel_beta_f2m"],
        rel_init_prev=REL_INIT_PREV,
        rel_sus_age=[(15, 25, "f", p["s_f_young"]),
                     (25, 50, "f", 1.0), (15, 50, "m", 1.0)])
    shift, sd_mult = p["age_gap_shift"], p["age_gap_sd_mult"]
    network_pars = dict(
        prop_f0=p["prop_f0"], prop_m0=p["prop_m0"],
        age_diff_pars={g: [(max(m + shift, 1.0), max(s * sd_mult, 0.2))
                           for m, s in v] for g, v in AGE_DIFF_BASE.items()},
        **{k: v * CONC_MULT for k, v in CONC_BASE.items()})
    return hiv_pars, network_pars


# --- the cascade axis (x) ------------------------------------------------------
# Nested and cumulative. (art_target, vls_target, testing_multiplier).
#
# art_target is coverage among PLHIV, but stisim fills it only from
# `diagnosed & ~on_art` -- so asking for 0.95 does NOT deliver 0.95, it delivers
# every diagnosed person. 031 obs 4 measured exactly that: on-ART-given-aware
# hit 1.000 in the bound strata. That is why this rung is called LINKAGE and not
# "95% coverage": 0.95 is a way of saying "link everyone you have found".
#
# vls_target must EXCEED the observed 0.959 (f) / 0.967 (m) or _ramp_long raises,
# by design -- it refuses to build an arm that would do nothing. 0.99 is about
# the best any SHIMS3 stratum achieves (50+ reads 0.989 / 0.992).
#
# testing multiplies the FITTED rates (0.5 / 0.6). The ramp plateau is 0.5, so
# TEST_RATE_M * 4 = 2.0 clips to a probability of 1.0 -- every undiagnosed man
# tested each year. C5 is therefore the saturation point, not an arbitrary end.
#
# *** 033 EXTENDS 032'S LADDER. ***
# 032 stopped at an ART target of 0.95 and concluded the cascade "saturates" at
# 0.948 suppression. That conclusion was wrong: decomposing the residual showed
# 78% of it was people DIAGNOSED BUT NOT ON ART, and delivered coverage was
# 0.958 against a target of 0.950 -- the TARGET was binding, not any ceiling.
# Worse, `aware_not_on_art` got WORSE as testing improved (0.021 -> 0.040),
# because once awareness reaches 0.998 a 0.95-of-PLHIV target only demands 95.2%
# of the diagnosed, so linkage falls back.
#
# Rungs C6-C8 raise the ART target itself, then the VLS target, to find where
# the model actually stops. C0-C5 are re-run unchanged so the surface is one
# continuous object rather than two spliced ones.
#
# *** 034: A LADDER ESWATINI COULD ACTUALLY CLIMB. ***
# 033 obs 9 showed its own top rungs demand 0.999 aware / 0.998 linked / 0.999
# suppressed from EVERY age-sex group, against SHIMS3's measured 0.937 / 0.973 /
# 0.962 and a best-ever single stratum of 0.983 / 0.998 / 0.992. Those rungs are
# a bound, not policy. These are built from targets that have been observed.
#
# Every target is applied PER STRATUM. _ramp_long takes a per-stratum max(), so
# a target raises the laggards and leaves the leaders alone -- which is the
# point, because Eswatini's gap is equity, not the aggregate. The national
# figures already clear the second and third 95s; men 25-34 do not.
#
# S1_testing_only exists to isolate 033 obs 3: with the ART target untouched,
# scaling testing cut the undiagnosed by 1.18 pp but raised diagnosed-but-
# untreated by 0.94 pp, leaving suppression essentially flat. Diagnosis converts
# into suppression only if treatment expands to absorb it.
#
# ART target 0.9025 = 0.95 aware x 0.95 linked. The model's ART input is
# coverage among ALL PLHIV, not among the diagnosed, so that product is the
# honest translation of UNAIDS 95-95.
#
# VLS is left at baseline in S1 and S2: Eswatini's measured 0.962 already
# exceeds the third 95, and raising it would credit the scenario with a gain
# that is not on the table.
CASCADE = {
    "S0_status_quo":     (None,   None,  1.0),
    "S1_testing_only":   (None,   None,  3.0),
    "S2_unaids_95":      (0.9025, None,  3.0),
    "S3_99_96_98":       (0.96,   0.98,  4.0),
}
# The theoretical bound (ART 0.999 / VLS 0.999) lived here as S4 and has been
# removed: exp 033 already established it, it is not a scenario anybody could
# pursue, and carrying it invited it being read as one.

# --- the PrEP axis (y) ---------------------------------------------------------
# Nested and DISJOINT. scenarios.agyw() is `female & 15-24` and does NOT exclude
# FSW, so the shipped helpers overlap. Each rung below subtracts what the lower
# rungs already cover, and each Prep instance gets its own name -- two instances
# enrolling the same agent would double-count person-years and corrupt the
# efficiency metric, which was 031's most robust output.


def _fsw(sim):
    return S.fsw(sim)


def _agyw_risk_not_fsw(sim):
    return S.agyw_risk(sim) & ~S.fsw(sim)


def _agyw_rest(sim):
    return S.agyw(sim) & ~S.fsw(sim) & ~S.agyw_risk(sim)


def _women_25_34_rest(sim):
    ppl = sim.people
    return (ppl.female & (ppl.age >= 25) & (ppl.age < 35) & ~S.fsw(sim))


# (eligibility fn, coverage, label) added at each rung
PREP_RUNGS = {
    "P0_none":       [],
    "P1_fsw":        [(_fsw, 0.60, "fsw")],
    "P2_agyw_risk":  [(_fsw, 0.60, "fsw"), (_agyw_risk_not_fsw, 0.30, "agywrisk")],
    "P3_agyw_all":   [(_fsw, 0.60, "fsw"), (_agyw_risk_not_fsw, 0.30, "agywrisk"),
                      (_agyw_rest, 0.30, "agywrest")],
    "P4_women_25_34": [(_fsw, 0.60, "fsw"), (_agyw_risk_not_fsw, 0.30, "agywrisk"),
                       (_agyw_rest, 0.30, "agywrest"),
                       (_women_25_34_rest, 0.30, "w2534")],
}


def make_cell_kwargs(casc, prep):
    art_t, vls_t, test_mult = CASCADE[casc]
    kw = {}

    if art_t is not None or vls_t is not None:
        art_tbl, vls_tbl = S.baseline_tables()
        if art_t is not None:
            kw["art_coverage"] = S.art_coverage_target(
                art_tbl, art_t, SCEN_START, SCEN_REACH, stop=STOP)
        if vls_t is not None:
            kw["art_vls_coverage"] = S.cascade_target(
                vls_tbl, vls_t, SCEN_START, SCEN_REACH, stop=STOP)

    if test_mult != 1.0:
        # test_boost, NOT test_rate_m/f. Multiplying the fitted rates scales
        # testing back to 1990 and changes the epidemic BEFORE the scenario
        # starts -- this experiment's first run did that, and its cascade arms
        # entered 2026 with up to 13% lower incidence than status quo before
        # their scenario had done anything. Those results were discarded.
        kw["test_boost"] = test_mult
        kw["test_boost_start"] = SCEN_START
        kw["test_boost_reach"] = SCEN_REACH

    rungs = PREP_RUNGS[prep]
    if rungs:
        kw["extra_interventions"] = [
            S.lenacapavir(cov, elig, start=SCEN_START, scale_to=SCEN_REACH,
                          name=f"prep_{lab}")
            for elig, cov, lab in rungs]
    return kw


def cell_fingerprint(casc, prep):
    art_t, vls_t, test_mult = CASCADE[casc]
    parts = [casc, prep, repr((art_t, vls_t, test_mult)),
             f"len_eff={S.LEN_EFF}", f"len_dur={S.LEN_DUR_MONTHS}",
             f"scen={SCEN_START}-{SCEN_REACH}", f"stop={STOP}",
             f"base_test={TEST_RATE_M},{TEST_RATE_F}",
             ";".join(f"{lab}:{cov}" for _, cov, lab in PREP_RUNGS[prep])]
    kw = make_cell_kwargs(casc, prep)
    for k in sorted(kw):
        v = kw[k]
        if isinstance(v, pd.DataFrame):
            parts.append(f"{k}:{pd.util.hash_pandas_object(v, index=False).sum()}")
        elif k == "extra_interventions":
            for iv in v:
                parts.append(f"{k}:{iv.name}:{sorted(iv.pars.items(), key=str)}")
        else:
            parts.append(f"{k}:{v!r}")
    return hashlib.sha1("|".join(parts).encode()).hexdigest()[:8]


def run_key(casc, prep, seed):
    return f"{casc}__{prep}__{cell_fingerprint(casc, prep)}__{seed:03d}"


def _one(casc, prep, seed):
    path = SIM_DIR / f"{run_key(casc, prep, seed)}.parquet"
    if path.exists():
        return
    hiv_pars, network_pars = build_pars(BEST)
    t0 = time.perf_counter()
    sim = make_sim(seed=seed, start=START, stop=STOP, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars,
                   analyzers=[PopByAgeSex(), Cascade(), CascadeByAge()],
                   **make_cell_kwargs(casc, prep))
    sim.run()
    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns
            if c == "timevec" or c.startswith(KEEP_PREFIX)]
    out = df[keep].copy()
    if len(keep) < 20:
        raise RuntimeError(
            f"only {len(keep)} columns matched {KEEP_PREFIX} -- the export "
            f"naming changed; grid outputs would be empty.")
    out["cascade"], out["prep"], out["seed"] = casc, prep, seed
    out.to_parquet(path, index=False)
    print(f"  {run_key(casc, prep, seed)}  {time.perf_counter()-t0:.0f}s",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_seeds", type=int, default=N_SEEDS)
    ap.add_argument("--n_workers", type=int, default=110)
    ap.add_argument("--cascades", nargs="*", default=None,
                    help="subset of cascade rungs to run, e.g. S0_status_quo "
                         "S2_unaids_95. Cells already cached are skipped, so "
                         "narrowing a run never discards completed work.")
    ap.add_argument("--smoke", action="store_true",
                    help="corners of the grid, 1 seed -- proves the harness")
    a = ap.parse_args()

    cascades, preps = list(CASCADE), list(PREP_RUNGS)
    if a.cascades:
        unknown = [c for c in a.cascades if c not in CASCADE]
        if unknown:
            raise SystemExit(f"unknown cascade(s) {unknown}; "
                             f"available: {list(CASCADE)}")
        cascades = [c for c in cascades if c in a.cascades]
    seeds = list(range(1, a.n_seeds + 1))
    if a.smoke:
        cascades = ["S0_status_quo", "S3_99_96_98"]
        preps = ["P0_none", "P4_women_25_34"]
        seeds = [1]

    todo = [(c, p, s) for c in cascades for p in preps for s in seeds]
    print(f"{len(todo)} runs: {len(cascades)} cascade x {len(preps)} PrEP x "
          f"{len(seeds)} seeds, {a.n_workers} workers, {START}-{STOP}")
    for c in cascades:
        art_t, vls_t, tm = CASCADE[c]
        print(f"    {c:18} art={art_t} vls={vls_t} test_x{tm}")
    for p in preps:
        print(f"    {p:18} {[lab for _, _, lab in PREP_RUNGS[p]]}")

    t0 = sc.tic()
    if a.n_workers <= 1:
        for t in todo:
            _one(*t)
    else:
        sc.parallelize(_one, iterarg=todo, ncpus=a.n_workers)
    sc.toc(t0, label="all runs")
    print(f"\n{len(list(SIM_DIR.glob('*.parquet')))} parquet in {SIM_DIR}")


if __name__ == "__main__":
    main()
