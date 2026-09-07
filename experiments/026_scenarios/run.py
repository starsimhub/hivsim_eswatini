"""Exp 026 -- decision scenarios: long-acting PrEP vs the treatment cascade.

Arms x parameter sets x seeds. Per-run parquet, so an interrupted job resumes
free and a reclaimed spot VM costs nothing already done.

Usage (from the repo ROOT -- run_sims reads data/ by relative path)
  python experiments/026_scenarios/run.py --smoke              # 2 arms, 1 seed
  python experiments/026_scenarios/run.py --n_workers 120
  python experiments/026_scenarios/run.py --params nroy        # once 025 lands

Why PARAM_SETS is a list
------------------------
The calibration will not finish before the 1 Oct abstract deadline, so arms run
at 024's design row 868 -- the single draw that fit 48/48 targets. Making the
interface a LIST from the start means swapping in the wave-2 NROY sample is a
`--params nroy` flag rather than a rewrite. That is the only reason this looks
over-engineered for one parameter set.

With one set, differences between arms carry SEED noise only, not parameter
uncertainty. The abstract must say so and must not report a credible interval.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import sciris as sc

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from run_sims import make_sim                     # noqa: E402
from analyzers import PopByAgeSex, Cascade        # noqa: E402
import scenarios as S                             # noqa: E402

OUT = HERE / "outputs"
SIM_DIR = OUT / "sims"
for d in (OUT, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

START, STOP = 1985, 2041      # tables run to 2040; +1 so 2040 is complete
SCEN_START, SCEN_REACH = 2026, 2030
N_SEEDS = 10
KEEP_PREFIX = ("popagesex.", "cascade.", "hiv.")  # verified at runtime below

# 024 design row 868. Kept in natural units; identical to
# experiments/025_hm_wave2/run.py::TRUTH_NATURAL.
BEST = dict(beta_m2f=0.014161835611057863, rel_beta_f2m=0.3097261030208422,
            s_f_young=2.224000314771824, age_gap_shift=0.0699801603950227,
            age_gap_sd_mult=1.5001144712607821, prop_f0=0.608751696061895,
            prop_m0=0.4742150061435999)
PARAM_SETS = {"best": [BEST]}

AGE_DIFF_BASE = {"teens": [(7, 3), (6, 3), (5, 1)],
                 "young": [(8, 3), (7, 3), (5, 2)],
                 "adult": [(8, 3), (7, 3), (5, 2)]}
CONC_BASE = {"f1_conc": 0.15, "f2_conc": 0.25, "m1_conc": 0.15, "m2_conc": 0.5}
REL_INIT_PREV, CONC_MULT = 0.2, 1.0


def build_pars(p):
    """Identical to 025's _one_point. Copied, not re-derived: if the two
    diverge, scenario arms stop being comparable to the calibration that
    justified their parameters."""
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


# --- Arms ----------------------------------------------------------------------
# Each entry is (prep_spec, cascade_spec). prep_spec is None or
# (coverage, eligibility); cascade_spec is None or (art_target, vls_target).
#
# The cascade axis is ART COVERAGE, not suppression: Eswatini's suppression
# among the treated is already 0.959/0.967, so a "third 95" arm is a no-op --
# scenarios._ramp_long raises on it. See README.
ARMS = {
    "baseline":          (None,                    None),
    "prep_agyw":         ((0.30, "agyw"),          None),
    "prep_agyw_risk":    ((0.30, "agyw_risk"),     None),
    "prep_fsw":          ((0.60, "fsw"),           None),
    "art_95":            (None,                    (0.95, None)),
    "both":              ((0.30, "agyw"),          (0.95, None)),
    "prep_at_high_art":  ((0.30, "agyw"),          (0.95, None)),
}
# prep_at_high_art is `both` with the cascade already achieved when PrEP starts,
# rather than ramping alongside it. The headline quantity is
# prep_at_high_art - art_95: what prevention buys AFTER the coverage gap closes.
ARM_OVERRIDES = {"prep_at_high_art": dict(cascade_start=2020, cascade_reach=2025)}


def make_arm_kwargs(arm):
    prep_spec, casc = ARMS[arm]
    ov = ARM_OVERRIDES.get(arm, {})
    kw = {}
    if prep_spec is not None:
        cov, elig = prep_spec
        kw["extra_interventions"] = [
            S.lenacapavir(cov, elig, start=SCEN_START, scale_to=SCEN_REACH)]
    if casc is not None:
        art_t, vls_t = casc
        art_tbl, vls_tbl = S.baseline_tables()
        c_start = ov.get("cascade_start", SCEN_START)
        c_reach = ov.get("cascade_reach", SCEN_REACH)
        if art_t is not None:
            kw["art_coverage"] = S.art_coverage_target(
                art_tbl, art_t, c_start, c_reach, stop=STOP)
        if vls_t is not None:
            kw["art_vls_coverage"] = S.cascade_target(
                vls_tbl, vls_t, c_start, c_reach, stop=STOP)
    return kw


def run_key(arm, pset, pidx, seed):
    return f"{arm}__{pset}{pidx:03d}__{seed:03d}"


def _one(arm, pset, pidx, seed):
    path = SIM_DIR / f"{run_key(arm, pset, pidx, seed)}.parquet"
    if path.exists():
        return
    p = PARAM_SETS[pset][pidx]
    hiv_pars, network_pars = build_pars(p)
    t0 = time.perf_counter()
    sim = make_sim(seed=seed, start=START, stop=STOP, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars,
                   analyzers=[PopByAgeSex(), Cascade()],
                   **make_arm_kwargs(arm))
    sim.run()
    # resample/use_years/sep are load-bearing, not cosmetic: without sep="."
    # the column names do not match KEEP_PREFIX and the filter silently keeps
    # nothing, and without use_years=True timevec stays a datetime. Identical to
    # 024's export so the two are directly comparable.
    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns
            if c == "timevec" or c.startswith(KEEP_PREFIX)]
    out = df[keep].copy()
    if len(keep) < 20:
        raise RuntimeError(
            f"only {len(keep)} columns matched {KEEP_PREFIX} -- the export "
            f"naming changed; scenario outputs would be empty. Got: "
            f"{sorted(df.columns)[:12]}")
    out["arm"], out["pset"], out["pidx"], out["seed"] = arm, pset, pidx, seed
    out.to_parquet(path, index=False)
    print(f"  {run_key(arm, pset, pidx, seed)}  {time.perf_counter()-t0:.0f}s",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default="best", choices=sorted(PARAM_SETS))
    ap.add_argument("--n_seeds", type=int, default=N_SEEDS)
    ap.add_argument("--n_workers", type=int, default=8)
    ap.add_argument("--arms", nargs="*", default=None)
    ap.add_argument("--smoke", action="store_true",
                    help="2 arms x 1 param set x 1 seed, to prove the harness")
    a = ap.parse_args()

    arms = a.arms or list(ARMS)
    seeds = list(range(1, a.n_seeds + 1))
    if a.smoke:
        arms, seeds = ["baseline", "art_95"], [1]

    todo = [(arm, a.params, i, s)
            for arm in arms
            for i in range(len(PARAM_SETS[a.params]))
            for s in seeds]
    print(f"{len(todo)} runs: {len(arms)} arms x "
          f"{len(PARAM_SETS[a.params])} param set(s) x {len(seeds)} seed(s), "
          f"{a.n_workers} workers, {START}-{STOP}")
    for arm in arms:
        print(f"    {arm:20} {make_arm_kwargs(arm).keys()}")

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
