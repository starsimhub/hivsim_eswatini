"""Exp 029 -- the care cascade decomposed into its three conditional steps.

Runs the calibration baseline at one parameter point and records, per 5-year
age band and sex, the four counts the 95-95-95 cascade is built from:
infected, diagnosed, on ART, virally suppressed.

Why one point and not the NROY: exp 028 showed the NROY's marginals are
essentially the prior, so sampling it would widen every interval here without
making any of them more informative about the cascade. The cascade is driven by
the testing ramps and the ART coverage table -- neither of which is a calibrated
parameter -- so a point run answers the question and an ensemble would not.

Seeds, not draws, are therefore the only uncertainty shown.

Usage
  python experiments/029_cascade_audit/run.py
  python experiments/029_cascade_audit/run.py --n_seeds 4 --n_workers 4   # smoke
Run from the repo ROOT -- run_sims reads data/ by relative path.
"""

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import sciris as sc

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from run_sims import make_sim                                   # noqa: E402
from analyzers import PopByAgeSex, Cascade, CascadeByAge        # noqa: E402

OUT = HERE / "outputs"
SIM_DIR = OUT / "sims"
for d in (OUT, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

# 20 000 agents is exp 020's floor and exp 028's setting. Held here because the
# small strata this experiment reports on -- men 15-24 living with HIV are ~50
# unweighted in SHIMS3 -- are exactly the ones that go noisy below it.
N_AGENTS = 20_000
START, STOP = 1985, 2031      # 2030 is 026's scenario horizon; +1 so it completes
N_SEEDS = 20                  # up from 026's 10: age x sex x cascade-step is a
                              # finer stratification than 026 ever reported on

# 024 design row 868 -- the draw that fit 48/48 targets, and the point exp 026
# built its scenarios on. Identical to 026/capacity_check.py::BEST so the two
# experiments describe the same model.
BEST = dict(beta_m2f=0.014161835611057863, rel_beta_f2m=0.3097261030208422,
            s_f_young=2.224000314771824, age_gap_shift=0.0699801603950227,
            age_gap_sd_mult=1.5001144712607821, prop_f0=0.608751696061895,
            prop_m0=0.4742150061435999)

AGE_DIFF_BASE = {"teens": [(7, 3), (6, 3), (5, 1)],
                 "young": [(8, 3), (7, 3), (5, 2)],
                 "adult": [(8, 3), (7, 3), (5, 2)]}
REL_INIT_PREV = 0.2
CONC_BASE = {"f1_conc": 0.15, "f2_conc": 0.25, "m1_conc": 0.15, "m2_conc": 0.5}

KEEP_PREFIX = ("popagesex.", "cascade.", "cascadeage.", "hiv.")


def build_pars(p):
    """Translate a parameter set into hiv_pars / network_pars.

    Copied from 026/capacity_check.py rather than imported: that file is an
    experiment artefact and 026 is closed. If these diverge, 029's cascade stops
    describing the model 026 costed its scenarios on.
    """
    hiv_pars = dict(
        beta_m2f=p["beta_m2f"], rel_beta_f2m=p["rel_beta_f2m"],
        rel_init_prev=REL_INIT_PREV,
        rel_sus_age=[(15, 25, 'f', p["s_f_young"]),
                     (25, 50, 'f', 1.0), (15, 50, 'm', 1.0)])
    shift, sd_mult = p["age_gap_shift"], p["age_gap_sd_mult"]
    network_pars = dict(
        prop_f0=p["prop_f0"], prop_m0=p["prop_m0"],
        age_diff_pars={g: [(max(m + shift, 1.0), max(s * sd_mult, 0.2))
                           for m, s in v] for g, v in AGE_DIFF_BASE.items()},
        **CONC_BASE)
    return hiv_pars, network_pars


def _one(seed):
    """One seed. Writes its own parquet so a killed run loses only the tail."""
    path = SIM_DIR / f"seed{seed:03d}.parquet"
    if path.exists():
        return str(path)
    t0 = time.perf_counter()
    hiv_pars, network_pars = build_pars(BEST)
    sim = make_sim(seed=seed, start=START, stop=STOP, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars,
                   analyzers=[PopByAgeSex(), Cascade(), CascadeByAge()])
    sim.pars.n_agents = N_AGENTS
    sim.run()
    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns
            if c == "timevec" or c.startswith(KEEP_PREFIX)]
    # Same guard as 026: a rename upstream would silently produce empty output
    # rather than an error, and the figure would be built from nothing.
    n_casc = sum(c.startswith("cascadeage.") for c in keep)
    if n_casc < 100:
        raise RuntimeError(
            f"only {n_casc} cascadeage.* columns survived the filter -- the "
            f"analyzer did not register or export naming changed. "
            f"Sample: {sorted(df.columns)[:12]}")
    out = df[keep].copy()
    out["seed"] = seed
    out.to_parquet(path, index=False)
    print(f"  seed {seed:3d}  {time.perf_counter() - t0:5.0f}s  "
          f"({n_casc} cascade cols)", flush=True)
    return str(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_seeds", type=int, default=N_SEEDS)
    ap.add_argument("--n_workers", type=int, default=10)
    args = ap.parse_args()

    seeds = list(range(args.n_seeds))
    print(f"exp 029: {len(seeds)} seeds, N={N_AGENTS:,}, {START}-{STOP}, "
          f"{args.n_workers} workers", flush=True)
    t0 = time.perf_counter()
    sc.parallelize(_one, seeds, ncpus=args.n_workers)

    # Consolidate before anything downstream touches the results -- per-seed
    # files are an intermediate, not a deliverable.
    frames = [pd.read_parquet(p) for p in sorted(SIM_DIR.glob("seed*.parquet"))]
    allsims = pd.concat(frames, ignore_index=True)
    allsims.to_parquet(OUT / "results.parquet", index=False)
    print(f"done in {time.perf_counter() - t0:.0f}s -> "
          f"outputs/results.parquet  {allsims.shape}", flush=True)


if __name__ == "__main__":
    main()
