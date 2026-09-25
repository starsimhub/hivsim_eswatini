"""Exp 030 -- scan the sex-specific routine testing rate against SHIMS3 awareness.

Phase 1 (`--phase scan`): two 1-D scans, k_m with k_f=1.0 and k_f with k_m=1.0,
3 seeds each, to 2022. Near-separable because each sex's testing rate reaches
the other's awareness only through transmission, which is second-order for a
2021 stock. k=(1.0, 1.0) appears in both and is the regression check against 029.

Phase 2 (`--phase confirm --km K --kf K`): the selected pair, 10 seeds, to 2031,
so 029's ceiling and coverage figures can be regenerated on the new model.

Usage
  python experiments/030_testing_rate/run.py --phase scan
  python experiments/030_testing_rate/run.py --phase confirm --km 0.4 --kf 0.7
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

N_AGENTS = 20_000
START = 1985
STOP_FIT, STOP_CONFIRM = 2022, 2031
SEEDS_SCAN, SEEDS_CONFIRM = 3, 10

# Identical to 029/026 -- the draw that fit 48/48 targets. The epidemiological
# parameters are NOT re-fitted here; only the two testing multipliers move.
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

K_M_GRID = [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
K_F_GRID = [0.3, 0.4, 0.5, 0.6, 0.8, 1.0]


def build_pars(p):
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


def _one(km, kf, seed, stop):
    """One point. Writes its own parquet; resumable.

    Takes four positional args rather than a tuple: sc.parallelize splats each
    entry of the job list into the call, so a tuple arrives already unpacked.
    """
    tag = f"km{km:.2f}_kf{kf:.2f}_s{seed:03d}_to{stop}"
    path = SIM_DIR / f"{tag}.parquet"
    if path.exists():
        return
    t0 = time.perf_counter()
    hiv_pars, network_pars = build_pars(BEST)
    sim = make_sim(seed=seed, start=START, stop=stop, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars,
                   test_rate_m=km, test_rate_f=kf,
                   analyzers=[PopByAgeSex(), Cascade(), CascadeByAge()])
    sim.pars.n_agents = N_AGENTS
    sim.run()
    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns
            if c == "timevec" or c.startswith(KEEP_PREFIX)]
    n_casc = sum(c.startswith("cascadeage.") for c in keep)
    if n_casc < 100:
        raise RuntimeError(f"only {n_casc} cascadeage.* columns -- analyzer "
                           f"missing or export renamed")
    out = df[keep].copy()
    out["k_m"], out["k_f"], out["seed"] = km, kf, seed
    out.to_parquet(path, index=False)
    print(f"  {tag}  {time.perf_counter() - t0:5.0f}s", flush=True)


def consolidate(pattern, outfile):
    """Per-seed files are an intermediate; consolidate before anything reads them."""
    frames = [pd.read_parquet(p) for p in sorted(SIM_DIR.glob(pattern))]
    if not frames:
        raise RuntimeError(f"no sims matched {pattern}")
    allsims = pd.concat(frames, ignore_index=True)
    allsims.to_parquet(OUT / outfile, index=False)
    print(f"-> outputs/{outfile}  {allsims.shape}", flush=True)
    return allsims


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["scan", "confirm"], default="scan")
    ap.add_argument("--km", type=float, default=1.0)
    ap.add_argument("--kf", type=float, default=1.0)
    ap.add_argument("--n_workers", type=int, default=10)
    args = ap.parse_args()

    if args.phase == "scan":
        # Two 1-D scans sharing the (1.0, 1.0) corner, which is also the
        # regression check against 029.
        jobs = {(km, 1.0, s, STOP_FIT) for km in K_M_GRID
                for s in range(SEEDS_SCAN)}
        jobs |= {(1.0, kf, s, STOP_FIT) for kf in K_F_GRID
                 for s in range(SEEDS_SCAN)}
        jobs = sorted(jobs)
        print(f"exp 030 scan: {len(jobs)} sims, N={N_AGENTS:,}, "
              f"{START}-{STOP_FIT}, {args.n_workers} workers", flush=True)
        t0 = time.perf_counter()
        sc.parallelize(_one, jobs, ncpus=args.n_workers)
        consolidate(f"*_to{STOP_FIT}.parquet", "scan.parquet")
        print(f"scan done in {time.perf_counter() - t0:.0f}s", flush=True)
    else:
        jobs = [(args.km, args.kf, s, STOP_CONFIRM)
                for s in range(SEEDS_CONFIRM)]
        print(f"exp 030 confirm: k_m={args.km}, k_f={args.kf}, "
              f"{len(jobs)} seeds to {STOP_CONFIRM}", flush=True)
        t0 = time.perf_counter()
        sc.parallelize(_one, jobs, ncpus=args.n_workers)
        consolidate(f"km{args.km:.2f}_kf{args.kf:.2f}_*_to{STOP_CONFIRM}.parquet",
                    "confirm.parquet")
        print(f"confirm done in {time.perf_counter() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
