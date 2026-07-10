"""Exp 015 — fit + network dashboards for the fixed (prevalence-target) VMMC arm.

Reuses plot_dashboard.py's run/collect/plot pipeline, but swaps in the in-repo
VMMCPrevalenceTarget class so panel F (VMMC coverage by age) shows the corrected,
age-differentiated coverage rather than the ~100% overshoot of upstream VMMC.

Usage (from repo root):
    python experiments/015_vmmc_prevalence_target/run_dashboard.py --n_seeds 10
"""

import argparse
import pathlib
import sys

import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
import plot_dashboard as pd_dash
from plot_dashboard import DATA_DIR, CALIB_DIR

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))
from vmmc_prevalence_target import VMMCPrevalenceTarget

FIG = str(HERE / "figures")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_seeds", type=int, default=10)
    parser.add_argument("--label", type=str, default="015_fixed_vmmc")
    args = parser.parse_args()

    print(f"Running {args.n_seeds} sims with prevalence-target VMMC...")
    sims = pd_dash.run_sims(args.n_seeds, vmmc_class=VMMCPrevalenceTarget)
    print("Collecting data...")
    data = pd_dash.collect(sims)

    targets = pd.read_csv(f"{DATA_DIR}/eswatini_hiv_calib.csv")
    phia_art = pd.read_csv(f"{CALIB_DIR}/art_coverage_by_age_sex.csv")
    inc_calib = pd.read_csv(f"{CALIB_DIR}/incidence_by_sex.csv")
    prev_calib = pd.read_csv(f"{CALIB_DIR}/prevalence_by_age_sex.csv")
    prev_2021 = pd.read_csv(f"{CALIB_DIR}/prevalence_2021_VALIDATION_ONLY.csv")
    prev_calib = pd.concat([prev_calib, prev_2021], ignore_index=True)
    condom_data = pd.read_csv(f"{DATA_DIR}/condom_use.csv")
    vmmc_data = pd.read_csv(f"{DATA_DIR}/vmmc_coverage.csv")

    pd_dash.plot_fit_dashboard(data, targets, phia_art, inc_calib, prev_calib,
                               vmmc_data, args.label, args.n_seeds, outdir=FIG)
    pd_dash.plot_network_dashboard(data, condom_data, args.label, args.n_seeds, outdir=FIG)
