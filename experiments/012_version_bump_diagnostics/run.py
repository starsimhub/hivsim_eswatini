"""Exp 012 — version-bump diagnostics on the updated stack.

Runs the model unchanged on starsim 3.5.0 / stisim 1.5.8 for N seeds and
regenerates the two standard diagnostics (fit dashboard + network dashboard)
by reusing the plotting functions in the repo-root plot_dashboard.py.

Figures -> figures/ ; collected multi-seed data + version stamp -> outputs/.

Usage (from repo root):
    python experiments/012_version_bump_diagnostics/run.py --n_seeds 10
"""

import argparse
import pathlib
import json
import sys

import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

# Put the repo root on sys.path so the repo-root plot_dashboard / run_sims
# modules import when this script is launched from a subdirectory.
REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from plot_dashboard import (
    run_sims, collect, plot_fit_dashboard, plot_network_dashboard,
)

HERE = pathlib.Path(__file__).parent
FIG = HERE / "figures"
OUT = HERE / "outputs"
DATA_DIR = "data"
CALIB_DIR = "calibration_data"
LABEL = "012_new_stack"


def main(n_seeds):
    FIG.mkdir(exist_ok=True)
    OUT.mkdir(exist_ok=True)

    # Stamp the exact versions this baseline was produced on.
    stamp = {"starsim": ss.__version__, "stisim": sti.__version__, "n_seeds": n_seeds}
    (OUT / "version_stamp.json").write_text(json.dumps(stamp, indent=2))
    print("Version stamp:", stamp)

    print(f"Running {n_seeds} sims on the updated stack...")
    sims = run_sims(n_seeds)
    print("Collecting data...")
    data = collect(sims)

    # Persist collected multi-seed data so figures regenerate without re-running.
    sc.save(str(OUT / "dashboard_data.obj"), data)

    # Targets (same loading as plot_dashboard.__main__).
    targets = pd.read_csv(f"{DATA_DIR}/eswatini_hiv_calib.csv")
    phia_art = pd.read_csv(f"{CALIB_DIR}/art_coverage_by_age_sex.csv")
    inc_calib = pd.read_csv(f"{CALIB_DIR}/incidence_by_sex.csv")
    prev_calib = pd.read_csv(f"{CALIB_DIR}/prevalence_by_age_sex.csv")
    prev_2021 = pd.read_csv(f"{CALIB_DIR}/prevalence_2021_VALIDATION_ONLY.csv")
    prev_calib = pd.concat([prev_calib, prev_2021], ignore_index=True)
    condom_data = pd.read_csv(f"{DATA_DIR}/condom_use.csv")
    vmmc_data = pd.read_csv(f"{DATA_DIR}/vmmc_coverage.csv")

    print("Plotting fit dashboard...")
    plot_fit_dashboard(data, targets, phia_art, inc_calib, prev_calib, vmmc_data,
                       LABEL, n_seeds, outdir=str(FIG))
    print("Plotting network dashboard...")
    plot_network_dashboard(data, condom_data, LABEL, n_seeds, outdir=str(FIG))

    print("Done. Figures ->", FIG)
    for p in sorted(FIG.glob("*.png")):
        print("   ", p.name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_seeds", type=int, default=10)
    args = parser.parse_args()
    main(args.n_seeds)
