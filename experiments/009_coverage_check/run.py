"""
Experiment 009 — Prior predictive coverage check.

Draw 50 parameter sets from the implicit Optuna prior, simulate one
trajectory per draw, and check whether the calibration targets (frozen in
experiment 008) fall inside the ensemble.

See README.md for the full plan and success criteria.
"""

from __future__ import annotations

# Pin BLAS threading for parallel safety, mirroring run_calibrations.py
import os
os.environ.update(
    OMP_NUM_THREADS="1",
    OPENBLAS_NUM_THREADS="1",
    NUMEXPR_NUM_THREADS="1",
    MKL_NUM_THREADS="1",
)

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sciris as sc

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]  # hivsim_eswatini/
OUT_DIR = EXP_DIR / "outputs"
OUT_DIR.mkdir(exist_ok=True)

# Imports from the model repo (which lives at REPO_ROOT)
sys.path.insert(0, str(REPO_ROOT))
from run_sims import make_sim  # noqa: E402

# --- Configuration ---------------------------------------------------------

N_DRAWS = 50
SIM_STOP = 2026  # No need to project past calibration target horizon
DRAW_SEED = 20260507  # Reproducible prior draw

# Implicit uniform prior (from run_calibrations.py)
PRIOR = {
    "hiv.beta_m2f": (0.002, 0.014),
    "hiv.eff_condom": (0.5, 0.9),
    "hiv.rel_dur_on_art": (1.0, 20.0),
    "structuredsexual.prop_f0": (0.55, 0.9),
    "structuredsexual.prop_m0": (0.55, 0.8),
    "structuredsexual.m1_conc": (0.05, 0.3),
}

# Targets file from experiment 008
TARGETS_PATH = REPO_ROOT / "experiments" / "008_calibration_targets" / "outputs" / "calibration_targets.csv"

# PHIA-year age bin → model output column name pattern.
# The model exposes 5-year bins from 15:20 through 30:35 on `hiv.prevalence_*_*_*`,
# and 35:40 through 60:65 on the `hiv_epi.prevalence_*_*_*` analyzer columns.
PREV_COL_TEMPLATE = {
    "hiv": "hiv.prevalence_{sex}_{lo}_{hi}",
    "hiv_epi": "hiv_epi.prevalence_{sex}_{lo}_{hi}",
}


def model_prev_col(sex: str, age_low: int, age_high: int) -> str:
    """Pick the right prevalence column for a (sex, age_low, age_high) triple."""
    s = sex.lower()
    if age_low < 35:
        return PREV_COL_TEMPLATE["hiv"].format(sex=s, lo=age_low, hi=age_high)
    return PREV_COL_TEMPLATE["hiv_epi"].format(sex=s, lo=age_low, hi=age_high)


# --- Prior drawing ---------------------------------------------------------

def draw_prior(n_draws: int, seed: int) -> pd.DataFrame:
    """Uniform LHS-style draws from the prior, returned as a tidy frame."""
    rng = np.random.default_rng(seed)
    cols = {"par_idx": np.arange(n_draws)}
    for name, (lo, hi) in PRIOR.items():
        cols[name] = rng.uniform(lo, hi, size=n_draws)
    return pd.DataFrame(cols)


def split_pars(row: pd.Series) -> tuple[dict, dict]:
    """Split a parameter row into (hiv_pars, network_pars) for make_sim."""
    hiv, net = {}, {}
    for col, val in row.items():
        if col == "par_idx":
            continue
        module, key = col.split(".", 1)
        target = hiv if module == "hiv" else net
        target[key] = float(val)
    return hiv, net


# --- Single sim worker -----------------------------------------------------

def run_one(par_idx: int, hiv_pars: dict, network_pars: dict) -> pd.DataFrame:
    """Run one sim and return a slim per-year dataframe."""
    sim = make_sim(seed=par_idx, stop=SIM_STOP, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars)
    sim.run()
    df = sim.to_df(resample="year", use_years=True, sep=".")
    # Trim to target-relevant columns to keep memory + I/O light
    keep_prefixes = ("timevec", "hiv.new_deaths", "hiv.prevalence_", "hiv_epi.prevalence_")
    keep = [c for c in df.columns if any(c.startswith(p) or c == p for p in keep_prefixes)]
    out = df[keep].copy()
    out["par_idx"] = par_idx
    return out


# --- Coverage assessment ---------------------------------------------------

def assess_coverage(targets: pd.DataFrame, ensemble: pd.DataFrame) -> pd.DataFrame:
    """For each target row, compute simulated 5–95% envelope and in/out flag."""
    rows = []

    # Index ensemble by year for fast lookup
    ens_by_year = {y: g for y, g in ensemble.groupby("timevec")}

    for _, t in targets.iterrows():
        year = int(t["year"])
        if year not in ens_by_year:
            continue
        sub = ens_by_year[year]

        if t["quantity"] == "prevalence":
            col = model_prev_col(t["sex"], int(t["age_low"]), int(t["age_high"]))
            if col not in sub.columns:
                continue
            sim_vals = sub[col].dropna().values
        elif t["quantity"] == "aids_deaths":
            sim_vals = sub["hiv.new_deaths"].dropna().values
        else:
            continue

        if len(sim_vals) == 0:
            continue

        sim_lo = np.percentile(sim_vals, 5)
        sim_hi = np.percentile(sim_vals, 95)
        sim_med = np.percentile(sim_vals, 50)
        in_envelope = sim_lo <= t["value"] <= sim_hi

        rows.append({
            "quantity": t["quantity"],
            "year": year,
            "sex": t["sex"],
            "age_low": t["age_low"],
            "age_high": t["age_high"],
            "obs_value": t["value"],
            "obs_lower": t["lower"],
            "obs_upper": t["upper"],
            "sim_p05": sim_lo,
            "sim_p50": sim_med,
            "sim_p95": sim_hi,
            "in_envelope": in_envelope,
        })

    return pd.DataFrame(rows)


# --- Plots -----------------------------------------------------------------

def plot_prevalence_coverage(targets: pd.DataFrame, ensemble: pd.DataFrame,
                             path: Path) -> None:
    """3 PHIA years × 2 sexes panel grid; ensemble band + observed CIs."""
    prev = targets[targets["quantity"] == "prevalence"].copy()
    phia_years = sorted(prev["year"].unique())
    sexes = ["M", "F"]

    fig, axes = plt.subplots(2, len(phia_years),
                             figsize=(4.5 * len(phia_years), 8),
                             sharey=True)
    if len(phia_years) == 1:
        axes = axes.reshape(2, 1)

    for j, year in enumerate(phia_years):
        for i, sex in enumerate(sexes):
            ax = axes[i, j]
            obs = prev[(prev["year"] == year) & (prev["sex"] == sex)].copy()
            obs["age_mid"] = (obs["age_low"] + obs["age_high"]) / 2
            obs = obs.sort_values("age_mid")

            # Compute ensemble band per age bin
            year_ens = ensemble[ensemble["timevec"] == year]
            ages, p05, p50, p95 = [], [], [], []
            for _, r in obs.iterrows():
                col = model_prev_col(sex, int(r["age_low"]), int(r["age_high"]))
                if col in year_ens.columns:
                    vals = year_ens[col].dropna().values
                    if len(vals):
                        ages.append((r["age_low"] + r["age_high"]) / 2)
                        p05.append(np.percentile(vals, 5))
                        p50.append(np.percentile(vals, 50))
                        p95.append(np.percentile(vals, 95))

            if ages:
                ax.fill_between(ages, p05, p95, color="C0", alpha=0.25,
                                label="Sim 5–95 %")
                ax.plot(ages, p50, color="C0", lw=1.5, label="Sim median")

            err = [obs["value"] - obs["lower"], obs["upper"] - obs["value"]]
            ax.errorbar(obs["age_mid"], obs["value"], yerr=err, fmt="o",
                        color="C3", capsize=3, label="PHIA")

            ax.set_title(f"PHIA {year} — {'Male' if sex == 'M' else 'Female'}")
            ax.set_xlabel("Age (mid-bin)")
            if j == 0:
                ax.set_ylabel("HIV prevalence")
            ax.set_ylim(0, max(0.7, prev["upper"].max() * 1.05))
            ax.grid(alpha=0.3)
            if i == 0 and j == len(phia_years) - 1:
                ax.legend(fontsize=8, loc="upper right")

    fig.suptitle("Coverage check — PHIA prevalence vs. prior ensemble",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_deaths_coverage(targets: pd.DataFrame, ensemble: pd.DataFrame,
                         path: Path) -> None:
    deaths = targets[targets["quantity"] == "aids_deaths"].sort_values("year")

    by_year = ensemble.groupby("timevec")["hiv.new_deaths"]
    yrs = sorted(by_year.groups.keys())
    p05 = [np.percentile(by_year.get_group(y).dropna(), 5) for y in yrs]
    p50 = [np.percentile(by_year.get_group(y).dropna(), 50) for y in yrs]
    p95 = [np.percentile(by_year.get_group(y).dropna(), 95) for y in yrs]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.fill_between(yrs, p05, p95, color="C0", alpha=0.25, label="Sim 5–95 %")
    ax.plot(yrs, p50, color="C0", lw=1.5, label="Sim median")
    ax.errorbar(deaths["year"], deaths["value"],
                yerr=[deaths["value"] - deaths["lower"],
                      deaths["upper"] - deaths["value"]],
                fmt="o", color="C3", capsize=3, ms=4, label="UNAIDS (±15% placeholder)")

    ax.set_title("Coverage check — UNAIDS HIV deaths vs. prior ensemble")
    ax.set_xlabel("Year")
    ax.set_ylabel("New AIDS deaths (annual)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


# --- Main ------------------------------------------------------------------

def main() -> None:
    print(f"Drawing {N_DRAWS} parameter sets from prior...")
    draws = draw_prior(N_DRAWS, DRAW_SEED)
    draws.to_csv(OUT_DIR / "draws.csv", index=False)

    work = []
    for _, row in draws.iterrows():
        hiv_pars, network_pars = split_pars(row)
        work.append((int(row["par_idx"]), hiv_pars, network_pars))

    print(f"Running {len(work)} sims (parallel)...")
    t0 = sc.tic()
    dfs = sc.parallelize(run_one, iterarg=work)
    sc.toc(t0, label="ensemble run")

    ensemble = pd.concat(dfs, ignore_index=True)
    ensemble.to_parquet(OUT_DIR / "ensemble.parquet", index=False)
    print(f"Ensemble: {len(ensemble)} rows × {len(ensemble.columns)} columns "
          f"from {ensemble['par_idx'].nunique()} sims")

    print("Loading targets and assessing coverage...")
    targets = pd.read_csv(TARGETS_PATH)
    coverage = assess_coverage(targets, ensemble)
    coverage.to_csv(OUT_DIR / "coverage_summary.csv", index=False)

    n_total = len(coverage)
    n_in = int(coverage["in_envelope"].sum())
    print(f"\nCoverage: {n_in}/{n_total} target rows inside 5–95% envelope")
    by_q = coverage.groupby("quantity")["in_envelope"].agg(["sum", "count"])
    print(by_q)

    plot_prevalence_coverage(targets, ensemble, OUT_DIR / "coverage_prevalence.png")
    plot_deaths_coverage(targets, ensemble, OUT_DIR / "coverage_deaths.png")

    print(f"\nWrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
