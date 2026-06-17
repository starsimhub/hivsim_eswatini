"""
Step 1 — Diagnostic: does age_diff_pars drive realized partner-age gaps?

Run 5 seeds × 3 sim configs (age_diff_pars zeros / defaults / doubled),
extract realized partner-age distribution from network edges at sim end,
plot per-config histograms. If they're nearly identical we've confirmed
the rank-based matching does not respect the parameter.
"""

from __future__ import annotations

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
REPO_ROOT = EXP_DIR.parents[1]
OUT_DIR = EXP_DIR / "outputs"
OUT_DIR.mkdir(exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
from run_sims import make_sim  # noqa: E402

import json

N_SEEDS = 5
SIM_STOP = 2016  # Full epidemic + ART scale-up; gives mature network

# Variant label suffix for output files (set via EXP011_LABEL env var by
# orchestration scripts; empty string means use canonical filenames).
LABEL = os.environ.get("EXP011_LABEL", "")

# Three sim configurations — same across (woman age × risk level) cells per config.
# Format: dict of (mu, sd) per (age_group, risk_level).
def make_age_diff_pars(mu: float, sd: float) -> dict:
    return dict(
        teens=[(mu, sd), (mu, sd), (mu, sd)],
        young=[(mu, sd), (mu, sd), (mu, sd)],
        adult=[(mu, sd), (mu, sd), (mu, sd)],
    )

# Optional override: EXP011_CONFIGS_JSON env var lets an orchestration script
# substitute custom CONFIGS (e.g., DHS-derived per-woman-age × risk-level pars
# instead of the uniform mu/sd diagnostic sweep).
_configs_override = os.environ.get("EXP011_CONFIGS_JSON")
if _configs_override:
    CONFIGS = {k: {gk: [tuple(p) for p in pairs] for gk, pairs in v.items()}
               for k, v in json.loads(_configs_override).items()}
else:
    CONFIGS = {
        "A_zeros":    make_age_diff_pars(mu=0,  sd=3),
        "B_defaults": make_age_diff_pars(mu=7,  sd=3),
        "C_doubled":  make_age_diff_pars(mu=14, sd=3),
    }


def realized_pairs_from_sim(sim) -> pd.DataFrame:
    """Pull current network edges with stored age_p1/age_p2 metadata."""
    sn = sim.networks.structuredsexual
    edges = sn.edges
    age_p1 = np.asarray(edges.age_p1)
    age_p2 = np.asarray(edges.age_p2)
    valid = ~(np.isnan(age_p1) | np.isnan(age_p2))
    df = pd.DataFrame({
        "age_man":   age_p1[valid],   # p1 is male in StructuredSexual edges
        "age_woman": age_p2[valid],
    })
    df["gap"] = df["age_man"] - df["age_woman"]
    return df


def run_one(config_name: str, age_diff_pars: dict, seed: int) -> pd.DataFrame:
    sim = make_sim(seed=seed, stop=SIM_STOP, verbose=-1,
                   network_pars={"age_diff_pars": age_diff_pars})
    sim.run()
    df = realized_pairs_from_sim(sim)
    df["config"] = config_name
    df["seed"] = seed
    return df


def main() -> None:
    work = []
    for config_name, age_diff_pars in CONFIGS.items():
        for seed in range(N_SEEDS):
            work.append((config_name, age_diff_pars, seed))

    print(f"Running {len(work)} sims ({len(CONFIGS)} configs × {N_SEEDS} seeds)...")
    t0 = sc.tic()
    dfs = sc.parallelize(run_one, iterarg=work)
    sc.toc(t0, label="diagnostic ensemble")

    edges = pd.concat(dfs, ignore_index=True)
    suffix = f"_{LABEL}" if LABEL else ""
    edges.to_csv(OUT_DIR / f"rank_test_realized_gaps_raw{suffix}.csv", index=False)

    # Aggregate: per-config mean/std of gap, plus by woman age bin
    rows = []
    woman_bins = [(15, 25), (25, 35), (35, 50)]
    for config in CONFIGS:
        sub = edges[edges["config"] == config]
        rows.append({
            "config": config,
            "woman_age_bin": "15-49",
            "n": len(sub),
            "mean_gap": sub["gap"].mean(),
            "std_gap":  sub["gap"].std(),
            "median_gap": sub["gap"].median(),
        })
        for lo, hi in woman_bins:
            ssub = sub[(sub["age_woman"] >= lo) & (sub["age_woman"] < hi)]
            rows.append({
                "config": config,
                "woman_age_bin": f"{lo}-{hi-1}",
                "n": len(ssub),
                "mean_gap": ssub["gap"].mean() if len(ssub) else np.nan,
                "std_gap":  ssub["gap"].std()  if len(ssub) else np.nan,
                "median_gap": ssub["gap"].median() if len(ssub) else np.nan,
            })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / f"rank_test_realized_gaps{suffix}.csv", index=False)
    print("\nSummary (mean realized gap, by config and woman age bin):")
    pivot = summary.pivot(index="woman_age_bin", columns="config", values="mean_gap")
    print(pivot.round(2))

    # Plot: histograms overlay per config (gap distribution)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    bins = np.arange(-15, 31, 1)
    default_colors = {"A_zeros": "C2", "B_defaults": "C0", "C_doubled": "C3"}
    default_labels = {"A_zeros": "A: mu=0",
                      "B_defaults": "B: mu=7 (default)",
                      "C_doubled": "C: mu=14"}
    fallback_palette = [f"C{i}" for i in range(10)]
    for idx, config in enumerate(CONFIGS):
        sub = edges[edges["config"] == config]
        mu_actual = sub["gap"].mean()
        color = default_colors.get(config, fallback_palette[idx % len(fallback_palette)])
        label = default_labels.get(config, config)
        ax.hist(sub["gap"], bins=bins, density=True, alpha=0.45,
                color=color,
                label=f"{label}  →  realized mu={mu_actual:.1f}")
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel("Partner age − woman age (years)")
    ax.set_ylabel("Density")
    ax.set_title(f"Realized partner-age gap by config\n"
                 f"({N_SEEDS} seeds × {len(CONFIGS)} configs at sim end {SIM_STOP})")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Mean gap by woman age bin
    ax = axes[1]
    fine_bins = [(15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 50)]
    for idx, config in enumerate(CONFIGS):
        sub = edges[edges["config"] == config]
        mids, mus = [], []
        for lo, hi in fine_bins:
            s2 = sub[(sub["age_woman"] >= lo) & (sub["age_woman"] < hi)]
            if len(s2):
                mids.append((lo + hi) / 2)
                mus.append(s2["gap"].mean())
        color = default_colors.get(config, fallback_palette[idx % len(fallback_palette)])
        label = default_labels.get(config, config)
        ax.plot(mids, mus, marker="o", color=color, label=label, lw=2)

    # DHS reference
    dhs_summary = pd.read_csv(OUT_DIR / "dhs_partner_age_summary.csv")
    dhs_5yr = dhs_summary[dhs_summary["grouping"] == "5yr"]
    dhs_mids = (dhs_5yr["age_lo"] + dhs_5yr["age_hi"]) / 2
    ax.plot(dhs_mids, dhs_5yr["mu_gap"], marker="s", ls="--",
            color="black", label="DHS Eswatini 2006-07", lw=2)

    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel("Woman's age (years)")
    ax.set_ylabel("Mean realized partner-age gap (years)")
    ax.set_title("Mean realized gap by woman age bin\n(does config change the output?)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT_DIR / f"rank_test_realized_gaps{suffix}.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote {OUT_DIR / f'rank_test_realized_gaps{suffix}.png'}")


if __name__ == "__main__":
    main()
