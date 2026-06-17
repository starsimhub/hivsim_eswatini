"""
Build a consolidated comparison plot from the variant CSVs produced by
run_pr477_comparison.py. Reads every available variant; gracefully skips
any whose CSV is missing (so partial runs still produce useful output).

Two figures:
    figures/algorithm_comparison_summary.png — overlay of mean realized gap
        by woman age bin across all variants, plus DHS Eswatini overlay.
    figures/algorithm_comparison_grid.png — small-multiples histogram grid
        (one panel per variant) of the realized partner-age gap distribution.
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EXP_DIR = Path(__file__).resolve().parent
OUT_DIR = EXP_DIR / "outputs"
FIG_DIR = EXP_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)


# (label, human_label, color)
VARIANTS = [
    ("baseline_no_fix",  "Baseline (no fix)",          "#888888"),
    ("E_gaussian_patch", "Gaussian patch",             "#0066cc"),
    ("A_pr477_md1",      "PR #477 (md=1)",             "#cc0033"),
    ("B_pr477_md3",      "PR #477 (md=3)",             "#dd7733"),
    ("C_pr477_md5",      "PR #477 (md=5)",             "#cc9900"),
    ("D_pr477_dhs",      "PR #477 + DHS pars (md=1)",  "#7733aa"),
]


def load_raw(label: str) -> pd.DataFrame | None:
    """Variant raw CSV has one row per realized edge."""
    if label == "baseline_no_fix":
        # Special-case: baseline's raw filename was the original "rank_test_realized_gaps_raw.csv"
        candidates = [
            OUT_DIR / "rank_test_realized_gaps_raw_baseline_no_fix.csv",
            OUT_DIR / "rank_test_realized_gaps_raw.csv",
        ]
        for p in candidates:
            if p.exists():
                return pd.read_csv(p)
        return None
    p = OUT_DIR / f"rank_test_realized_gaps_raw_{label}.csv"
    if p.exists():
        return pd.read_csv(p)
    return None


def load_summary(label: str) -> pd.DataFrame | None:
    """Variant summary CSV has aggregated mean/std/median by woman age bin."""
    if label == "baseline_no_fix":
        candidates = [
            OUT_DIR / "rank_test_realized_gaps_baseline_no_fix.csv",
            OUT_DIR / "rank_test_realized_gaps.csv",
        ]
        for p in candidates:
            if p.exists():
                return pd.read_csv(p)
        return None
    p = OUT_DIR / f"rank_test_realized_gaps_{label}.csv"
    if p.exists():
        return pd.read_csv(p)
    return None


def main() -> int:
    # Load DHS reference
    dhs_path = OUT_DIR / "dhs_partner_age_summary.csv"
    dhs_5yr = None
    if dhs_path.exists():
        dhs = pd.read_csv(dhs_path)
        dhs_5yr = dhs[dhs["grouping"] == "5yr"].copy()
        dhs_5yr["mid"] = (dhs_5yr["age_lo"] + dhs_5yr["age_hi"]) / 2

    available = []
    for label, human, color in VARIANTS:
        raw = load_raw(label)
        summary = load_summary(label)
        if raw is None and summary is None:
            print(f"skip {label}: no data")
            continue
        available.append((label, human, color, raw, summary))

    if not available:
        print("No variant data found in outputs/.")
        return 1

    # ---- Figure 1: mean realized gap by woman age, per variant ----
    fig, ax = plt.subplots(figsize=(11, 6.5))
    fine_bins = [(15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 50)]

    # For each variant, focus on the "B_defaults" / "default" config if present,
    # else aggregate across all configs in the variant.
    for label, human, color, raw, summary in available:
        if raw is None:
            continue
        # Prefer B_defaults config (mu=7) where present
        if "config" in raw.columns and "B_defaults" in set(raw["config"]):
            sub_full = raw[raw["config"] == "B_defaults"]
        else:
            sub_full = raw  # DHS variant has its own config name
        mids, mus = [], []
        for lo, hi in fine_bins:
            s = sub_full[(sub_full["age_woman"] >= lo) & (sub_full["age_woman"] < hi)]
            if len(s) > 10:  # avoid noisy points
                mids.append((lo + hi) / 2)
                mus.append(s["gap"].mean())
        ax.plot(mids, mus, marker="o", color=color, label=human, lw=2.0, alpha=0.9)

    if dhs_5yr is not None and len(dhs_5yr):
        ax.plot(dhs_5yr["mid"], dhs_5yr["mu_gap"], marker="s", ls="--",
                color="black", lw=2.5, label="DHS Eswatini 2006-07")

    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel("Woman's age (years)")
    ax.set_ylabel("Mean realized partner-age gap (years)")
    ax.set_title("Algorithm comparison: mean realized partner-age gap by woman age\n"
                 "(B_defaults config where present; DHS variant uses its own pars)")
    ax.legend(fontsize=9, loc="best")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out1 = FIG_DIR / "algorithm_comparison_summary.png"
    fig.savefig(out1, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out1}")

    # ---- Figure 2: small-multiples histogram grid ----
    n = len(available)
    ncols = 2
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 3.4 * nrows), sharex=True, sharey=True)
    axes = np.atleast_2d(axes).reshape(-1)
    bins = np.arange(-15, 31, 1)
    for i, (label, human, color, raw, summary) in enumerate(available):
        ax = axes[i]
        if raw is None:
            ax.set_visible(False)
            continue
        for cfg, cfg_color, alpha in [
            ("A_zeros", "#2ca02c", 0.45),
            ("B_defaults", "#1f77b4", 0.55),
            ("C_doubled", "#d62728", 0.45),
            ("DHS_eswatini", color, 0.7),
        ]:
            if "config" not in raw.columns:
                continue
            sub = raw[raw["config"] == cfg]
            if not len(sub):
                continue
            mu = sub["gap"].mean()
            ax.hist(sub["gap"], bins=bins, density=True, alpha=alpha,
                    color=cfg_color, label=f"{cfg}: μ={mu:.1f}")
        ax.axvline(0, color="grey", lw=0.5, ls="--")
        ax.set_title(human, fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    for ax in axes[-ncols:]:
        ax.set_xlabel("Partner age − woman age (years)")
    for r in range(nrows):
        axes[r * ncols].set_ylabel("Density")
    fig.suptitle("Realized partner-age gap distribution — by algorithm variant", y=1.00)
    fig.tight_layout()
    out2 = FIG_DIR / "algorithm_comparison_grid.png"
    fig.savefig(out2, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out2}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
