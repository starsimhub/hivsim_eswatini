"""
Extract empirical partner-age distribution from DHS Eswatini 2006-07.

Source: references/SZIR51FL/szir51fl.dat (DHS V Individual Recode, women 15-49).
Variables used:
    v005 (cols 38-45) — sampling weight (divide by 1e6 per DHS convention)
    v012 (cols 66-67) — woman's age in years
    v730 (cols 4608-4609) — partner's age in years (current/cohabiting/most-recent;
                            missing for never-married women without recent partner)
    v821a (col 4776) — last partner younger/same/older (backup if v730 sparse)

Outputs:
    outputs/dhs_partner_age_raw.csv — per-respondent (woman_age, partner_age, gap, weight)
    outputs/dhs_partner_age_summary.csv — (mu, std, N, weighted_N) per woman age bin
    outputs/dhs_partner_age_distribution.png — boxplot/density of gap by woman age bin
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
OUT_DIR = EXP_DIR / "outputs"
OUT_DIR.mkdir(exist_ok=True)

DHS_DAT = REPO_ROOT / "references" / "SZIR51FL" / "szir51fl.dat"

# Fixed-width column positions (1-indexed in .dct, convert to 0-indexed slice)
COLS = {
    "v005": (37, 45),     # weight × 1e6
    "v012": (65, 67),     # respondent age
    "v730": (4607, 4609), # partner age
    "v821a": (4775, 4776),
}


def parse_dat(path: Path) -> pd.DataFrame:
    """Parse the DHS .dat file with the named variables only."""
    rows = []
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            line = line.rstrip("\n")
            row = {}
            for name, (lo, hi) in COLS.items():
                raw = line[lo:hi].strip()
                row[name] = pd.to_numeric(raw, errors="coerce")
            rows.append(row)
    df = pd.DataFrame(rows)
    df["weight"] = df["v005"] / 1e6
    return df


def gap_summary(df: pd.DataFrame, bins: list[tuple[int, int]]) -> pd.DataFrame:
    """Weighted mean/std of partner-age gap per woman age bin."""
    rows = []
    for lo, hi in bins:
        sub = df[(df["v012"] >= lo) & (df["v012"] < hi) & df["gap"].notna()]
        if len(sub) == 0:
            rows.append({"age_lo": lo, "age_hi": hi, "n": 0,
                         "weighted_n": 0, "mu_gap": np.nan, "sd_gap": np.nan})
            continue
        w = sub["weight"].values
        x = sub["gap"].values
        mu = np.average(x, weights=w)
        var = np.average((x - mu) ** 2, weights=w)
        rows.append({
            "age_lo": lo, "age_hi": hi,
            "n": len(sub),
            "weighted_n": float(w.sum()),
            "mu_gap": mu,
            "sd_gap": float(np.sqrt(var)),
        })
    return pd.DataFrame(rows)


def plot_distribution(df: pd.DataFrame, summary: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Box-style: gap distribution per 5-yr bin
    ax = axes[0]
    bins = [(15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 50)]
    data, labels = [], []
    for lo, hi in bins:
        sub = df[(df["v012"] >= lo) & (df["v012"] < hi) & df["gap"].notna()]
        if len(sub):
            data.append(sub["gap"].values)
            labels.append(f"{lo}-{hi-1}\n(n={len(sub)})")
    ax.boxplot(data, tick_labels=labels, showfliers=False)
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel("Woman's age (years)")
    ax.set_ylabel("Partner age − woman age (years)")
    ax.set_title("DHS Eswatini 2006-07 — partner-age gap by woman age\n(unweighted boxplot; whiskers=10/90 pct)")
    ax.grid(alpha=0.3)

    # Mean ± std overlay
    ax = axes[1]
    mids = [(r.age_lo + r.age_hi) / 2 for r in summary.itertuples()]
    mus = summary["mu_gap"]
    sds = summary["sd_gap"]
    ax.errorbar(mids, mus, yerr=sds, fmt="o-", capsize=5, color="C0",
                label="Weighted mean ± 1 SD")
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    for x, y, lo, hi in zip(mids, mus, summary["age_lo"], summary["age_hi"]):
        ax.annotate(f"μ={y:.1f}", (x, y), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=8)
    ax.set_xlabel("Woman's age (years)")
    ax.set_ylabel("Partner age − woman age (years)")
    ax.set_title("Weighted summary (input for stisim age_diff_pars)")
    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    print(f"Reading {DHS_DAT}")
    df = parse_dat(DHS_DAT)
    print(f"  {len(df)} records loaded")

    # Filter to plausible partner ages and compute gap
    df["gap"] = df["v730"] - df["v012"]
    df.loc[(df["v730"] < 12) | (df["v730"] > 80), "gap"] = np.nan
    df.loc[df["gap"].abs() > 40, "gap"] = np.nan
    print(f"  {df['gap'].notna().sum()} records with valid partner-age gap")

    # Save raw extract
    raw = df[["v012", "v730", "gap", "weight"]].rename(
        columns={"v012": "woman_age", "v730": "partner_age"})
    raw.to_csv(OUT_DIR / "dhs_partner_age_raw.csv", index=False)
    print(f"  Wrote {OUT_DIR / 'dhs_partner_age_raw.csv'}")

    # Summary by stisim age groups + finer 5-yr bins
    bins_stisim = [(15, 20), (20, 25), (25, 50)]   # teens / young / adult
    bins_5yr = [(15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 50)]

    summary_stisim = gap_summary(df, bins_stisim).assign(grouping="stisim")
    summary_5yr = gap_summary(df, bins_5yr).assign(grouping="5yr")
    summary = pd.concat([summary_stisim, summary_5yr], ignore_index=True)
    summary.to_csv(OUT_DIR / "dhs_partner_age_summary.csv", index=False)
    print(f"\nstisim age groups (mu_gap +/- sd_gap):")
    for _, r in summary_stisim.iterrows():
        print(f"  {r['age_lo']}-{r['age_hi']}: mu={r['mu_gap']:.2f}, sd={r['sd_gap']:.2f} "
              f"(n={r['n']:.0f}, weighted_n={r['weighted_n']:.0f})")

    plot_distribution(df, summary_5yr,
                      OUT_DIR / "dhs_partner_age_distribution.png")
    print(f"\nWrote {OUT_DIR / 'dhs_partner_age_distribution.png'}")


if __name__ == "__main__":
    main()
