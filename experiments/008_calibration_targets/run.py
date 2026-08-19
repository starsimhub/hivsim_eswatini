"""
Experiment 008 — Calibration targets.

Load PHIA prevalence (age × sex × year) and UNAIDS HIV deaths into a single
tidy dataframe, audit uncertainty coverage, and plot the target set with CIs.
Output: outputs/calibration_targets.csv + outputs/calibration_targets.png.

See README.md for the full plan and success criteria.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]  # hivsim_eswatini/
OUT_DIR = EXP_DIR / "outputs"
OUT_DIR.mkdir(exist_ok=True)

CALIB_DATA = REPO_ROOT / "calibration_data"
DATA = REPO_ROOT / "data"

# Tidy schema for every target row
TARGET_COLS = [
    "quantity", "year", "sex", "age_low", "age_high",
    "value", "lower", "upper", "source", "notes",
]

# UNAIDS Spectrum estimates typically carry ~15% relative uncertainty on
# annual AIDS deaths. We apply this as a placeholder until we source the
# bounds directly from Spectrum/AIDSinfo.
UNAIDS_DEATHS_REL_UNCERTAINTY = 0.15


def _parse_age_bin(bin_str: str) -> tuple[int, int]:
    """'[15:20)' → (15, 20)."""
    m = re.match(r"\[(\d+):(\d+)\)", bin_str)
    if not m:
        raise ValueError(f"Unexpected age bin format: {bin_str!r}")
    return int(m.group(1)), int(m.group(2))


def load_phia_prevalence() -> pd.DataFrame:
    """PHIA prevalence by age × sex × year (2007, 2011, 2016)."""
    df = pd.read_csv(CALIB_DATA / "prevalence_by_age_sex.csv")
    sex_map = {0: "M", 1: "F"}
    age_pairs = df["AgeBin"].map(_parse_age_bin)
    out = pd.DataFrame({
        "quantity": "prevalence",
        "year": df["Year"].astype(int),
        "sex": df["Gender"].map(sex_map),
        "age_low": [a[0] for a in age_pairs],
        "age_high": [a[1] for a in age_pairs],
        "value": df["NationalPrevalence"],
        "lower": df["lb"],
        "upper": df["ub"],
        "source": "PHIA",
        "notes": "n=" + df["Count"].round().astype("Int64").astype(str),
    })
    return out[TARGET_COLS]


def load_unaids_deaths() -> pd.DataFrame:
    """UNAIDS annual AIDS deaths (point estimates) from eswatini_hiv_calib.csv.

    No uncertainty in source file — applying ±15% relative as a placeholder
    consistent with typical Spectrum estimates. Replace with sourced bounds
    in a follow-up.
    """
    df = pd.read_csv(DATA / "eswatini_hiv_calib.csv", usecols=["time", "hiv.new_deaths"])
    df = df.dropna(subset=["hiv.new_deaths"]).copy()
    out = pd.DataFrame({
        "quantity": "aids_deaths",
        "year": df["time"].astype(int),
        "sex": "both",
        "age_low": pd.NA,
        "age_high": pd.NA,
        "value": df["hiv.new_deaths"],
        "lower": df["hiv.new_deaths"] * (1 - UNAIDS_DEATHS_REL_UNCERTAINTY),
        "upper": df["hiv.new_deaths"] * (1 + UNAIDS_DEATHS_REL_UNCERTAINTY),
        "source": "UNAIDS Spectrum (point); uncertainty = ±15% placeholder",
        "notes": "uncertainty placeholder — source bounds from AIDSinfo as follow-up",
    })
    return out[TARGET_COLS]


def audit(targets: pd.DataFrame) -> pd.DataFrame:
    """Per-row audit of uncertainty coverage."""
    return targets.assign(
        has_lower=targets["lower"].notna(),
        has_upper=targets["upper"].notna(),
        sourced_uncertainty=~targets["notes"].str.contains("placeholder", na=False),
    )


def plot_targets(targets: pd.DataFrame, path: Path) -> None:
    """Two-panel figure: PHIA prevalence (faceted by sex) and UNAIDS deaths."""
    prev = targets[targets["quantity"] == "prevalence"].copy()
    deaths = targets[targets["quantity"] == "aids_deaths"].copy()

    prev["age_mid"] = (prev["age_low"] + prev["age_high"]) / 2
    years = sorted(prev["year"].unique())
    sexes = ["M", "F"]
    colors = plt.cm.viridis([i / max(len(years) - 1, 1) for i in range(len(years))])

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # PHIA prevalence — one panel per sex
    for ax, sex in zip(axes[:2], sexes):
        sub = prev[prev["sex"] == sex]
        for color, year in zip(colors, years):
            yr = sub[sub["year"] == year].sort_values("age_mid")
            if yr.empty:
                continue
            err = [yr["value"] - yr["lower"], yr["upper"] - yr["value"]]
            ax.errorbar(yr["age_mid"], yr["value"], yerr=err, marker="o",
                        capsize=3, label=str(year), color=color)
        ax.set_title(f"PHIA prevalence — {'Male' if sex == 'M' else 'Female'}")
        ax.set_xlabel("Age (mid-bin)")
        ax.set_ylabel("HIV prevalence")
        ax.set_ylim(0, max(prev["upper"].max() * 1.05, 0.6))
        ax.grid(alpha=0.3)
        ax.legend(title="PHIA year", fontsize=8)

    # UNAIDS deaths — line + shaded ±15% band
    ax = axes[2]
    deaths_sorted = deaths.sort_values("year")
    ax.plot(deaths_sorted["year"], deaths_sorted["value"], marker=".",
            color="C3", label="UNAIDS point")
    ax.fill_between(deaths_sorted["year"], deaths_sorted["lower"],
                    deaths_sorted["upper"], color="C3", alpha=0.2,
                    label="±15% (placeholder)")
    ax.set_title("UNAIDS HIV deaths (annual)")
    ax.set_xlabel("Year")
    ax.set_ylabel("New AIDS deaths")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle("Eswatini calibration targets — experiment 008", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parts = [
        load_phia_prevalence(),
        load_unaids_deaths(),
    ]
    targets = pd.concat(parts, ignore_index=True)

    audit_df = audit(targets)
    n = len(audit_df)
    n_with_unc = int(audit_df["has_lower"].sum())
    n_sourced = int(audit_df["sourced_uncertainty"].sum())
    print(f"Loaded {n} target rows.")
    print(f"  {n_with_unc}/{n} rows have lower/upper bounds populated.")
    print(f"  {n_sourced}/{n} rows have uncertainty sourced (not placeholder).")
    by_q = targets.groupby("quantity").size()
    for q, count in by_q.items():
        print(f"  {q}: {count} rows")

    targets_path = OUT_DIR / "calibration_targets.csv"
    audit_path = OUT_DIR / "calibration_targets_audit.csv"
    fig_path = OUT_DIR / "calibration_targets.png"

    targets.to_csv(targets_path, index=False)
    audit_df.to_csv(audit_path, index=False)
    plot_targets(targets, fig_path)

    print(f"\nWrote {targets_path}")
    print(f"Wrote {audit_path}")
    print(f"Wrote {fig_path}")


if __name__ == "__main__":
    main()
