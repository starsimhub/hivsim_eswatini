"""Exp 016 — counterfactual-shape comparison and the excess-death tables.

Kept separate from run.py because it answers a different question: not "what
does the fix do to the model" but "how much does the fix depend on the assumed
shape of the non-AIDS trend". The shape is not identifiable from this data —
two AIDS-free anchors, every year between them contaminated — so the honest
treatment is a sensitivity, not a single choice.

Usage (from repo root):
    python experiments/016_double_counted_mortality/analysis_methods.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
OUT, FIG = EXP_DIR / "outputs", EXP_DIR / "figures"
sys.path.insert(0, str(EXP_DIR))
from mortality_construction import build_hiv_deleted, deleted_fraction  # noqa: E402

METHODS = [
    ("loglinear", 1.0, "exponential (log-linear) — as run", "C0", "-"),
    ("linear",    1.0, "linear",                            "C2", "--"),
    ("power",     2.0, "delayed decline (power, p=2)",      "C4", "-."),
    ("sigmoid",   8.0, "S-curve (sigmoid, k=8)",            "C1", ":"),
]
AGES = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
BANDS = [(0, 15, "0-14"), (15, 50, "15-49"), (50, 200, "50+")]


def counterfactuals(all_cause):
    return {name: build_hiv_deleted(all_cause, method=name, par=par)
            for name, par, *_ in METHODS}


def plot_methods(all_cause, cfs, path, xmax=2035):
    """Observed all-cause vs every candidate counterfactual, per age and sex."""
    fig, axes = plt.subplots(3, 6, figsize=(19, 9), sharex=True)
    for ax, age in zip(axes.flat, AGES):
        for sex, mark in (("Female", "o"), ("Male", "^")):
            obs = all_cause[(all_cause.AgeStart == age) & (all_cause.Sex == sex)
                            & (all_cause.Time <= xmax)].sort_values("Time")
            ax.plot(obs.Time, obs.Value, color="k", lw=2.0, marker=mark, ms=3.5,
                    alpha=0.85 if sex == "Female" else 0.45,
                    label=f"{sex} observed" if age == 0 else None)
        for name, par, lab, col, ls in METHODS:
            g = cfs[name]
            g = g[(g.AgeStart == age) & (g.Sex == "Female") & (g.Time <= xmax)].sort_values("Time")
            ax.plot(g.Time, g.Value, color=col, ls=ls, lw=1.5,
                    label=lab if age == 0 else None)
        ax.set_yscale("log")
        ax.set_title(f"age {age}", fontsize=9)
        ax.grid(alpha=0.3, which="both")
        ax.tick_params(labelsize=7)
    for ax in axes.flat[len(AGES):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Year", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("Mortality rate (log)", fontsize=8)
    axes.flat[0].legend(fontsize=6, loc="lower left")
    fig.suptitle("Exp 016 — candidate non-AIDS counterfactuals (female shown as lines; "
                 "black = observed all-cause, both sexes).\nThe gap to the black line is "
                 "what each shape would attribute to AIDS.", fontsize=12)
    fig.tight_layout(); fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)


def excess_deaths(all_cause, cf, pop, years):
    """Excess deaths = (observed rate - counterfactual) x population, by age/sex/year."""
    d = deleted_fraction(all_cause, cf)[["Time", "Sex", "AgeStart", "deleted_rate"]].copy()
    d["bin"] = (d.AgeStart // 5) * 5                      # match population bins
    d = d.groupby(["Time", "Sex", "bin"], as_index=False).deleted_rate.mean()
    out = []
    for sex, key in (("Female", "f"), ("Male", "m")):
        for b in sorted(d["bin"].unique()):
            col = f"popagesex.n_alive_{key}_{b}_{b+5}"
            if col not in pop.columns:
                continue
            g = d[(d.Sex == sex) & (d["bin"] == b)].sort_values("Time")
            for y in years:
                rate = np.interp(y, g.Time.values, g.deleted_rate.values)
                n = pop.loc[pop.timevec == y, col].mean()
                out.append(dict(sex=sex, age=b, year=y, deaths=rate * n))
    return pd.DataFrame(out)


def main():
    all_cause = pd.read_csv(REPO_ROOT / "data" / "eswatini_deaths.csv")
    res = pd.read_parquet(OUT / "results.parquet")
    pop = res[(res.arm == "all_cause") & (res.pset == "default")]
    unaids = pd.read_csv(REPO_ROOT / "data" / "eswatini_hiv_calib.csv",
                         usecols=["time", "hiv.new_deaths"]).dropna()
    years = [1995, 2000, 2005, 2010, 2015, 2020]

    cfs = counterfactuals(all_cause)
    plot_methods(all_cause, cfs, FIG / "counterfactual_methods.png")

    # Excess deaths by age and sex for the method actually used
    ex = excess_deaths(all_cause, cfs["loglinear"], pop, years)
    tab = ex.pivot_table(index=["sex", "age"], columns="year", values="deaths").round(0)
    tab.to_csv(OUT / "excess_deaths_by_age_sex.csv")
    print("=== Implied AIDS deaths by age and sex (exponential/log-linear counterfactual) ===\n")
    print(tab.astype(int).to_string())

    print("\n=== Totals by age band ===\n")
    ex["band"] = pd.cut(ex.age, [-1, 14, 49, 200], labels=["0-14", "15-49", "50+"])
    band = ex.pivot_table(index="band", columns="year", values="deaths", aggfunc="sum", observed=True).round(0)
    print(band.astype(int).to_string())

    print("\n=== All methods vs UNAIDS ===\n")
    rows = {}
    for name, par, lab, *_ in METHODS:
        e = excess_deaths(all_cause, cfs[name], pop, years)
        rows[lab] = e.groupby("year").deaths.sum()
    comp = pd.DataFrame(rows).T
    u = unaids.set_index("time")["hiv.new_deaths"].reindex(years)
    comp.loc["UNAIDS Spectrum"] = u.values
    print(comp.round(0).astype(int).to_string())
    print("\nratio to UNAIDS:")
    print((comp.drop(index="UNAIDS Spectrum") / u.values).round(2).to_string())
    comp.to_csv(OUT / "method_comparison_vs_unaids.csv")
    print(f"\nWrote {OUT/'excess_deaths_by_age_sex.csv'} and {OUT/'method_comparison_vs_unaids.csv'}")


if __name__ == "__main__":
    main()
