"""Exp 030 analysis -- pick (k_m, k_f) from the scans, then audit the consequences.

  --phase scan     reads outputs/scan.parquet
                   -> outputs/scan_awareness.csv, outputs/scan_fit.csv
                   -> figures/testing_rate_scan.png
  --phase confirm  reads outputs/confirm.parquet
                   -> the full cascade + the ART reachability audit, every year
                   -> figures/{cascade_bars,art_reachability,residual_by_age}.png

Run from the repo ROOT.
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from cascade_analysis import (aggregate, conditionals, load_awareness_targets,
                              SHIMS_BINS, ART_BINS, PLOT_BINS, COLS, SEXLAB,
                              INK, MUTED)                       # noqa: E402

OUT = HERE / "outputs"
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------- scan phase

def scan():
    df = pd.read_parquet(OUT / "scan.parquet")
    df = df[df.timevec == 2021]
    long = aggregate(df, SHIMS_BINS, group_cols=("seed", "k_m", "k_f"))
    cond = conditionals(long, keys=("k_m", "k_f", "year", "sex", "bin"))
    aw = cond[(cond.measure == "aware_of_all_plhiv")
              & cond["bin"].isin(PLOT_BINS)].copy()

    tg = load_awareness_targets(REPO)
    aw = aw.merge(tg[["sex", "bin", "survey"]], on=["sex", "bin"], how="left")
    aw["err_pp"] = (aw.model - aw.survey) * 100
    aw.to_csv(OUT / "scan_awareness.csv", index=False)

    # Each sex is scanned with the OTHER held at 1.0, so select the right arm.
    rows = []
    for sex, kcol, other in (("m", "k_m", "k_f"), ("f", "k_f", "k_m")):
        arm = aw[(aw.sex == sex) & (np.isclose(aw[other], 1.0))]
        for k, g in arm.groupby(kcol):
            rows.append(dict(sex=sex, k=k, n_bins=len(g),
                             sse_pp2=float((g.err_pp ** 2).sum()),
                             rmse_pp=float(np.sqrt((g.err_pp ** 2).mean())),
                             mean_err_pp=float(g.err_pp.mean()),
                             max_abs_err_pp=float(g.err_pp.abs().max())))
    fit = pd.DataFrame(rows).sort_values(["sex", "k"])
    fit.to_csv(OUT / "scan_fit.csv", index=False)

    best = {s: float(g.loc[g.sse_pp2.idxmin(), "k"])
            for s, g in fit.groupby("sex")}
    _fig_scan(aw, fit, best)

    print(fit.round(2).to_string(index=False))
    print(f"\nbest by SSE:  k_m = {best['m']},  k_f = {best['f']}")
    print("\nregression check -- k=(1.0, 1.0) awareness vs exp 029 (20 seeds):")
    base = aw[np.isclose(aw.k_m, 1.0) & np.isclose(aw.k_f, 1.0)]
    e029 = pd.read_csv(REPO / "experiments" / "029_cascade_audit" / "outputs"
                       / "cascade_vs_shims3.csv")
    e029 = e029[(e029.measure == "aware_of_all_plhiv")
                & e029["bin"].isin(PLOT_BINS)][["sex", "bin", "model"]]
    chk = base.merge(e029, on=["sex", "bin"], suffixes=("_030", "_029"))
    chk["delta_pp"] = (chk.model_030 - chk.model_029) * 100
    print(chk[["sex", "bin", "model_030", "model_029", "delta_pp"]]
          .round(3).to_string(index=False))
    print(f"max |delta| = {chk.delta_pp.abs().max():.2f} pp "
          f"(seed noise expected: 030 uses 3 seeds, 029 used 20)")
    return best


def _fig_scan(aw, fit, best):
    """Awareness against the testing multiplier, one panel per sex.

    Age bands share a panel because they share an x-axis and a unit, and the
    question is whether ONE k can satisfy all four -- which is exactly what
    overlaying them shows.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    styles = dict(zip(PLOT_BINS, ["-o", "-s", "-^", "-D"]))
    shades = dict(zip(PLOT_BINS, [0.35, 0.58, 0.8, 1.0]))

    for ax, (sex, kcol, other) in zip(axes, (("m", "k_m", "k_f"),
                                             ("f", "k_f", "k_m"))):
        arm = aw[(aw.sex == sex) & (np.isclose(aw[other], 1.0))]
        for b in PLOT_BINS:
            g = arm[arm["bin"] == b].sort_values(kcol)
            if not len(g):
                continue
            c = plt.matplotlib.colors.to_rgba(COLS[sex], shades[b])
            ax.plot(g[kcol], g.model, styles[b], color=c, lw=1.8, ms=6,
                    label=b, zorder=3)
            tgt = g.survey.iloc[0]
            ax.axhline(tgt, color=c, lw=1.0, ls=(0, (3, 3)), zorder=1)
            ax.text(ax.get_xlim()[1], tgt, f" {b}", fontsize=7.5, color=c,
                    va="center", ha="left")
        k = best[sex]
        ax.axvline(k, color=INK, lw=1.3, ls=(0, (4, 3)), zorder=2)
        ax.text(k, 0.455, f" best SSE\n k={k:g}", fontsize=8.5, color=INK,
                ha="left", va="bottom")
        ax.set_title(f"{SEXLAB[sex]} — {kcol} (other sex held at 1.0)",
                     fontsize=11, color=INK)
        ax.set_xlabel("Testing-rate multiplier on the general-population ramp",
                      fontsize=9.5, color=MUTED)
        ax.grid(color="#ececec", lw=0.8)
        ax.set_axisbelow(True)
        ax.set_ylim(0.45, 1.03)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(frameon=False, fontsize=8.5, loc="lower right",
                  title="age band", title_fontsize=8.5)

    axes[0].set_ylabel("Aware of status | living with HIV", fontsize=9.5,
                       color=MUTED)
    fig.suptitle("Can one testing rate per sex reproduce all four age bands?",
                 fontsize=13, color=INK, y=1.0)
    fig.text(0.5, -0.06, "Solid = model (3 seeds, N=20,000, 2021). Dashed "
             "horizontal = SHIMS3 2021 target for the band of the same shade. "
             "If the solid curves cross their targets at different k, no single "
             "rate fits the age profile.", ha="center", fontsize=8.5, color=MUTED)
    fig.tight_layout()
    fig.savefig(FIG / "testing_rate_scan.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------- confirm phase

def confirm():
    df = pd.read_parquet(OUT / "confirm.parquet")
    km, kf = df.k_m.iloc[0], df.k_f.iloc[0]

    long = aggregate(df, SHIMS_BINS, group_cols=("seed",))
    cond = conditionals(long, keys=("year", "sex", "bin"))
    tg = load_awareness_targets(REPO)
    comp = (cond[(cond.year == 2021) & cond["bin"].isin(PLOT_BINS)]
            .merge(tg, on=["sex", "bin", "measure"], how="inner"))
    comp["diff_pp"] = (comp.model - comp.survey) * 100
    comp.to_csv(OUT / "confirm_awareness.csv", index=False)

    # --- the primary diagnostic: is the ART target reachable, EVERY year? ---
    art_cond = conditionals(aggregate(df, ART_BINS, group_cols=("seed",)),
                            keys=("year", "sex", "bin"))
    piv = art_cond[art_cond.measure.isin(["aware_of_all_plhiv",
                                          "on_art_of_all_plhiv"])]
    piv = piv.pivot_table(index=["year", "sex", "bin"], columns="measure",
                          values="model").reset_index()
    tgt = pd.read_csv(REPO / "data" / "art_coverage.csv")
    # Gender 0 = FEMALE here (stisim's convention, per
    # art_coverage_construction.py). This is the OPPOSITE of
    # calibration_data/art_coverage_by_age_sex.csv. Reversing it silently swaps
    # the sexes and still produces a plausible figure.
    tgt["sex"] = tgt.Gender.map({0: "f", 1: "m"})
    tgt = tgt.rename(columns={"Year": "year", "AgeBin": "bin", "p_art": "target"})

    # The table is tabulated at a handful of years (2004/2011/2016/2021) and
    # stisim INTERPOLATES between them, so merging on the tabulated years only
    # would check 32 strata-years instead of every year the model ran -- and the
    # binding constraint could easily fall in an interpolated year. Rebuild the
    # annual series per stratum the way parse_coverage does at smoothness=0
    # (linear), held flat outside the tabulated range.
    years = np.sort(piv.year.unique())
    filled = []
    for (sex, b), g in tgt.groupby(["sex", "bin"]):
        g = g.sort_values("year")
        filled.append(pd.DataFrame({
            "year": years, "sex": sex, "bin": b,
            "target": np.interp(years, g.year.values, g.target.values),
            "interpolated": ~np.isin(years, g.year.values)}))
    tgt_annual = pd.concat(filled, ignore_index=True)

    reach = piv.merge(tgt_annual, on=["year", "sex", "bin"], how="inner")
    reach["shortfall_pp"] = (reach.on_art_of_all_plhiv - reach.target) * 100
    reach["ceiling_slack_pp"] = (reach.aware_of_all_plhiv - reach.target) * 100
    reach.to_csv(OUT / "art_reachability.csv", index=False)

    _fig_reach(reach, km, kf)
    _fig_residual(comp, km, kf)

    print(f"confirm: k_m={km}, k_f={kf}, {df.seed.nunique()} seeds\n")
    print("awareness vs SHIMS3 2021:")
    print(comp[["sex", "bin", "model", "survey", "diff_pp"]]
          .sort_values(["sex", "bin"]).round(3).to_string(index=False))
    worst = reach.loc[reach.shortfall_pp.idxmin()]
    print(f"\nART reachability, {int(reach.year.min())}-{int(reach.year.max())}:")
    print(f"  worst shortfall: {worst.shortfall_pp:+.2f} pp "
          f"({worst.sex} {worst['bin']} in {int(worst.year)})")
    print(f"  strata-years short by >1 pp: "
          f"{int((reach.shortfall_pp < -1).sum())} of {len(reach)}")
    tight = reach.loc[reach.ceiling_slack_pp.idxmin()]
    print(f"  tightest ceiling: {tight.ceiling_slack_pp:+.2f} pp slack "
          f"({tight.sex} {tight['bin']} in {int(tight.year)})")


def _fig_reach(reach, km, kf):
    """ART coverage achieved vs target, and the awareness ceiling, every year."""
    bins = list(ART_BINS)
    fig, axes = plt.subplots(2, len(bins), figsize=(15, 7.2), sharey=True,
                             sharex=True)
    for row, sex in enumerate(("f", "m")):
        for col, b in enumerate(bins):
            ax = axes[row, col]
            g = reach[(reach.sex == sex) & (reach["bin"] == b)].sort_values("year")
            ax.fill_between(g.year, g.target, g.aware_of_all_plhiv,
                            where=g.aware_of_all_plhiv >= g.target,
                            color="#cfe3cf", lw=0, zorder=1, label="headroom")
            ax.fill_between(g.year, g.target, g.aware_of_all_plhiv,
                            where=g.aware_of_all_plhiv < g.target,
                            color="#f2c7c7", lw=0, zorder=1,
                            label="target above ceiling")
            ax.plot(g.year, g.aware_of_all_plhiv, "-", color=MUTED, lw=1.4,
                    zorder=3, label="awareness (ceiling)")
            ax.plot(g.year, g.target, "--", color=INK, lw=1.4, zorder=4,
                    label="ART target")
            ax.plot(g.year, g.on_art_of_all_plhiv, "-", color=COLS[sex], lw=2,
                    zorder=5, label="ART achieved")
            ax.set_ylim(0, 1.05)
            ax.grid(color="#ececec", lw=0.8)
            ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if row == 0:
                ax.set_title(b, fontsize=10.5, color=INK)
            if row == 1:
                ax.set_xlabel("Year", fontsize=9, color=MUTED)
        axes[row, 0].set_ylabel(f"{SEXLAB[sex]}\nshare of PLHIV", fontsize=9.5,
                                color=COLS[sex])
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.04))
    fig.suptitle(f"Is the ART coverage target reachable once testing is slowed? "
                 f"(k_m={km:g}, k_f={kf:g})", fontsize=13, color=INK, y=1.0)
    fig.text(0.5, -0.085, "Awareness is a hard ceiling: stisim fills the ART "
             "target only from already-diagnosed agents. Red shading is where "
             "the target sits above the ceiling and therefore cannot fill.",
             ha="center", fontsize=8.5, color=MUTED)
    fig.tight_layout()
    fig.savefig(FIG / "art_reachability.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def _fig_residual(comp, km, kf):
    """Awareness residual by age band -- the ramp-SHAPE flag.

    Systematic ordering by age band means the assumed scale-up shape is wrong,
    not just its level. CONFOUNDED with the model's age-at-infection
    distribution, which is unconstrained here (the age-banded SHIMS incidence
    rows are fit=False for false-recency bias). Report, do not diagnose.
    """
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    x = np.arange(len(PLOT_BINS))
    for sex in ("f", "m"):
        g = comp[comp.sex == sex].set_index("bin").reindex(PLOT_BINS)
        ax.plot(x, g.diff_pp, "-o", color=COLS[sex], lw=2, ms=8, mec="white",
                mew=1.2, label=SEXLAB[sex])
    ax.axhline(0, color=INK, lw=1.2, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(PLOT_BINS)
    ax.set_xlabel("Age band", fontsize=9.5, color=MUTED)
    ax.set_ylabel("Model - SHIMS3 (percentage points)", fontsize=9.5, color=MUTED)
    ax.set_title(f"Awareness residual by age (k_m={km:g}, k_f={kf:g})\n"
                 f"Systematic ordering by age would flag the ramp SHAPE, "
                 f"not its level", fontsize=11, color=INK)
    ax.grid(axis="y", color="#ececec", lw=0.8)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=9.5)
    fig.tight_layout()
    fig.savefig(FIG / "residual_by_age.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["scan", "confirm"], default="scan")
    a = ap.parse_args()
    scan() if a.phase == "scan" else confirm()
