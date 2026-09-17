"""Exp 029 analysis -- model cascade vs SHIMS3, and the ART coverage ceiling.

Produces:
  outputs/cascade_vs_shims3.csv   model vs survey, every step x bin x sex, 2021
  outputs/cascade_timeseries.csv  the five cascade quantities, all years/bins
  outputs/art_ceiling.csv         achievable ART coverage per stratum vs 0.95
  figures/cascade_95s.png         the three conditional steps, model vs SHIMS3
  figures/cascade_bars.png        the conventional cascade, share of all PLHIV
  figures/coverage_vs_target.png  model ART coverage vs its input, all years
  figures/art_ceiling.png         awareness as the ceiling on ART coverage
  figures/prevalence_fit_vs_phia.png   the standard figure (workflow requirement)

Run from the repo ROOT.
"""

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

OUT = HERE / "outputs"
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# House colours, carried from standard_figures.py so 029 reads as the same
# project as 019-024. Validated: adjacent-pair CVD dE 21.8 (protan), 28.8 normal.
COLS = {"f": "#c0392b", "m": "#2c6fbb"}
SEXLAB = {"f": "Women", "m": "Men"}
INK, MUTED = "#222222", "#6b6b6b"

# The analyzer's 5-year bands, summed into whichever grouping is needed.
BANDS = [(a, a + 5) for a in range(15, 80, 5)] + [(80, 200)]
SHIMS_BINS = {"15-24": (15, 25), "25-34": (25, 35), "35-49": (35, 50),
              "50+": (50, 200), "15-49": (15, 50), "15+": (15, 200)}
PLOT_BINS = ["15-24", "25-34", "35-49", "50+"]
# data/art_coverage.csv drives the model on these. [45,100) is read as 45+ --
# the model has no agents above 100, so the two are the same set.
ART_BINS = {"[15,25)": (15, 25), "[25,35)": (25, 35),
            "[35,45)": (35, 45), "[45,100)": (45, 200)}

STATES = ("n_infected", "n_diagnosed", "n_on_art", "n_effective_art")


def load():
    return pd.read_parquet(OUT / "results.parquet")


def aggregate(df, bins):
    """Sum 5-year bands into `bins`; returns long (year, seed, sex, bin, state)."""
    rows = []
    for name, (lo_b, hi_b) in bins.items():
        members = [(lo, hi) for lo, hi in BANDS if lo >= lo_b and hi <= hi_b]
        if not members:
            raise ValueError(f"no 5-year bands inside {name}")
        for sex in ("f", "m"):
            acc = {}
            for state in STATES:
                cols = [f"cascadeage.{state}_{sex}_{lo}_{hi}" for lo, hi in members]
                acc[state] = df[cols].sum(axis=1)
            sub = pd.DataFrame(acc)
            sub["year"] = df["timevec"].values
            sub["seed"] = df["seed"].values
            sub["sex"] = sex
            sub["bin"] = name
            rows.append(sub)
    return pd.concat(rows, ignore_index=True)


def conditionals(long):
    """Pooled point estimate + across-seed spread for the five cascade measures.

    The point estimate pools numerator and denominator across seeds before
    dividing. Averaging per-seed ratios instead would bias small strata, where
    a seed with few PLHIV carries the same weight as one with many -- and in
    men 15-24 that is a real difference, not a hypothetical.
    """
    defs = {  # measure: (numerator, denominator)
        "aware_of_all_plhiv":  ("n_diagnosed", "n_infected"),
        "on_art_given_aware":  ("n_on_art", "n_diagnosed"),
        "vls_given_art":       ("n_effective_art", "n_on_art"),
        "on_art_of_all_plhiv": ("n_on_art", "n_infected"),
        "vls_of_all_plhiv":    ("n_effective_art", "n_infected"),
    }
    pooled = long.groupby(["year", "sex", "bin"], as_index=False)[list(STATES)].sum()
    out = []
    for measure, (num, den) in defs.items():
        p = pooled.copy()
        p["measure"] = measure
        p["model"] = p[num] / p[den].replace(0, np.nan)
        per = long.copy()
        per["r"] = per[num] / per[den].replace(0, np.nan)
        sd = (per.groupby(["year", "sex", "bin"])["r"]
                 .agg(model_sd="std", n_seeds="count", model_seedmean="mean")
                 .reset_index())
        merged = p.merge(sd, on=["year", "sex", "bin"], how="left")
        merged = merged.rename(columns={num: "numerator", den: "denominator"})
        out.append(merged[["year", "sex", "bin", "measure", "model", "model_sd",
                           "model_seedmean", "n_seeds", "numerator", "denominator"]])
    return pd.concat(out, ignore_index=True)


def main():
    df = load()
    n_seeds = df.seed.nunique()
    long = aggregate(df, SHIMS_BINS)
    cond = conditionals(long)
    cond.to_csv(OUT / "cascade_timeseries.csv", index=False)

    # ---- model vs SHIMS3 2021 -------------------------------------------
    shims = pd.read_csv(REPO / "data" / "eswatini_cascade_95s.csv")
    # The file's own `survey` column is the survey NAME ("SHIMS3 2021"); rename
    # it out of the way before `value` takes that slot, or the merge silently
    # produces two columns called `survey` and every later `.survey` breaks.
    shims = (shims[shims.sex.isin(["f", "m"])]
             .rename(columns={"survey": "survey_name", "age": "bin",
                              "value": "survey"}))
    comp = (cond[cond.year == 2021]
            .merge(shims[["sex", "bin", "measure", "survey", "n_unweighted",
                          "small_denominator"]],
                   on=["sex", "bin", "measure"], how="inner"))
    comp["diff_pp"] = (comp.model - comp.survey) * 100
    comp.to_csv(OUT / "cascade_vs_shims3.csv", index=False)

    # ---- the ceiling ----------------------------------------------------
    # Awareness is a hard cap on ART coverage: stisim fills a stratified ART
    # target only from `diagnosed & ~on_art` in the stratum. So whatever target
    # 026 set, this is what the stratum could actually deliver.
    art_cond = conditionals(aggregate(df, ART_BINS))
    ceil = art_cond[art_cond.measure.isin(["aware_of_all_plhiv",
                                           "on_art_of_all_plhiv"])
                    & art_cond.year.isin([2021, 2030])]
    ceil = ceil.pivot_table(index=["year", "sex", "bin"], columns="measure",
                            values="model").reset_index()
    ceil["headroom_to_95"] = 0.95 - ceil["aware_of_all_plhiv"]
    ceil["art_95_reachable"] = ceil["aware_of_all_plhiv"] >= 0.95
    ceil.to_csv(OUT / "art_ceiling.csv", index=False)

    _fig_cascade(comp, n_seeds)
    _fig_cascade_bars(comp, n_seeds)
    _fig_ceiling(ceil)
    cov = _fig_coverage_vs_target(art_cond)
    cov.to_csv(OUT / "coverage_vs_target.csv", index=False)
    _fig_standard(df)

    show = comp[comp.measure.isin(["aware_of_all_plhiv", "on_art_given_aware",
                                   "vls_given_art"])]
    print(show.pivot_table(index=["measure", "bin"], columns="sex",
                           values=["model", "survey", "diff_pp"])
              .round(3).to_string())
    print("\nART ceiling (2030):")
    print(ceil[ceil.year == 2030].round(3).to_string(index=False))


def _fig_cascade(comp, n_seeds):
    """Three conditional steps, one panel each, shared 0-1 axis.

    Small multiples rather than one crowded panel: the three steps share a unit
    (a proportion of the previous step) but not a meaning, and stacking them on
    one axis would invite reading a trend across steps that does not exist.
    """
    steps = [("aware_of_all_plhiv", "1st 95\nAware of status | living with HIV"),
             ("on_art_given_aware", "2nd 95\nOn ART | aware"),
             ("vls_given_art", "3rd 95\nSuppressed | on ART")]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), sharey=True)
    x = np.arange(len(PLOT_BINS))

    for ax, (measure, title) in zip(axes, steps):
        sub = comp[comp.measure == measure]
        for sex in ("f", "m"):
            s = sub[sub.sex == sex].set_index("bin").reindex(PLOT_BINS)
            ax.fill_between(x, s.model - s.model_sd, s.model + s.model_sd,
                            color=COLS[sex], alpha=0.15, lw=0, zorder=1)
            ax.plot(x, s.model, "-o", color=COLS[sex], lw=2, ms=8, zorder=3,
                    mec="white", mew=1.2)
            ax.plot(x, s.survey, "s", ms=9, mfc="none", mew=2.0,
                    color=COLS[sex], zorder=4)
            # Women's gap labelled to the right of the stem, men's to the left.
            # Without the split the two sexes' labels overlap wherever the gaps
            # are similar -- which is most of the 15-24 column.
            dx, ha = (0.10, "left") if sex == "f" else (-0.10, "right")
            last = len(PLOT_BINS) - 1
            for xi, (m, v, small) in enumerate(zip(s.model, s.survey,
                                                   s.small_denominator)):
                # At the two end columns an outward label lands in the gutter
                # between panels, where it reads as belonging to the neighbour.
                # Flip it inward instead.
                ldx, lha = dx, ha
                if (xi == last and dx > 0) or (xi == 0 and dx < 0):
                    ldx, lha = -dx, ("right" if ha == "left" else "left")
                if np.isfinite(m) and np.isfinite(v) and abs(m - v) > 0.03:
                    ax.annotate("", xy=(xi, v), xytext=(xi, m),
                                arrowprops=dict(arrowstyle="-", ls=(0, (2, 2)),
                                                color=COLS[sex], lw=1.1), zorder=2)
                    ax.text(xi + ldx, (m + v) / 2, f"{(m - v) * 100:+.0f}pp",
                            fontsize=8.5, color=MUTED, va="center", ha=lha)
                if small:
                    ax.text(xi, v - 0.028, "small n", fontsize=7, color=MUTED,
                            ha="center", style="italic")
        ax.axhline(0.95, color=MUTED, lw=1.1, ls=(0, (1, 2)), zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels(PLOT_BINS)
        ax.set_title(title, fontsize=10.5, color=INK)
        ax.set_xlabel("Age (years)", fontsize=9.5, color=MUTED)
        # Data spans 0.72-1.00; a floor at 0.45 left half the canvas empty and
        # flattened the one gap the figure exists to show.
        ax.set_ylim(0.68, 1.035)
        ax.grid(axis="y", color="#e8e8e8", lw=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    axes[0].set_ylabel("Proportion", fontsize=9.5, color=MUTED)
    # Left panel, left edge: the only corner with no marks near the 0.95 line.
    axes[0].text(-0.42, 0.957, "95% target", fontsize=8, color=MUTED, ha="left")

    handles = [plt.Line2D([], [], color=COLS["f"], lw=2, marker="o", ms=8,
                          mec="white", label="Model - women"),
               plt.Line2D([], [], color=COLS["m"], lw=2, marker="o", ms=8,
                          mec="white", label="Model - men"),
               plt.Line2D([], [], color=COLS["f"], lw=0, marker="s", ms=9,
                          mfc="none", mew=2.0, label="SHIMS3 2021 - women"),
               plt.Line2D([], [], color=COLS["m"], lw=0, marker="s", ms=9,
                          mfc="none", mew=2.0, label="SHIMS3 2021 - men")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("The conditional 95-95-95 cascade in 2021: model against SHIMS3, "
                 "by age and sex", fontsize=12.5, color=INK, y=1.0)
    fig.text(0.5, -0.10, f"Model: exp 029, 024 row 868, N=20,000, {n_seeds} seeds; "
             "band is +/-1 SD across seeds. Survey: SHIMS3 2021 Table 9.1.B "
             "(no published CI). Labels show model - survey in percentage points.",
             ha="center", fontsize=8.5, color=MUTED)
    fig.tight_layout()
    fig.savefig(FIG / "cascade_95s.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def _fig_ceiling(ceil):
    """Awareness as the cap on ART coverage, per ART-target stratum."""
    c = ceil[ceil.year == 2030]
    bins = list(ART_BINS)
    x = np.arange(len(bins))
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    for sex in ("f", "m"):
        s = c[c.sex == sex].set_index("bin").reindex(bins)
        ax.plot(x, s["aware_of_all_plhiv"], "-o", color=COLS[sex], lw=2, ms=8,
                mec="white", mew=1.2, label=f"{SEXLAB[sex]} - aware of status")
        ax.plot(x, s["on_art_of_all_plhiv"], "--^", color=COLS[sex], lw=1.6,
                ms=7, mfc="none", label=f"{SEXLAB[sex]} - on ART")
    ax.axhline(0.95, color="#444444", lw=1.4, ls=(0, (4, 3)), zorder=0)
    # Left edge: at [15,25) nothing sits above 0.95, whereas on the right the
    # label lands on top of the women's on-ART line.
    ax.text(-0.08, 0.957, "exp 026 art_95 target", fontsize=8.5, color=INK,
            ha="left")
    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_xlabel("ART coverage stratum (data/art_coverage.csv bins)",
                  fontsize=9.5, color=MUTED)
    ax.set_ylabel("Proportion of people living with HIV", fontsize=9.5, color=MUTED)
    ax.set_title("Awareness caps ART coverage: the model cannot treat whom it has "
                 "not diagnosed\nBaseline, 2030", fontsize=11.5, color=INK)
    ax.grid(axis="y", color="#e8e8e8", lw=0.8)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "art_ceiling.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def _fig_cascade_bars(comp, n_seeds):
    """The cascade in its conventional form: each stage as a share of all PLHIV.

    This is the same data as `cascade_95s.png` but on the survey's *overall*
    basis (SHIMS3 Table 9.1.A) rather than the conditional one (Table 9.1.B).
    Both are worth having: the conditional view isolates which step is wrong,
    the unconditional view is the familiar descending cascade and is what the
    95-95-95 targets are usually quoted against.

    Model is a solid bar, survey a hollow one -- fill, not hue, separates them,
    so sex keeps the project's red/blue and identity is never colour-alone.
    """
    stages = [("aware_of_all_plhiv", "Aware\nof status"),
              ("on_art_of_all_plhiv", "On\nART"),
              ("vls_of_all_plhiv", "Virally\nsuppressed")]
    fig, axes = plt.subplots(2, len(PLOT_BINS), figsize=(14.5, 7.4),
                             sharey=True, sharex=True)
    xs = np.arange(len(stages))
    w = 0.34

    for row, sex in enumerate(("f", "m")):
        for col, b in enumerate(PLOT_BINS):
            ax = axes[row, col]
            s = (comp[(comp.sex == sex) & (comp["bin"] == b)]
                 .set_index("measure").reindex([m for m, _ in stages]))
            ax.bar(xs - w / 2, s.model, width=w, color=COLS[sex], zorder=3)
            ax.bar(xs + w / 2, s.survey, width=w, facecolor="none",
                   edgecolor=COLS[sex], lw=1.8, zorder=3)
            # Seed spread on the model bar only -- the survey has no published CI.
            ax.errorbar(xs - w / 2, s.model, yerr=s.model_sd, fmt="none",
                        ecolor="#444444", elinewidth=1.1, capsize=3, zorder=4)
            for xi, (m, v) in enumerate(zip(s.model, s.survey)):
                if np.isfinite(m) and np.isfinite(v) and abs(m - v) >= 0.03:
                    ax.text(xi, max(m, v) + 0.035, f"{(m - v) * 100:+.0f}pp",
                            ha="center", fontsize=8.5, color=MUTED)
            ax.axhline(0.95, color=MUTED, lw=1.0, ls=(0, (1, 2)), zorder=1)
            ax.set_ylim(0, 1.16)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 0.95])
            ax.set_yticklabels(["0", "25%", "50%", "75%", "95%"])
            ax.set_xticks(xs)
            ax.set_xticklabels([lbl for _, lbl in stages], fontsize=8.5)
            ax.grid(axis="y", color="#ececec", lw=0.8)
            ax.set_axisbelow(True)
            for sp in ("top", "right", "left"):
                ax.spines[sp].set_visible(False)
            if row == 0:
                ax.set_title(f"{b}", fontsize=11, color=INK)
        axes[row, 0].set_ylabel("Share of all people\nliving with HIV",
                                fontsize=9, color=MUTED)
        # Sex as a row label in its own colour; the ylabel stays neutral ink so
        # the two are not competing for the same slot.
        axes[row, 0].text(-0.30, 0.5, SEXLAB[sex], transform=axes[row, 0].transAxes,
                          rotation=90, va="center", ha="center", fontsize=12,
                          color=COLS[sex])

    # Neutral swatches: fill vs outline is what separates model from survey, and
    # a red legend key above a blue row would imply hue carries that meaning.
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor="#8a8a8a", edgecolor="none",
                             label="Model"),
               plt.Rectangle((0, 0), 1, 1, facecolor="none", edgecolor="#8a8a8a",
                             lw=1.8, label="SHIMS3 2021")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=10, bbox_to_anchor=(0.5, -0.035))
    fig.suptitle("HIV treatment cascade, 2021: model against SHIMS3, "
                 "by age and sex", fontsize=13, color=INK, y=1.0)
    fig.text(0.5, -0.075, f"Each bar is a share of ALL people living with HIV in "
             f"that age/sex group (SHIMS3 2021 Table 9.1.A basis). Model: exp 029, "
             f"024 row 868, N=20,000, {n_seeds} seeds; error bars +/-1 SD across "
             f"seeds. Survey has no published CI. Gaps under 3 pp are unlabelled.",
             ha="center", fontsize=8.5, color=MUTED)
    fig.tight_layout()
    fig.savefig(FIG / "cascade_bars.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def _fig_coverage_vs_target(art_cond):
    """Does the model hit its own ART coverage input, in every survey year?

    This is the second 95 on the survey's own (unconditional) basis, and it is
    the only cascade quantity with data in all three survey years. It is also
    the shortfall test: `data/art_coverage.csv` is an *input*, so a gap here is
    not a fit failure but the coverage target failing to fill -- which stisim
    does silently.

    Sex convention: `data/art_coverage.csv` is Gender 0 = FEMALE (stisim's, as
    documented in art_coverage_construction.py). Note this is the OPPOSITE of
    `calibration_data/art_coverage_by_age_sex.csv`, which follows the PHIA
    source's 0 = male. Getting this backwards silently swaps the two sexes and
    still produces a plausible-looking figure.
    """
    tgt = pd.read_csv(REPO / "data" / "art_coverage.csv")
    tgt["sex"] = tgt.Gender.map({0: "f", 1: "m"})
    tgt = tgt.rename(columns={"Year": "year", "AgeBin": "bin", "p_art": "target"})

    model = art_cond[art_cond.measure == "on_art_of_all_plhiv"]
    cov = model.merge(tgt[["year", "sex", "bin", "target"]],
                      on=["year", "sex", "bin"], how="inner")
    cov["shortfall_pp"] = (cov.model - cov.target) * 100

    years = [2011, 2016, 2021]
    bins = list(ART_BINS)
    x = np.arange(len(bins))
    fig, axes = plt.subplots(1, len(years), figsize=(14.5, 4.8), sharey=True)
    for ax, yr in zip(axes, years):
        sub = cov[cov.year == yr]
        for sex in ("f", "m"):
            s = sub[sub.sex == sex].set_index("bin").reindex(bins)
            ax.fill_between(x, s.model - s.model_sd, s.model + s.model_sd,
                            color=COLS[sex], alpha=0.15, lw=0, zorder=1)
            ax.plot(x, s.model, "-o", color=COLS[sex], lw=2, ms=8, zorder=3,
                    mec="white", mew=1.2)
            ax.plot(x, s.target, "s", ms=9, mfc="none", mew=2.0,
                    color=COLS[sex], zorder=4)
        ax.set_xticks(x)
        ax.set_xticklabels(bins, fontsize=8.5)
        ax.set_title(str(yr), fontsize=11, color=INK)
        ax.set_xlabel("Age stratum", fontsize=9.5, color=MUTED)
        ax.grid(axis="y", color="#e8e8e8", lw=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("On ART among people living with HIV",
                       fontsize=9.5, color=MUTED)
    handles = [plt.Line2D([], [], color=COLS[s], lw=2, marker="o", ms=8,
                          mec="white", label=f"Model - {SEXLAB[s].lower()}")
               for s in ("f", "m")]
    handles += [plt.Line2D([], [], color=COLS[s], lw=0, marker="s", ms=9,
                           mfc="none", mew=2.0,
                           label=f"Input target - {SEXLAB[s].lower()}")
                for s in ("f", "m")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("ART coverage: does the model reach the target it was given?",
                 fontsize=12.5, color=INK, y=1.0)
    fig.tight_layout()
    fig.savefig(FIG / "coverage_vs_target.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return cov


def _fig_standard(df):
    """The standard prevalence-fit figure every experiment owes (CLAUDE.md)."""
    try:
        import standard_figures as SF
        d = df.copy()
        d["arm"] = "baseline"
        sc = SF.plot_prevalence_fit(d, "029 baseline",
                                    FIG / "prevalence_fit_vs_phia.png")
        print("\nprevalence scorecard:", sc)
    except Exception as e:  # noqa: BLE001
        print(f"\n[warn] standard prevalence figure failed: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
