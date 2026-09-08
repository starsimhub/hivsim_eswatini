"""Summary figures for the Eswatini calibration presentation.

Four charts that no single experiment produces, because each one is a
cross-experiment view:

  figures/fit_progression.png     fit vs PHIA across the milestone configurations
  figures/tradeoff.png            the deaths-vs-prevalence trade-off, measured 3 ways
  figures/coverage.png            the three coverage checks
  figures/constrained_dims.png    wave 1's constrained directions, compact
  figures/journey.png             the five-month arc against the three coverage checks
  figures/sigma_scan.png          024's discrepancy scan, built large for the story deck
  figures/mortality_fix.png       016's double-counted mortality: the problem and the fix
  figures/vls.png                 021's viral-suppression input, assumed vs measured vs survey

Every number is read from the experiment record: `mae`/`bias`/`n_within_ci` come
from `hivsim_eswatini/standard_figures.scorecard` applied to the surviving
parquet outputs, so they are the same quantities each experiment's own SUMMARY
reports. Values that cannot be recomputed (peak-death shares, the package's PCA)
are transcribed from the SUMMARY tables and marked as such below.

Run: `python presentations/make_figures.py [--theme neutral|bmgf]` from the
HIVsim root. `bmgf` writes to `figures/bmgf/` and uses the Gates Foundation
brand palette pulled from the `BMGF Layouts` theme; `neutral` writes to
`figures/` and uses the dataviz reference palette.
"""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def _up(start, test):
    """Nearest ancestor of `start` (inclusive) satisfying `test`."""
    for d in [start, *start.parents]:
        if test(d):
            return d
    raise FileNotFoundError(f"no ancestor of {start} matched")


PRES = Path(__file__).resolve().parent
# The model repo is whichever ancestor holds experiments/ -- resolved rather
# than counted, so moving presentations/ does not break the scripts again.
REPO = _up(PRES, lambda d: (d / "experiments").is_dir())
EXP = REPO / "experiments"
ROOT = REPO

sys.path.insert(0, str(REPO))
from standard_figures import load_targets, scorecard  # noqa: E402

# --- Themes -------------------------------------------------------------------
# Both palettes are validated all-pairs on a light surface with the dataviz
# validator. Three categorical slots each; do not add a fourth -- past three,
# the all-pairs CVD and normal-vision floors stop clearing.
#
# neutral: dataviz reference slots 1-3.  CVD dE 9.2, normal-vision dE 24.0.
# bmgf:    accent4 / accent1 / accent5 of the `BMGF Layouts` theme, read out of
#          the foundation template rather than recalled. Blue-Orange-Turquoise
#          gives CVD dE 16.2, normal-vision dE 23.1. Turquoise sits below 3:1
#          contrast, so every chart using it carries direct labels -- the
#          dataviz "relief rule".
THEMES = {
    "neutral": dict(
        blue="#2a78d6", orange="#eb6834", aqua="#1baf7a", red="#e34948",
        ink="#22221f", ink2="#6b6a63", grid="#b8b7b0", edge="#c9c8c1",
        family=["Segoe UI", "DejaVu Sans"], subdir="",
    ),
    "bmgf": dict(
        blue="#248AF9", orange="#F85C02", aqua="#3AC9B1", red="#D93027",
        ink="#303A44", ink2="#6E7681", grid="#C7C4BA", edge="#C7C4BA",
        family=["Calibri", "Segoe UI", "DejaVu Sans"], subdir="bmgf",
    ),
}

_p = argparse.ArgumentParser()
_p.add_argument("--theme", choices=sorted(THEMES), default="neutral")
THEME = _p.parse_args().theme
T = THEMES[THEME]

OUT = PRES / "figures" / T["subdir"]
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, AQUA, RED = T["blue"], T["orange"], T["aqua"], T["red"]
INK, INK2 = T["ink"], T["ink2"]
GRID = dict(alpha=0.25, linewidth=0.7, color=T["grid"])

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": T["edge"], "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "font.family": "sans-serif", "font.sans-serif": T["family"],
    "font.size": 11, "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False,
})


# ----------------------------------------------------------------------------
# 1. Fit progression
# ----------------------------------------------------------------------------
# Each entry names the experiment's adopted or headline configuration. `kind`
# is the *object* being scored, and it is not cosmetic: an ensemble median over
# draws and a mean over seeds at one parameter point are different things, so
# the chart colours by kind rather than pretending one axis of progress.
SPECS = [
    dict(exp="009_coverage_check", file="ensemble.parquet", kind="ensemble",
         label="009\nprior", note="50-draw prior, 6 par"),
    dict(exp="014_prior_expansion", file="ensemble.parquet", kind="ensemble",
         label="014\nprior", note="50-draw prior, 9 par"),
    dict(exp="018_adopt_and_size", file="adopt.parquet", kind="arm",
         key="tag", select="high_transmission",
         label="018\nv1.1", note="fixed point"),
    dict(exp="021_vls_input", file="results.parquet", kind="arm",
         key="arm", select="B_vls_phia",
         label="021\nv1.2", note="fixed point"),
    dict(exp="022_survival_and_vls", file="results.parquet", kind="arm",
         key="arm", select="A_base",
         label="022\nv1.3", note="fixed point"),
    dict(exp="024_hm_wave1", file="ensemble.parquet", kind="ensemble",
         label="024\nprior", note="1000-draw prior, 7 par"),
    dict(exp="024_hm_wave1", file="ensemble.parquet", kind="best",
         key="point", select=868,
         label="024\nbest", note="calibrated point"),
]

KIND_COLOR = {"ensemble": BLUE, "arm": ORANGE, "best": AQUA}
KIND_LABEL = {"ensemble": "prior ensemble (median over draws)",
              "arm": "fixed parameter point (mean over seeds)",
              "best": "best joint draw of 1000"}


def collect_scorecards():
    tg = load_targets()
    rows = []
    for s in SPECS:
        path = EXP / s["exp"] / "outputs" / s["file"]
        if not path.exists():
            print(f"  skip {s['exp']}: {s['file']} not on disk")
            continue
        df = pd.read_parquet(path)
        if s.get("select") is not None:
            df = df[df[s["key"]] == s["select"]]
        # 'best' is a single point scored like an arm: mean over its rows.
        kind = "ensemble" if s["kind"] == "ensemble" else "arm"
        sc = scorecard(df, kind=kind, tg=tg)
        rows.append(dict(label=s["label"], note=s["note"], kind=s["kind"], **sc))
        print(f"  {s['exp']:<28} {s['kind']:<9} "
              f"MAE {sc.get('mae', float('nan')):.4f} "
              f"bias {sc.get('bias', float('nan')):+.4f} "
              f"within-CI {sc.get('n_within_ci', 0)}/54")
    return pd.DataFrame(rows)


def fig_fit_progression(sc):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.1))
    x = np.arange(len(sc))
    colors = [KIND_COLOR[k] for k in sc.kind]

    def bars(ax, col, title, ylabel, fmt="{:.3f}", invert=False):
        v = sc[col].to_numpy(dtype=float)
        ax.bar(x, v, color=colors, width=0.68, zorder=3)
        for xi, vi in zip(x, v):
            va = "top" if vi < 0 else "bottom"
            off = -3 if vi < 0 else 3
            ax.annotate(fmt.format(vi), (xi, vi), textcoords="offset points",
                        xytext=(0, off), ha="center", va=va, fontsize=9.5,
                        color=INK)
        ax.set_xticks(x)
        ax.set_xticklabels(sc.label, fontsize=10, linespacing=1.3)
        ax.set(title=title, ylabel=ylabel)
        ax.grid(axis="y", zorder=0, **GRID)
        if invert:
            ax.axhline(0, color="#8c8b84", linewidth=0.9, zorder=4)

    bars(axes[0], "mae",
         "Mean absolute error vs PHIA\n54 age x sex x year strata",
         "MAE, prevalence")
    bars(axes[1], "bias",
         "Mean bias vs PHIA\nnegative = model below the surveys",
         "bias, prevalence", fmt="{:+.3f}", invert=True)
    # Clipped so the 018 -> 024 story (-0.040 -> -0.007) stays legible. 014's
    # -0.233 is annotated as off-scale rather than dropped.
    axes[1].set_ylim(-0.10, 0.007)
    axes[1].annotate("014 bias is −0.233 —\nclipped, off scale",
                     (1.45, -0.079), ha="left", va="center", fontsize=9.5,
                     color=INK2, linespacing=1.4, zorder=5)
    bars(axes[2], "n_within_ci",
         "Strata inside the PHIA 95% CI\nout of 54",
         "strata", fmt="{:.0f}")
    axes[2].set_ylim(0, 54)

    handles = [plt.Rectangle((0, 0), 1, 1, color=KIND_COLOR[k])
               for k in ("ensemble", "arm", "best")]
    fig.legend(handles, [KIND_LABEL[k] for k in ("ensemble", "arm", "best")],
               loc="lower center", ncol=3, frameon=False, fontsize=10,
               bbox_to_anchor=(0.5, -0.015))
    fig.suptitle("Fit to PHIA prevalence across the milestone configurations",
                 fontsize=14, fontweight="bold", y=0.99)
    fig.text(0.5, 0.055,
             "Ensembles and fixed points are different objects and are coloured "
             "as such; 016 and 017 are absent because their outputs kept only "
             "aggregate infected counts, so per-band prevalence is unrecoverable.",
             ha="center", fontsize=8.5, color=INK2)
    fig.tight_layout(rect=(0, 0.095, 1, 0.96))
    fig.savefig(OUT / "fit_progression.png", dpi=200)
    plt.close(fig)
    print("  -> figures/fit_progression.png")


# ----------------------------------------------------------------------------
# 2. The deaths / prevalence trade-off
# ----------------------------------------------------------------------------
# Transcribed from the SUMMARY tables of 019 (obs 2, obs 6) and 022 (obs 4, and
# the scorecard). Peak AIDS deaths as a share of the UNAIDS peak is not
# recomputable from the parquet without re-deriving the peak alignment, so it
# is quoted. Both experiments ran at the same high-transmission parameter point.
TRADEOFF = pd.DataFrame([
    # exp, arm label, % of UNAIDS peak deaths, PHIA MAE, label offset in points
    ("019", "A  flat 13.1 y (upstream)", 64.3, 0.0590, (14, 6), "left"),
    ("019", "B  flat 11.5 y", 73.8, 0.0657, (14, -2), "left"),
    ("019", "C  mild age gradient", 78.2, 0.0725, (-14, -16), "right"),
    ("019", "D  full ALPHA gradient", 78.1, 0.0782, (-14, 4), "right"),
    ("022", "A  baseline", 64.3, 0.0584, (14, -10), "left"),
    ("022", "B  EMOD pivoting gradient", 73.3, 0.0817, (14, 0), "left"),
], columns=["exp", "arm", "pct_unaids", "mae", "offset", "ha"])


def fig_tradeoff():
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    for exp, color, marker, lab in (
            ("019", BLUE, "o", "019 — shorten survival for everyone, then by age"),
            ("022", ORANGE, "s", "022 — longer when young, shorter when older")):
        g = TRADEOFF[TRADEOFF.exp == exp].sort_values("pct_unaids")
        ax.plot(g.pct_unaids, g.mae, marker=marker, markersize=10,
                color=color, linewidth=2, markeredgecolor="white",
                markeredgewidth=1.6, zorder=3, label=lab)
    for r in TRADEOFF.itertuples():
        ax.annotate(r.arm, (r.pct_unaids, r.mae), textcoords="offset points",
                    xytext=r.offset, ha=r.ha, va="center", fontsize=13.8,
                    color=INK2)

    ax.axhline(0.0584, color="#8c8b84", linestyle=":", linewidth=1.2, zorder=1)
    ax.annotate("where we started", (81.6, 0.0584),
                textcoords="offset points", xytext=(0, 5), ha="right",
                fontsize=13.0, color=INK2)
    ax.set(xlabel="AIDS deaths at the peak, % of the reported figure  →  better",
           ylabel="how far prevalence sits from the surveys  →  worse",
           title="Every way of recovering deaths costs us on prevalence",
           xlim=(61, 82), ylim=(0.0555, 0.0855))
    ax.grid(zorder=0, **GRID)
    ax.legend(frameon=False, fontsize=14.5, loc="upper left")
    fig.text(0.5, 0.015,
             "Three different ways of shortening survival all land on the same line — that makes it a real\n"
             "feature of how the model works, not a quirk of one set of numbers.",
             ha="center", fontsize=12.3, color=INK2, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    fig.savefig(OUT / "tradeoff.png", dpi=200)
    plt.close(fig)
    print("  -> figures/tradeoff.png")


# ----------------------------------------------------------------------------
# 3. The three coverage checks
# ----------------------------------------------------------------------------
# From the SUMMARY headlines of 009, 014 and 024. The target sets differ (89
# rows for 009/014, 48 registered features for 024) and that is stated on the
# chart rather than hidden by plotting only percentages.
COVERAGE = pd.DataFrame([
    ("009\n2026-05-07", 30, 89,
     "6 parameters\nbroken network\nbroken VMMC"),
    ("014\n2026-08-11", 4, 89,
     "9 parameters, incl. mortality\nprior sampled below the\nestablishment threshold"),
    ("024\n2026-09-03", 45, 48,
     "7 pruned parameters\nmodel-v1.3\nnetwork + VMMC fixed"),
], columns=["label", "inside", "total", "note"])


def fig_coverage():
    fig, ax = plt.subplots(figsize=(10.2, 5.8))
    pct = 100 * COVERAGE.inside / COVERAGE.total
    x = np.arange(len(COVERAGE))
    ax.bar(x, pct, color=[RED, RED, AQUA], width=0.58, zorder=3)
    for xi, p, r in zip(x, pct, COVERAGE.itertuples()):
        ax.annotate(f"{p:.0f}%", (xi, p), textcoords="offset points",
                    xytext=(0, 22), ha="center", fontsize=17,
                    fontweight="bold", color=INK)
        ax.annotate(f"{r.inside}/{r.total} targets", (xi, p),
                    textcoords="offset points", xytext=(0, 6), ha="center",
                    fontsize=10.5, color=INK2)
        # Notes sit below the axis so they never depend on bar height.
        ax.annotate(r.note, (xi, 0), xycoords=("data", "axes fraction"),
                    textcoords="offset points", xytext=(0, -52), ha="center",
                    va="top", fontsize=9.5, color=INK2, linespacing=1.45,
                    annotation_clip=False)
    ax.set_xticks(x)
    ax.set_xticklabels(COVERAGE.label, fontsize=12, linespacing=1.4)
    ax.set(ylabel="inside the 5–95% prior envelope (%)", ylim=(0, 112))
    ax.set_title("The prior predictive check, three attempts", pad=34)
    ax.grid(axis="y", zorder=0, **GRID)
    fig.text(0.5, 0.038,
             "009 and 014 score 89 target rows (PHIA prevalence + UNAIDS deaths); 024 scores 48 registered\n"
             "features on a settled target set, so the percentages are not strictly like-for-like. 024's best\n"
             "single draw hits 48/48 within 3σ — the model can reach the data jointly, not just target by target.",
             ha="center", fontsize=8.5, color=INK2, linespacing=1.6)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.87, bottom=0.36)
    fig.savefig(OUT / "coverage.png", dpi=200)
    plt.close(fig)
    print("  -> figures/coverage.png")


# ----------------------------------------------------------------------------
# 4. Constrained directions, compact
# ----------------------------------------------------------------------------
# Transcribed from the history_matching package's own wave-1 diagnostic
# (024/outputs/hm/wave1/wave1/constrained_dims.png), reproduced compactly for a
# slide. Values match 024's SUMMARY observation 2.
PC_REDUCTION = [35.7, 9.3, 5.5, 0.0, 0.0, 0.0, 0.0]
PC1_LOADINGS = pd.DataFrame([
    ("beta_m2f", 0.67), ("rel_beta_f2m", 0.54), ("s_f_young", 0.35),
    ("age_gap_shift", 0.18), ("age_gap_sd_mult", 0.13),
    ("prop_m0", -0.25), ("prop_f0", -0.02),
], columns=["par", "loading"])


def fig_constrained_dims():
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.0),
                             gridspec_kw=dict(width_ratios=[1, 1.25]))

    ax = axes[0]
    x = np.arange(7)
    ax.bar(x, PC_REDUCTION, color=BLUE, width=0.66, zorder=3)
    for xi, v in zip(x, PC_REDUCTION):
        ax.annotate(f"{v:.1f}%", (xi, v), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=10, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([f"PC{i+1}" for i in range(7)], fontsize=10)
    ax.set(ylabel="variance reduction, prior → NROY (%)", ylim=(0, 44),
           title="One observation constrains one direction")
    ax.grid(axis="y", zorder=0, **GRID)

    ax = axes[1]
    g = PC1_LOADINGS
    y = np.arange(len(g))[::-1]
    colors = [RED if v > 0 else BLUE for v in g.loading]
    ax.barh(y, g.loading, color=colors, height=0.62, zorder=3)
    for yi, v, p in zip(y, g.loading, g.par):
        ha = "left" if v > 0 else "right"
        ax.annotate(f"{v:+.2f}", (v, yi), textcoords="offset points",
                    xytext=(5 if v > 0 else -5, 0), ha=ha, va="center",
                    fontsize=10, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(g.par, fontsize=10.5)
    ax.axvline(0, color="#8c8b84", linewidth=0.9, zorder=4)
    ax.set(xlabel="loading on PC1", xlim=(-0.42, 0.92),
           title="PC1 is composite transmission intensity")
    ax.grid(axis="x", zorder=0, **GRID)

    fig.suptitle("Wave 1's cut is diagonal in the box — the marginals barely move",
                 fontsize=14, fontweight="bold", y=0.995)
    fig.text(0.5, 0.02,
             "beta_m2f, rel_beta_f2m and s_f_young load on PC1 with the same sign: they enter the emulated feature as a product.\n"
             "PC4–PC7 are untouched, so reporting per-parameter marginal intervals would have made this wave look like it did nothing.",
             ha="center", fontsize=8.5, color=INK2, linespacing=1.6)
    fig.tight_layout(rect=(0, 0.10, 1, 0.94))
    fig.savefig(OUT / "constrained_dims.png", dpi=200)
    plt.close(fig)
    print("  -> figures/constrained_dims.png")


# ----------------------------------------------------------------------------
# 5. The journey: five months of experiments against the checkpoint results
# ----------------------------------------------------------------------------
# Dates are each experiment's own SUMMARY date (001-007 from experiments/log.md).
# Categories are a reading of what each one delivered, not a field in the record:
#   change  - a model change that was adopted
#   closed  - a hypothesis tested and closed, with no gain in fit
#   setup   - targets, tooling, sizing, parameter selection; no model change
# 015 sits before 014 deliberately: it ran a month earlier than its number.
JOURNEY = [
    ("001", "2026-04-07", "change"), ("002", "2026-04-14", "change"),
    ("003", "2026-04-14", "change"), ("004", "2026-05-04", "change"),
    ("005", "2026-05-06", "change"), ("006", "2026-05-06", "change"),
    ("007", "2026-05-06", "change"), ("008", "2026-05-07", "setup"),
    ("009", "2026-05-07", "closed"), ("011", "2026-06-16", "change"),
    ("012", "2026-06-17", "setup"),  ("013", "2026-07-09", "closed"),
    ("015", "2026-07-10", "change"), ("014", "2026-08-11", "closed"),
    ("016", "2026-08-19", "change"), ("017", "2026-08-26", "closed"),
    ("018", "2026-08-27", "setup"),  ("019", "2026-09-01", "closed"),
    ("020", "2026-09-01", "setup"),  ("021", "2026-09-01", "change"),
    ("022", "2026-09-02", "closed"), ("023", "2026-09-02", "setup"),
    ("024", "2026-09-03", "wave"),
]
# The three prior predictive checks -- the only checkpoint measured the same way
# at three points in time, so the only honest "did it get better" series here.
CHECKS = [("009", "2026-05-07", 30, 89), ("014", "2026-08-11", 4, 89),
          ("024", "2026-09-03", 45, 48)]

CAT = {
    "change": ("Changed the model\nand kept the change", AQUA, "o"),
    "closed": ("Tested an explanation\nand ruled it out", ORANGE, "s"),
    # Hexagon rather than a triangle: a triangle is too narrow at its vertical
    # centre to hold a three-digit label, which silently clipped 008/012/018/
    # 020/023 down to a single character.
    "setup":  ("Set the ground rules\ndata, model size, settings", BLUE, "h"),
}


def _spread(d, g, run, days, gap):
    """Place the members of one collision run at even spacing about their mean."""
    if len(run) < 2:
        return
    centre = days[run].mean()
    offs = (np.arange(len(run)) - (len(run) - 1) / 2) * gap
    for k, off in zip(run, offs):
        d.loc[g.index[k], "x"] = pd.Timestamp.fromordinal(
            int(round(centre + off)))


def fig_journey():
    d = pd.DataFrame(JOURNEY, columns=["exp", "date", "cat"])
    d["date"] = pd.to_datetime(d.date)
    c = pd.DataFrame(CHECKS, columns=["exp", "date", "inside", "total"])
    c["date"] = pd.to_datetime(c.date)
    c["pct"] = 100 * c.inside / c.total

    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=(14.2, 5.5), sharex=True,
        gridspec_kw=dict(height_ratios=[1.15, 1], hspace=0.16))

    # --- top: every experiment on one row per category -----------------------
    # A marker is about three days wide on a five-month axis, so experiments a
    # day or two apart still collide. Group each category's points into runs
    # that fall within a marker width and spread each run symmetrically about
    # its mean date. Positions are therefore nudged by a few days; the footnote
    # says so.
    GAP = 4.2
    d = d.sort_values("date").copy()
    d["x"] = d.date
    for cat, g in d.groupby("cat"):
        g = g.sort_values("date")
        days = g.date.map(pd.Timestamp.toordinal).to_numpy(dtype=float)
        run = [0]
        for k in range(1, len(days)):
            if days[k] - days[run[-1]] < GAP:
                run.append(k)
                continue
            _spread(d, g, run, days, GAP)
            run = [k]
        _spread(d, g, run, days, GAP)

    order = ["change", "closed", "setup"]
    for i, cat in enumerate(order):
        g = d[d.cat == cat]
        label, colour, marker = CAT[cat]
        ax.scatter(g.x, [i] * len(g), s=330, marker=marker, color=colour,
                   zorder=3, edgecolor="white", linewidth=1.4, label=label)
        for r in g.itertuples():
            ax.annotate(r.exp, (r.x, i), ha="center", va="center",
                        fontsize=9.8, color="white", fontweight="bold",
                        zorder=4)

    # Wave 1 gets a rule through both panels rather than a marker in a row --
    # it is what the whole timeline is heading towards, and the rule ties it to
    # the 94% below.
    w = d[d.cat == "wave"].iloc[0]
    for a in (ax, bx):
        a.axvline(w.date, color=INK, linestyle=(0, (3, 3)), linewidth=1.3,
                  zorder=1, alpha=0.6)
    ax.annotate("024 — first calibration run", (w.date, 2.42),
                textcoords="offset points", xytext=(-9, 0), ha="right",
                va="center", fontsize=13.7, color=INK, fontweight="bold")

    ax.set_yticks(range(len(order)))
    counts = d.cat.value_counts()
    ax.set_yticklabels(
        [f"{CAT[c_][0]}   ({counts.get(c_, 0)})" for c_ in order],
        fontsize=12.3, linespacing=1.4)
    ax.set_ylim(-0.7, 2.7)
    ax.grid(axis="x", zorder=0, **GRID)
    ax.set_title("Twenty-three experiments, April to September 2026",
                 loc="left")
    for spine in ("left",):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)

    # --- bottom: the checkpoint, measured three times ------------------------
    bx.plot(c.date, c.pct, "-", color="#B9B6AC", linewidth=1.6, zorder=2)
    bx.scatter(c.date, c.pct, s=240, zorder=3, edgecolor="white",
               linewidth=1.5, color=[RED, RED, AQUA])
    for r in c.itertuples():
        bx.annotate(f"{r.pct:.0f}%", (r.date, r.pct), textcoords="offset points",
                    xytext=(0, 13), ha="center", fontsize=16.9, fontweight="bold",
                    color=INK)
        bx.annotate(f"{r.exp} · {r.inside}/{r.total}", (r.date, r.pct),
                    textcoords="offset points", xytext=(0, -20), ha="center",
                    fontsize=12.3, color=INK2)
    bx.set(ylabel="data points the model\ncould reach (%)",
           ylim=(-14, 124))
    bx.grid(axis="both", zorder=0, **GRID)
    bx.set_title("The checkpoint: can the model produce the data at all?",
                 loc="left")

    import matplotlib.dates as mdates
    bx.xaxis.set_major_locator(mdates.MonthLocator())
    bx.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    bx.tick_params(axis="x", labelsize=10.5)

    fig.text(0.5, 0.015,
             "Four months of work between the first check and the second, and the second scored worse. "
             "The gain came at 024, the first run that let the transmission settings move.\n"
             "The first two checks score 89 data points, the third scores 48 summary numbers, so the "
             "percentages are not exactly like-for-like. Markers are nudged a few days where dates collide.",
             ha="center", fontsize=11.1, color=INK2, linespacing=1.5)
    fig.subplots_adjust(left=0.175, right=0.99, top=0.90, bottom=0.19)
    fig.savefig(OUT / "journey.png", dpi=200)
    plt.close(fig)
    print("  -> figures/journey.png")


# ----------------------------------------------------------------------------
# 6. The discrepancy-allowance scan -- the evidence for 024's reversal
# ----------------------------------------------------------------------------
# Straight from 024's outputs/best_point.json. Replaces the borrowed
# target_residuals.png on the story deck, whose per-target axis labels are
# unreadable at slide size and whose top panel is not the argument being made.
SIGMA_SCAN = {0.000: 3.98, 0.005: 3.01, 0.010: 1.99, 0.015: 1.43, 0.020: 1.12}


def fig_sigma_scan():
    x = np.array(list(SIGMA_SCAN))
    y = np.array(list(SIGMA_SCAN.values()))
    fig, ax = plt.subplots(figsize=(8.4, 5.4))

    ax.axhspan(0, 3, color=AQUA, alpha=0.10, zorder=0)
    ax.axhline(3, color=INK2, linestyle=":", linewidth=1.6, zorder=2)
    ax.annotate("close enough", (0.0203, 3), fontsize=13, color=INK2,
                va="center", ha="right", xytext=(0, -14),
                textcoords="offset points")

    ax.plot(x, y, "-o", color=BLUE, linewidth=2.6, markersize=11,
            markeredgecolor="white", markeredgewidth=1.8, zorder=3)
    for xi, yi in zip(x, y):
        ax.annotate(f"{yi:.2f}", (xi, yi), textcoords="offset points",
                    xytext=(0, 15), ha="center", fontsize=13.5,
                    fontweight="bold", color=INK)

    ax.scatter([0], [3.98], s=320, facecolor="none", edgecolor=ORANGE,
               linewidth=2.6, zorder=4)
    ax.annotate("no allowance at all — the best\nversion is still 4 errors short",
                (0, 3.98), textcoords="offset points", xytext=(26, -6),
                ha="left", va="center", fontsize=13, color=ORANGE,
                linespacing=1.35)

    ax.set(xlabel="allowance made for the model being wrong",
           ylabel="how far off the best version is\n(standard errors)",
           xlim=(-0.0016, 0.0225), ylim=(0, 4.7))
    ax.set_xticks(list(SIGMA_SCAN))
    ax.tick_params(labelsize=13)
    ax.xaxis.label.set_size(13.5)
    ax.yaxis.label.set_size(13.5)
    ax.set_title("The model was not built wrong", fontsize=16, loc="left")
    ax.grid(zorder=0, **GRID)

    fig.text(0.5, 0.02,
             "We set aside an allowance to absorb a shortfall seven experiments had "
             "blamed on the model being built wrong.\nRemove the allowance and the best "
             "version is still four standard errors short — it was the settings, not the model.",
             ha="center", fontsize=11, color=INK2, linespacing=1.5)
    fig.subplots_adjust(left=0.115, right=0.985, top=0.90, bottom=0.245)
    fig.savefig(OUT / "sigma_scan.png", dpi=200)
    plt.close(fig)
    print("  -> figures/sigma_scan.png")


# ----------------------------------------------------------------------------
# 7. Double-counted mortality: the problem and the fix
# ----------------------------------------------------------------------------
def fig_mortality_fix():
    exp = EXP / "016_double_counted_mortality" / "outputs"
    implied = pd.read_csv(exp / "implied_aids_deaths.csv")
    calib = pd.read_csv(REPO / "data" / "eswatini_hiv_calib.csv")
    res = pd.read_parquet(exp / "results.parquet")
    # Default parameters: the arm 016 headlines, and the one the slide quotes.
    res = res[res.pset == "default"]
    pop = (res.groupby(["arm", "timevec"])["popagesex.n_alive_total"]
              .mean().reset_index())

    # Stacked, not side by side: the slide puts explanatory text beside this,
    # so it needs to be tall rather than wide.
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.2, 6.0))

    # --- left: the problem ---------------------------------------------------
    unaids = calib.dropna(subset=["hiv.new_deaths"])
    ax.plot(implied.year, implied.implied_aids_deaths, "-", color=RED,
            linewidth=3.0, zorder=3,
            label="AIDS deaths hidden in the\nbackground death rates")
    ax.plot(unaids.time, unaids["hiv.new_deaths"], "--", color=INK,
            linewidth=2.4, zorder=4, label="AIDS deaths reported by UNAIDS")
    ax.set(xlabel="year", ylabel="AIDS deaths per year", xlim=(1985, 2025))
    ax.set_title("The problem: the death data already contained the epidemic",
                 fontsize=13.5, loc="left")
    ax.legend(frameon=False, fontsize=11.5, loc="upper left")
    ax.grid(zorder=0, **GRID)
    ax.tick_params(labelsize=12.5)
    ax.xaxis.label.set_size(13); ax.yaxis.label.set_size(13)

    # --- right: the fix ------------------------------------------------------
    style = {"all_cause": (RED, "AIDS deaths left in (what we had)"),
             "hiv_deleted": (AQUA, "AIDS deaths taken out (the fix)")}
    for arm, (colour, label) in style.items():
        g = pop[pop.arm == arm]
        bx.plot(g.timevec, g["popagesex.n_alive_total"] / 1e6, "-",
                color=colour, linewidth=3.0, zorder=3, label=label)
    tgt = calib.dropna(subset=["n_alive"])
    bx.plot(tgt.time, tgt.n_alive / 1e6, "o", color=INK, markersize=7,
            markeredgecolor="white", markeredgewidth=1.2, zorder=4,
            label="recorded population")
    bx.set(xlabel="year", ylabel="population (millions)", xlim=(1985, 2025))
    bx.set_title("The fix: taking them out puts the population back on track",
                 fontsize=13.5, loc="left")
    bx.legend(frameon=False, fontsize=11.5, loc="upper left")
    bx.grid(zorder=0, **GRID)
    bx.tick_params(labelsize=12.5)
    bx.xaxis.label.set_size(13); bx.yaxis.label.set_size(13)

    fig.text(0.5, 0.02,
             "Top: worked out from the death data alone, with no HIV information used.  Bottom: model population at default settings.",
             ha="center", fontsize=9.5, color=INK2, linespacing=1.5)
    fig.subplots_adjust(left=0.135, right=0.98, top=0.93, bottom=0.145,
                        hspace=0.62)
    fig.savefig(OUT / "mortality_fix.png", dpi=200)
    plt.close(fig)
    print("  -> figures/mortality_fix.png")


# ----------------------------------------------------------------------------
# 8. Viral suppression: what assuming everyone is suppressed costs
# ----------------------------------------------------------------------------
def fig_vls():
    src = (EXP / "021_vls_input" / "outputs" / "cascade_vs_phia.csv")
    d = pd.read_csv(src)
    d = d[d.quantity == "vls_among_plhiv"]
    keys = [(2016, "m"), (2016, "f"), (2021, "m"), (2021, "f")]
    labels = ["2016\nmen", "2016\nwomen", "2021\nmen", "2021\nwomen"]

    def series(arm):
        return [float(d[(d.arm == arm) & (d.year == y) & (d.sex == sx)]
                      .model.iloc[0]) * 100 for y, sx in keys]

    assumed = series("A_vls_1.0")
    supplied = series("B_vls_phia")
    survey = [float(d[(d.arm == "B_vls_phia") & (d.year == y) & (d.sex == sx)]
                    .phia.iloc[0]) * 100 for y, sx in keys]

    x = np.arange(len(keys))
    w = 0.27
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    bars = [
        (x - w, assumed, RED, "assumed everyone suppressed"),
        (x, supplied, AQUA, "using measured suppression"),
        (x + w, survey, INK, "what the surveys found"),
    ]
    for xs, vals, colour, label in bars:
        ax.bar(xs, vals, width=w * 0.92, color=colour, zorder=3, label=label)
        for xi, v in zip(xs, vals):
            ax.annotate(f"{v:.0f}", (xi, v), textcoords="offset points",
                        xytext=(0, 4), ha="center", fontsize=11.5, color=INK)

    # The gap that matters: assumed-minus-survey, in percentage points.
    for xi, a, sv in zip(x, assumed, survey):
        ax.annotate(f"+{a - sv:.1f} pts", (xi - w, a), textcoords="offset points",
                    xytext=(0, 20), ha="center", fontsize=12, color=RED,
                    fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=13, linespacing=1.3)
    ax.set(ylabel="% of people with HIV who are virally suppressed",
           ylim=(0, 112))
    ax.yaxis.label.set_size(13)
    ax.tick_params(labelsize=12.5)
    ax.set_title("Assuming everyone on treatment is suppressed overstated the "
                 "cascade", fontsize=15, loc="left")
    ax.grid(axis="y", zorder=0, **GRID)
    ax.legend(frameon=False, fontsize=12.5, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.13), handletextpad=0.5,
              columnspacing=1.8)

    fig.text(0.5, 0.02,
             "Suppression among everyone living with HIV. The error is largest in 2016, when it mattered most: "
             "a treatment scenario\nmeasured against a baseline that already assumes near-perfect suppression has "
             "no room left to show a benefit.",
             ha="center", fontsize=9.5, color=INK2, linespacing=1.5)
    fig.subplots_adjust(left=0.10, right=0.985, top=0.91, bottom=0.30)
    fig.savefig(OUT / "vls.png", dpi=200)
    plt.close(fig)
    print("  -> figures/vls.png")


if __name__ == "__main__":
    print("Scoring milestone configurations:")
    sc = collect_scorecards()
    # Theme-independent: the numbers are the same whatever the palette is.
    sc.to_csv(PRES / "fit_progression_scorecard.csv",
              index=False)
    fig_fit_progression(sc)
    fig_tradeoff()
    fig_coverage()
    fig_constrained_dims()
    fig_journey()
    fig_sigma_scan()
    fig_mortality_fix()
    fig_vls()
