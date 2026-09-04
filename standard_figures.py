"""The standard prevalence-fit figure, produced identically by every experiment.

Every experiment's `run.py` calls `plot_prevalence_fit(...)`, and the retrospective
`plot_fit_progression.py` calls the same function. One implementation, so the
figure cannot drift between experiments.

Why this exists
---------------
`plot_dashboard.py` produced `dashboard_fit_*.png` for experiments 003-015 and
was the standard. It runs its own sims, so once experiments became A/B harnesses
with their own arms and configs (016 onward) it could not consume their output
and silently fell out of use. 016, 017, 018 and 020 have no prevalence-fit
figure at all as a result, and 019 and 021 grew a second, incompatible format.
This module replaces the ad-hoc versions; `plot_dashboard.py` remains for
full model-level checks.

What the figure shows
---------------------
Age-stratified prevalence against PHIA by sex and survey year, plus a 15-49
subpanel. The age-stratified view is the one that matters: aggregating to 15-49
averages young under-infection against older over-infection and so *hides* the
model's real defect. The 15-49 panel is included because experiments 016 and 017
can only produce that resolution -- their local `PopByAgeSex` recorded
`n_alive` per band but only aggregate infected -- so it keeps the whole series
comparable.

Sex coding
----------
`calibration_data/prevalence_by_age_sex.csv` uses **0 = Male, 1 = Female**, per
experiments/008 and confirmed against SDHS 2006-07 (Gender=1 gives 0.101 at
15-19 in 2007, the figure for women). This is the OPPOSITE of stisim's internal
convention, where 0 = female. Four files in this project use two conventions;
see exp 021's config.yaml for the full map.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
TARGETS = REPO / "calibration_data" / "prevalence_by_age_sex.csv"
PHIA_SEX = {0: "m", 1: "f"}
SEX_LABEL = {"f": "Women", "m": "Men"}


def load_targets():
    """The 54 PHIA target strata: year x sex x 5-year band, 15-65."""
    t = pd.read_csv(TARGETS)
    t["sex"] = t.Gender.map(PHIA_SEX)
    t["age_low"] = t["start age"].astype(int)
    return t[["Year", "sex", "age_low", "Count", "NationalPrevalence",
              "lb", "ub"]].rename(columns={"Year": "year",
                                           "NationalPrevalence": "phia"})


def targets_15_49(tg=None):
    """PHIA 15-49 by sex, weighted by the survey's own denominators."""
    tg = load_targets() if tg is None else tg
    a = tg[tg.age_low.between(15, 45)]
    return (a.assign(num=a.phia * a.Count)
             .groupby(["year", "sex"])
             .apply(lambda g: g.num.sum() / g.Count.sum(), include_groups=False)
             .rename("phia").reset_index())


def _stratum_series(df, sex, lo):
    """Per-seed prevalence for one stratum, or None if not recoverable.

    Handles both conventions: `popagesex.n_infected_*` (018 onward) and the
    pre-018 `hiv.prevalence_*` / `hiv_epi.prevalence_*` split.
    """
    icol, acol = (f"popagesex.n_infected_{sex}_{lo}_{lo + 5}",
                  f"popagesex.n_alive_{sex}_{lo}_{lo + 5}")
    if icol in df.columns and acol in df.columns:
        return df[icol] / df[acol].replace(0, np.nan)
    stem = "hiv" if lo < 35 else "hiv_epi"   # pre-018: hiv.* is 15-35 only
    col = f"{stem}.prevalence_{sex}_{lo}_{lo + 5}"
    return df[col] if col in df.columns else None


def _series_15_49(df, sex):
    """Per-seed prevalence 15-49 for one sex, or None."""
    col = f"popagesex.prevalence_{sex}_15_49"
    if col in df.columns:
        return df[col]
    icols = [c for c in df.columns
             if c.startswith(f"popagesex.n_infected_{sex}_")
             and 15 <= int(c.split("_")[-2]) < 50]
    if not icols:
        return None
    acols = [c.replace("n_infected", "n_alive") for c in icols]
    return df[icols].sum(axis=1) / df[acols].sum(axis=1).replace(0, np.nan)


def fit_by_stratum(df, kind="arm", tg=None):
    """Model vs PHIA per target stratum. Empty if age bands are unrecoverable.

    kind='arm'      -> mean +/- SD over seeds
    kind='ensemble' -> median and 5-95% over draws
    """
    tg = load_targets() if tg is None else tg
    rows = []
    for _, t in tg.iterrows():
        g = df[df.timevec == t.year]
        if not len(g):
            continue
        s = _stratum_series(g, t.sex, int(t.age_low))
        if s is None:
            continue
        s = s.dropna()
        if not len(s):
            continue
        if kind == "ensemble":
            centre, lo_b, hi_b = s.median(), s.quantile(0.05), s.quantile(0.95)
        else:
            # ddof=1 is NaN for a single seed; a zero-width band is the honest
            # rendering there, and it keeps single-seed smoke tests plottable.
            sd = s.std(ddof=1)
            sd = 0.0 if not np.isfinite(sd) else sd
            centre, lo_b, hi_b = s.mean(), s.mean() - sd, s.mean() + sd
        rows.append(dict(year=int(t.year), sex=t.sex, age_low=int(t.age_low),
                         model=centre, model_lo=lo_b, model_hi=hi_b,
                         phia=t.phia, phia_lb=t.lb, phia_ub=t.ub))
    return pd.DataFrame(rows)


FIT_15_49_COLS = ["year", "sex", "model", "model_sd", "phia", "weighting"]


def fit_15_49(df, tg=None, weighting="model_own"):
    """Model vs PHIA at 15-49 by sex and survey year.

    `weighting` controls how the model's 15-49 aggregate is formed, and the
    choice is not cosmetic -- for exp 018 the two routes differ by 0.0056, which
    is larger than the entire measured change across 016-021. Mixing them across
    experiments would manufacture a trend, so a cross-experiment comparison must
    pin one route.

    'model_own'    sum(n_infected) / sum(n_alive) over the model's own 15-49
                   bands, or `popagesex.prevalence_{sex}_15_49` where that is
                   all that was saved (016, 017). Available for every
                   fixed-parameter experiment, so this is the default and the
                   one the progression chart uses.
    'phia_wt'      model band prevalences weighted by PHIA's denominators --
                   identical weighting on both sides of the comparison, but it
                   needs per-band data, which 016 and 017 lack.

    The returned `weighting` column records what was actually used, since the
    requested route is not always available.
    """
    tg = load_targets() if tg is None else tg
    p = targets_15_49(tg)
    rows = []
    for _, t in p.iterrows():
        g = df[df.timevec == t.year]
        if not len(g):
            continue
        s, route = None, None
        if weighting == "model_own":
            s = _series_15_49(g, t.sex)
            route = "model_own"
        if s is None:      # phia_wt requested, or no model denominators saved
            bands = tg[(tg.year == t.year) & (tg.sex == t.sex)
                       & (tg.age_low.between(15, 45))]
            num, wt = None, 0.0
            for _, b in bands.iterrows():
                bs = _stratum_series(g, t.sex, int(b.age_low))
                if bs is None:
                    num = None
                    break
                num = bs * b.Count if num is None else num + bs * b.Count
                wt += b.Count
            if num is None or wt <= 0:
                continue
            s, route = num / wt, "phia_wt"
        s = s.dropna()
        if not len(s):
            continue
        rows.append(dict(year=int(t.year), sex=t.sex, model=s.mean(),
                         model_sd=s.std(ddof=1), phia=t.phia, weighting=route))
    return pd.DataFrame(rows, columns=FIT_15_49_COLS)


def scorecard(df, kind="arm", tg=None):
    """MAE and bias at both resolutions. The 15-49 numbers flatter the model."""
    out = {}
    st = fit_by_stratum(df, kind=kind, tg=tg)
    if len(st):
        err = st.model - st.phia
        inside = ((st.model >= st.phia_lb) & (st.model <= st.phia_ub)).sum()
        out.update(n_strata=len(st), mae=err.abs().mean(), bias=err.mean(),
                   n_within_ci=int(inside))
    agg = fit_15_49(df, tg=tg)
    if len(agg):
        err = agg.model - agg.phia
        out.update(n_15_49=len(agg), mae_15_49=err.abs().mean(),
                   bias_15_49=err.mean())
    return out


def plot_prevalence_fit(df, label, outpath, kind="arm", tg=None, stamp=None):
    """The standard figure. `df` is one configuration's rows across seeds/draws.

    Returns the scorecard dict, so callers can print or store it.
    """
    tg = load_targets() if tg is None else tg
    st = fit_by_stratum(df, kind=kind, tg=tg)
    agg = fit_15_49(df, tg=tg)
    sc = scorecard(df, kind=kind, tg=tg)
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    years = sorted(st.year.unique()) if len(st) else sorted(agg.year.unique())
    ncol = max(len(years), 1)
    # Collapse to a single row when the age bands are unrecoverable (016, 017),
    # rather than leaving two thirds of the canvas as a placeholder.
    has_bands = len(st) > 0
    fig = plt.figure(figsize=(4.6 * ncol, 10.2 if has_bands else 4.4))
    gs = fig.add_gridspec(3 if has_bands else 1, ncol,
                          height_ratios=[1, 1, 0.85] if has_bands else [1],
                          hspace=0.35)
    row = 2 if has_bands else 0

    band_lbl = ("5-95% of draws" if kind == "ensemble"
                else "mean +/- 1 SD over seeds")

    # Rows 0-1: age-stratified, women then men
    if len(st):
        ymax = float(np.nanmax([st.model_hi.max(), st.phia_ub.max()])) * 1.05
        if not np.isfinite(ymax) or ymax <= 0:
            ymax = 1.0
        for i, sex in enumerate(("f", "m")):
            for j, yr in enumerate(years):
                ax = fig.add_subplot(gs[i, j])
                s = st[(st.year == yr) & (st.sex == sex)].sort_values("age_low")
                if len(s):
                    ax.plot(s.age_low, s.model, "-o", ms=3, color="#4c78a8",
                            label="model")
                    ax.fill_between(s.age_low, s.model_lo, s.model_hi,
                                    alpha=0.2, color="#4c78a8")
                    ax.errorbar(s.age_low, s.phia,
                                yerr=[s.phia - s.phia_lb, s.phia_ub - s.phia],
                                fmt="ko", ms=4, capsize=2, label="PHIA")
                ax.set(title=f"{SEX_LABEL[sex]} {yr}", ylim=(0, ymax))
                ax.grid(alpha=0.3)
                if j == 0:
                    ax.set_ylabel("HIV prevalence")
                if i == 1:
                    ax.set_xlabel("age")
                if i == 0 and j == 0:
                    ax.legend(fontsize=8)

    # Last row: the 15-49 subpanel plus the scorecard. The only row when the age
    # bands are unrecoverable.
    ax = fig.add_subplot(gs[row, 0:max(ncol - 1, 1)])
    if not len(agg):
        ax.text(0.5, 0.5, "prevalence 15-49 not recoverable", ha="center",
                va="center", fontsize=10)
        ax.axis("off")
    for sex, colr in (("f", "#e45756"), ("m", "#4c78a8")) if len(agg) else ():
        s = agg[agg.sex == sex].sort_values("year")
        if not len(s):
            continue
        ax.errorbar(s.year, s.model, yerr=s.model_sd, fmt="-o", ms=5,
                    color=colr, capsize=3, label=f"model {SEX_LABEL[sex]}")
        ax.plot(s.year, s.phia, "s--", ms=7, mfc="none", color=colr,
                label=f"PHIA {SEX_LABEL[sex]}")
    ax.set(title="Prevalence 15-49 by sex (this resolution exists for every "
                 "experiment)", xlabel="year", ylabel="prevalence 15-49")
    # Origin at zero: with an auto-scaled axis a 4-point gap on a 0.30 baseline
    # fills half the panel and reads as a catastrophic miss. Anchoring at 0 shows
    # the gap against the quantity's actual magnitude.
    if len(agg):
        top = float(np.nanmax([agg.model.max() + agg.model_sd.fillna(0).max(),
                               agg.phia.max()])) * 1.15
        if np.isfinite(top) and top > 0:
            ax.set_ylim(0, top)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    if ncol > 1:
        ax = fig.add_subplot(gs[row, ncol - 1])
        ax.axis("off")
        lines = [f"{label}", ""]
        if not has_bands:
            lines += ["Age-stratified view unavailable:",
                      "no per-band infected counts were",
                      "saved (016 and 017 only).", ""]
        if "mae" in sc:
            lines += [f"Age-stratified ({sc['n_strata']} strata)",
                      f"  MAE   {sc['mae']:.4f}",
                      f"  bias  {sc['bias']:+.4f}",
                      f"  within PHIA CI  {sc['n_within_ci']}/{sc['n_strata']}",
                      ""]
        if "mae_15_49" in sc:
            lines += ["Aggregated to 15-49",
                      f"  MAE   {sc['mae_15_49']:.4f}",
                      f"  bias  {sc['bias_15_49']:+.4f}",
                      "",
                      "NB the 15-49 figures flatter",
                      "the model: aggregating averages",
                      "young under-infection against",
                      "older over-infection."]
        ax.text(0.0, 0.98, "\n".join(lines), va="top", ha="left", fontsize=9,
                family="monospace")

    title = f"{label} — prevalence vs PHIA     band = {band_lbl}"
    if stamp:
        title += f"\n({stamp})"
    # Push the suptitle clear of the axes title in the collapsed one-row layout
    fig.suptitle(title, y=1.06 if not has_bands else 0.995, fontsize=12)
    fig.savefig(outpath, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return sc


# --- Incidence -----------------------------------------------------------------
#
# Added 2026-09-03 for exp 024. Shares this module with the prevalence figure for
# the same reason: one implementation, so the format cannot drift between
# experiments the way `dashboard_fit_*.png` did.
#
# TWO conventions collide here and both are load-bearing:
#   * incidence is a PERCENT per year (per 100 person-years), not a proportion,
#     because that is how SHIMS publishes it;
#   * the denominator is SUSCEPTIBLES (alive minus infected), not the whole
#     population. Using the whole population understates incidence by roughly
#     (1 - prevalence), which at Eswatini's prevalence is a ~30% error.

INC_TARGETS = "calibration_data/incidence_by_age_sex.csv"
INC_CUTOFF = 2017.0   # last year worth plotting: just past the 2016 target


def load_inc_targets(path=None):
    """Fitted incidence targets. Excludes the 2021 validation hold-out by design.

    `uninformative` marks rows whose published CI reaches zero -- four of eight,
    including all three male 2016 bands. They are plotted hollow so a reader does
    not mistake a point estimate that carries no information for a real miss.
    """
    root = Path(__file__).parent
    return pd.read_csv(root / (path or INC_TARGETS))


def _inc_series(df, sex, lo, hi):
    """Annual incidence % over time for one sex and age range, over susceptibles."""
    lo, hi = (lo // 5) * 5, -(-hi // 5) * 5
    new = alv = inf = None
    for b in range(lo, hi, 5):
        nc = f"popagesex.new_infections_{sex}_{b}_{b + 5}"
        if nc not in df.columns:
            return None
        new = df[nc] if new is None else new + df[nc]
        alv = df[f"popagesex.n_alive_{sex}_{b}_{b + 5}"] if alv is None \
            else alv + df[f"popagesex.n_alive_{sex}_{b}_{b + 5}"]
        inf = df[f"popagesex.n_infected_{sex}_{b}_{b + 5}"] if inf is None \
            else inf + df[f"popagesex.n_infected_{sex}_{b}_{b + 5}"]
    susc = alv - inf
    return 100.0 * new / susc.where(susc > 0)
AGG_INC = "data/eswatini_ppdv.csv"


def load_inc_aggregates(path=None):
    """Published 15-49 incidence aggregates, for the time-series panel.

    NOT fitting targets -- `incidence_by_age_sex.csv` carries the age bands, and
    registering the aggregate as well would double-count the same survey. They
    are exactly what a reader needs to judge the time-series curve, though, so
    they are plotted and labelled as a visual check.

    2011 is 18-49 (SHIMS1's published range), 2016 is 15-49. The 2021 row in this
    file is the validation hold-out and is dropped here.
    """
    root = Path(__file__).parent
    d = pd.read_csv(root / (path or AGG_INC)).dropna(subset=["Inc"])
    return d[d.Year.isin((2011, 2016))][
        ["Year", "AgeCat", "sex", "Inc", "Inc_lb", "Inc_ub"]]


def _band_stat(df, sex, lo, hi, year, kind):
    """(centre, lo, hi) of modelled incidence over one age band in one year."""
    s = _inc_series(df[np.floor(df.timevec) == year], sex, lo, hi)
    v = s.dropna() if s is not None else pd.Series(dtype=float)
    if not len(v):
        return np.nan, np.nan, np.nan
    if kind == "ensemble":
        return v.median(), v.quantile(0.05), v.quantile(0.95)
    sd = v.std(ddof=1)
    sd = 0.0 if not np.isfinite(sd) else sd
    return v.mean(), v.mean() - sd, v.mean() + sd


def plot_incidence_fit(df, label, outpath, kind="arm", tg=None, stamp=None,
                       cutoff=INC_CUTOFF):
    """Incidence against SHIMS by sex: a time series plus per-sex 2016 age profiles.

    Design notes, because the first version of this figure was unreadable
    ---------------------------------------------------------------------
    The 2016 age profile originally put both sexes on one axis, which made men's
    and women's markers land on the same band midpoints and left the reader
    guessing which glyph was model and which was data -- sex was encoded as
    colour on *both*, and filled-vs-hollow carried a third meaning on top.

    So: **colour is the model, black is the data, and the two sexes get their own
    panels.** The modelled band average is drawn as a horizontal segment spanning
    the band it actually averages over, rather than a point at the midpoint,
    because that is what the published estimate is.

    `cutoff` truncates the modelled curve just past the last calibration target.
    Plotting to 2026 invites reading the post-target years as a projection, which
    they are not -- nothing after 2016 constrains them, and the 2021 incidence
    hold-out is deliberately untouched.
    """
    tg = load_inc_targets() if tg is None else tg
    agg = load_inc_aggregates()
    d = df[df.timevec <= cutoff]
    cols = {"f": "#c0392b", "m": "#2c6fbb"}
    names = {"f": "Women", "m": "Men"}

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6),
                             gridspec_kw=dict(width_ratios=[1.75, 1, 1]))

    # --- Left: incidence over time, 15-49, by sex --------------------------------
    ax = axes[0]
    for sex in ("f", "m"):
        s = _inc_series(d, sex, 15, 50)
        if s is None:
            continue
        gg = pd.DataFrame({"t": d.timevec.values, "v": s.values}).dropna()
        if not len(gg):
            continue
        a = gg.groupby("t").v
        if kind == "ensemble":
            centre, lo_b, hi_b = a.median(), a.quantile(0.05), a.quantile(0.95)
        else:
            centre = a.mean()
            sd = a.std(ddof=1).fillna(0.0)
            lo_b, hi_b = centre - sd, centre + sd
        ax.plot(centre.index, centre.values, color=cols[sex], lw=2,
                label=f"model — {names[sex]}", zorder=3)
        if (hi_b - lo_b).abs().sum() > 0:
            ax.fill_between(centre.index, lo_b.values, hi_b.values,
                            color=cols[sex], alpha=0.15, lw=0, zorder=1)
    for _, r in agg.iterrows():
        ax.errorbar(r.Year, r.Inc,
                    yerr=[[r.Inc - r.Inc_lb], [r.Inc_ub - r.Inc]],
                    fmt="s", mfc=cols[r.sex], mec="black", mew=1.3, ms=9,
                    ecolor="black", elinewidth=1.4, capsize=4, zorder=6,
                    label="_nolegend_")
        ax.annotate(f"{r.Inc:.2f}", (r.Year, r.Inc), textcoords="offset points",
                    xytext=(9, -3), fontsize=7.5, color="black")
    ax.errorbar([], [], fmt="s", mfc="white", mec="black", mew=1.3, ms=9,
                ecolor="black", capsize=4, label="SHIMS estimate (95% CI)")
    ax.axvline(cutoff, ls=":", c="grey", lw=1.2)
    ax.annotate(f"curve stops at {cutoff:.0f}", (cutoff, 0),
                textcoords="offset points", xytext=(-5, 12), ha="right",
                fontsize=7.5, color="grey", rotation=90)
    ax.set_xlabel("year"); ax.set_ylabel("HIV incidence (% per year)")
    ax.set_title("Incidence 15-49 by sex\n2011 point is SHIMS1 18-49; 2016 is "
                 "SHIMS2 15-49", fontsize=9.5)
    ax.legend(fontsize=8, loc="upper right"); ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)

    # --- Right: the 2016 age profile, one panel per sex --------------------------
    bands = tg[tg.year == 2016]
    ymax = 0.0
    for k, sex in enumerate(("f", "m")):
        ax = axes[k + 1]
        b = bands[bands.sex == sex].sort_values("age_low")
        for _, t in b.iterrows():
            lo, hi = int(t.age_low), int(t.age_high)
            c, blo, bhi = _band_stat(df, sex, lo, hi, 2016, kind)
            # the model average drawn across the band it averages over
            ax.hlines(c, lo, hi, color=cols[sex], lw=2.6, zorder=4)
            if np.isfinite(blo) and bhi > blo:
                ax.add_patch(plt.Rectangle((lo, blo), hi - lo, bhi - blo,
                                           color=cols[sex], alpha=0.15, lw=0,
                                           zorder=1))
                ymax = max(ymax, bhi)
            mid = (lo + hi) / 2
            ax.errorbar(mid, t.incidence_pct,
                        yerr=[[t.incidence_pct - t.lb], [t.ub - t.incidence_pct]],
                        fmt="o", mfc="white" if t.uninformative else "black",
                        mec="black", mew=1.4, ms=8, ecolor="black",
                        elinewidth=1.4, capsize=4, zorder=6)
            ymax = max(ymax, t.ub)
        ax.hlines([], [], [], color=cols[sex], lw=2.6,
                  label=f"model — {names[sex]} (band average)")
        ax.errorbar([], [], fmt="o", mfc="black", mec="black", ms=8,
                    ecolor="black", capsize=4, label="SHIMS2 (95% CI)")
        if b.uninformative.any():
            ax.errorbar([], [], fmt="o", mfc="white", mec="black", mew=1.4,
                        ms=8, ecolor="black", capsize=4,
                        label="hollow: CI reaches 0, no information")
        ax.set_xlabel("age"); ax.set_xlim(13, 52)
        ax.set_title(f"{names[sex]} — 2016 age profile", fontsize=9.5)
        ax.legend(fontsize=7.5, loc="upper right"); ax.grid(alpha=0.3)
        if k == 0:
            ax.set_ylabel("HIV incidence (% per year)")
    for ax in axes[1:]:
        ax.set_ylim(0, ymax * 1.12)

    ttl = label if stamp is None else f"{label}\n{stamp}"
    fig.suptitle(ttl + "   —   incidence vs SHIMS   (colour = model, black = data; "
                 "2021 held out)", y=1.04, fontsize=11)
    fig.savefig(outpath, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_incidence_age_profile(df, label, outpath, years=(2011, 2016),
                               kind="arm", age_max=65, tg=None, stamp=None):
    """Modelled incidence across the full age range, by sex, one panel per year.

    The banded figure (`plot_incidence_fit`) can only show incidence at the three
    coarse bands SHIMS publishes, which hides the shape. This shows the model's
    own 5-year resolution over 15-`age_max`, which is where the age profile is
    actually legible.

    2011 has no published age-stratified incidence at all -- SHIMS1 reports only
    an 18-49 aggregate -- so that panel is model-only, and saying so is part of
    the point. Where SHIMS2 bands exist (2016) they are overlaid as grey
    segments spanning the band they average over, for reference only.
    """
    tg = load_inc_targets() if tg is None else tg
    cols = {"f": "#c0392b", "m": "#2c6fbb"}
    names = {"f": "Women", "m": "Men"}
    edges = list(range(15, age_max, 5))

    fig, axes = plt.subplots(1, len(years), figsize=(6.6 * len(years), 4.6),
                             squeeze=False)
    axes = axes[0]
    ymax = 0.0
    for ax, year in zip(axes, years):
        for sex in ("f", "m"):
            mid, c, lo_b, hi_b = [], [], [], []
            for a in edges:
                cc, l, h = _band_stat(df, sex, a, a + 5, year, kind)
                mid.append(a + 2.5); c.append(cc); lo_b.append(l); hi_b.append(h)
            ax.plot(mid, c, "-o", color=cols[sex], lw=2, ms=4,
                    label=f"model — {names[sex]}", zorder=4)
            ax.fill_between(mid, lo_b, hi_b, color=cols[sex], alpha=0.15, lw=0)
            ymax = max(ymax, np.nanmax(np.asarray(hi_b, dtype=float)))
        # An "age band" spanning 25+ years is an AGGREGATE, not a band. SHIMS1's
        # only 2011 row is 18-49, and drawing it like a band would invite reading
        # a whole-adult average as a measured age profile.
        b = tg[tg.year == year]
        span = b.age_high - b.age_low
        banded, aggregate = b[span < 25], b[span >= 25]
        for _, t in banded.iterrows():
            ax.hlines(t.incidence_pct, t.age_low, min(t.age_high, age_max),
                      color="black", lw=2.2,
                      ls="--" if t.uninformative else "-", zorder=6)
            ymax = max(ymax, t.incidence_pct)
        for _, t in aggregate.iterrows():
            ax.hlines(t.incidence_pct, t.age_low, min(t.age_high, age_max),
                      color="grey", lw=1.4, ls=(0, (1, 2)), zorder=5)
            ax.annotate(f"{names[t.sex]} {int(t.age_low)}-{int(t.age_high) - 1} "
                        f"aggregate = {t.incidence_pct:.2f}",
                        (min(t.age_high, age_max), t.incidence_pct),
                        textcoords="offset points", xytext=(-4, 4), ha="right",
                        fontsize=7, color="grey")
            ymax = max(ymax, t.incidence_pct)
        if len(banded):
            ax.hlines([], [], [], color="black", lw=2.2,
                      label="SHIMS2 band estimate")
            if banded.uninformative.any():
                ax.hlines([], [], [], color="black", lw=2.2, ls="--",
                          label="SHIMS2 band, CI reaches 0")
        if len(aggregate):
            ax.hlines([], [], [], color="grey", lw=1.4, ls=(0, (1, 2)),
                      label="SHIMS1 aggregate (no age detail published)")
        if not len(banded):
            ax.annotate("no published age-stratified incidence this year —\n"
                        "the model profile is unconstrained by data",
                        (0.5, 0.055), xycoords="axes fraction", ha="center",
                        fontsize=8, color="grey")
        ax.set_title(f"{year}", fontsize=11)
        ax.set_xlabel("age"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        ax.set_xlim(13, age_max + 2)
    axes[0].set_ylabel("HIV incidence (% per year)")
    for ax in axes:
        ax.set_ylim(0, ymax * 1.1)
    ttl = label if stamp is None else f"{label}\n{stamp}"
    fig.suptitle(ttl + "   —   modelled incidence by age and sex, 5-year bands",
                 y=1.03, fontsize=11)
    fig.savefig(outpath, dpi=130, bbox_inches="tight")
    plt.close(fig)
