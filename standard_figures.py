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
            sd = s.std(ddof=1)
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
        ymax = max(st.model_hi.max(), st.phia_ub.max()) * 1.05
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
    fig.suptitle(title, y=0.995, fontsize=12)
    fig.savefig(outpath, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return sc
