"""Prevalence fit vs PHIA, per experiment and across experiments.

Generates, for every experiment whose trajectories survive on disk:

  experiments/NNN_*/figures/prevalence_fit_vs_phia.png   one per experiment
  analysis/fit_progression.png                           cross-experiment view
  analysis/fit_scorecard.csv                             MAE / bias / in-CI count

Run: `python plot_fit_progression.py`

Retrofit, not re-run
--------------------
These figures are produced *retrospectively* from parquet already on disk. No
sims are re-run and no SUMMARY.md is edited. Each figure is stamped as a
retrofit so it cannot be mistaken for something the experiment itself produced.
`experiments/*/outputs/` is gitignored, so this only works where the outputs
still exist locally -- experiments 001-007 have none and cannot be recovered.

What is and is not comparable
-----------------------------
**009 and 014 are prior ensembles**: 50 draws spanning wide parameter space,
including regions where the epidemic never establishes. They are drawn as 5-95%
envelopes, because a mean over draws that include dead epidemics is not a
meaningful "fit".

**018, 019, 020 and 021 are fixed-parameter arms**: one hand-picked
high-transmission point, 10 seeds. They are drawn as mean +/- SD.

Putting both on a single "progress" axis would partly measure the change of
object rather than model improvement, so `fit_progression.png` separates them
and marks the three arms that are genuinely like-for-like (018 adopted, 019 arm
A, 021 arm B -- same parameters, same N, same seed count).

**016 and 017 are absent from the age-stratified figures.** Their local
`PopByAgeSex` recorded `n_alive` per band but only aggregate infected, so
prevalence by age band cannot be reconstructed. They contribute to the 15-49
by-sex panel only.

Sex coding
----------
`calibration_data/prevalence_by_age_sex.csv` uses **0 = Male, 1 = Female** --
per experiments/008 and confirmed against SDHS 2006-07. This is the *opposite*
of stisim's internal convention. See exp 021's config for the full map of which
file uses which.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
EXP = REPO / "experiments"
ANALYSIS = REPO / "analysis"
ANALYSIS.mkdir(exist_ok=True)

PHIA_SEX = {0: "m", 1: "f"}
STAMP = "retrofit: plot_fit_progression.py"

# kind: 'ensemble' -> 5-95% band over draws;  'arm' -> mean +/- SD over seeds
# cols: 'popagesex' -> n_infected/n_alive;    'hiv' -> hiv.* + hiv_epi.* prevalence
SPECS = [
    dict(exp="009_coverage_check", file="ensemble.parquet", kind="ensemble",
         cols="hiv", key="par_idx", select=None,
         label="009 first coverage check (prior ensemble)"),
    dict(exp="014_prior_expansion", file="ensemble.parquet", kind="ensemble",
         cols="hiv", key="par_idx", select=None,
         label="014 prior expansion (9-par ensemble)"),
    dict(exp="018_adopt_and_size", file="adopt.parquet", kind="arm",
         cols="popagesex", key="tag", select="high_transmission",
         label="018 adopted stack (model-v1.1)", likeforlike=True),
    dict(exp="019_age_dependent_survival", file="results.parquet", kind="arm",
         cols="popagesex", key="arm", select="A_flat_13",
         label="019 arm A (= upstream survival)", likeforlike=True),
    dict(exp="020_model_sizing", file="size.parquet", kind="arm",
         cols="popagesex", key="tag", select="n20000_high_transmission",
         label="020 N=20,000"),
    dict(exp="021_vls_input", file="results.parquet", kind="arm",
         cols="popagesex", key="arm", select="B_vls_phia",
         label="021 VLS input (model-v1.2)", likeforlike=True),
]

# 016/017 lack per-band infected counts; they carry prevalence 15-49 by sex only.
SEX_ONLY = [
    dict(exp="016_double_counted_mortality", file="results.parquet",
         key="arm", select="hiv_deleted", pset="high_transmission",
         label="016 HIV-deleted mortality"),
    dict(exp="017_version_bump", file="results.parquet",
         key="arm", select="C_1_5_11_hiv_deleted", pset="high_transmission",
         label="017 stisim 1.5.11"),
]


def targets():
    t = pd.read_csv(REPO / "calibration_data" / "prevalence_by_age_sex.csv")
    t["sex"] = t.Gender.map(PHIA_SEX)
    t["age_low"] = t["start age"].astype(int)
    return t[["Year", "sex", "age_low", "NationalPrevalence", "lb", "ub"]].rename(
        columns={"Year": "year", "NationalPrevalence": "phia"})


def _prev_series(df, spec, sex, lo):
    """Per-draw or per-seed prevalence for one stratum, or None if unavailable."""
    if spec["cols"] == "popagesex":
        icol = f"popagesex.n_infected_{sex}_{lo}_{lo + 5}"
        acol = f"popagesex.n_alive_{sex}_{lo}_{lo + 5}"
        if icol not in df.columns or acol not in df.columns:
            return None
        return df[icol] / df[acol].replace(0, np.nan)
    # pre-018 convention: hiv.* carries 15-35, hiv_epi.* carries 35-65
    stem = "hiv" if lo < 35 else "hiv_epi"
    col = f"{stem}.prevalence_{sex}_{lo}_{lo + 5}"
    if col not in df.columns:
        return None
    return df[col]


def collect(spec, tg):
    path = EXP / spec["exp"] / "outputs" / spec["file"]
    if not path.exists():
        print(f"  skip {spec['exp']}: {spec['file']} not on disk")
        return None
    df = pd.read_parquet(path)
    if spec["select"] is not None:
        if spec["key"] not in df.columns:
            print(f"  skip {spec['exp']}: no '{spec['key']}' column")
            return None
        df = df[df[spec["key"]] == spec["select"]]
        if not len(df):
            print(f"  skip {spec['exp']}: '{spec['select']}' not found")
            return None

    rows = []
    for _, t in tg.iterrows():
        g = df[df.timevec == t.year]
        if not len(g):
            continue
        s = _prev_series(g, spec, t.sex, int(t.age_low))
        if s is None:
            continue
        s = s.dropna()
        if not len(s):
            continue
        if spec["kind"] == "ensemble":
            centre, lo_b, hi_b = s.median(), s.quantile(0.05), s.quantile(0.95)
        else:
            sd = s.std(ddof=1)
            centre, lo_b, hi_b = s.mean(), s.mean() - sd, s.mean() + sd
        rows.append(dict(year=int(t.year), sex=t.sex, age_low=int(t.age_low),
                         model=centre, model_lo=lo_b, model_hi=hi_b,
                         phia=t.phia, phia_lb=t.lb, phia_ub=t.ub))
    if not rows:
        print(f"  skip {spec['exp']}: no strata resolved")
        return None
    out = pd.DataFrame(rows)
    out["label"], out["exp"] = spec["label"], spec["exp"]
    out["kind"] = spec["kind"]
    return out


def plot_one(d, spec, path):
    years = sorted(d.year.unique())
    fig, axes = plt.subplots(2, len(years), figsize=(4.5 * len(years), 7.4),
                             squeeze=False, sharey=True)
    band_kw = dict(alpha=0.2, color="#4c78a8")
    for j, yr in enumerate(years):
        for i, sex in enumerate(("f", "m")):
            ax = axes[i][j]
            s = d[(d.year == yr) & (d.sex == sex)].sort_values("age_low")
            if len(s):
                ax.plot(s.age_low, s.model, "-o", ms=3, color="#4c78a8",
                        label="model" if (i == 0 and j == 0) else None)
                ax.fill_between(s.age_low, s.model_lo, s.model_hi, **band_kw)
                ax.errorbar(s.age_low, s.phia,
                            yerr=[s.phia - s.phia_lb, s.phia_ub - s.phia],
                            fmt="ko", ms=4, capsize=2,
                            label="PHIA" if (i == 0 and j == 0) else None)
            ax.set_title(f"{'Women' if sex == 'f' else 'Men'} {yr}", fontsize=10)
            ax.grid(alpha=0.3)
            if j == 0:
                ax.set_ylabel("HIV prevalence")
            if i == 1:
                ax.set_xlabel("age")
            if i == 0 and j == 0:
                ax.legend(fontsize=8)
    band = "5-95% of prior draws" if spec["kind"] == "ensemble" else "mean +/- 1 SD over seeds"
    fig.suptitle(f"{spec['label']} — prevalence vs PHIA\n"
                 f"band = {band}   ({STAMP})", y=1.0, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def scorecard(all_d):
    rows = []
    for label, g in all_d.groupby("label", sort=False):
        err = g.model - g.phia
        inside = ((g.model >= g.phia_lb) & (g.model <= g.phia_ub)).sum()
        rows.append(dict(label=label, exp=g.exp.iloc[0], kind=g.kind.iloc[0],
                         n_strata=len(g), mae=err.abs().mean(), bias=err.mean(),
                         n_within_phia_ci=int(inside),
                         pct_within_ci=100 * inside / len(g)))
    return pd.DataFrame(rows)


def plot_progression(all_d, sc, sex_only, path):
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 4.8))

    ax = axes[0]
    arms = sc[sc.kind == "arm"]
    ens = sc[sc.kind == "ensemble"]
    x = np.arange(len(sc))
    colors = ["#f58518" if k == "ensemble" else "#4c78a8" for k in sc.kind]
    ax.bar(x, sc.mae, color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels([l.split()[0] for l in sc.label], rotation=0)
    ax.set(title="Mean absolute error vs PHIA\n(orange = prior ensemble, not "
                 "comparable to arms)", ylabel="MAE in prevalence")
    for xi, v in zip(x, sc.mae):
        ax.annotate(f"{v:.3f}", (xi, v), ha="center", va="bottom", fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1]
    ax.bar(x, sc.bias, color=colors)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([l.split()[0] for l in sc.label])
    ax.set(title="Bias vs PHIA (negative = model too low)", ylabel="mean error")
    for xi, v in zip(x, sc.bias):
        ax.annotate(f"{v:+.3f}", (xi, v), ha="center",
                    va="top" if v < 0 else "bottom", fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[2]
    for label, g in sex_only.groupby("label", sort=False):
        for sex, ls in (("f", "-"), ("m", "--")):
            s = g[g.sex == sex]
            ax.plot(s.year, s.model, ls, marker="o", ms=3, label=f"{label.split()[0]} {sex}")
    ax.set(title="Prevalence 15-49 by sex, all experiments\n"
                 "(016/017 recoverable only at this resolution)",
           xlabel="year", ylabel="prevalence 15-49")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)

    fig.suptitle(f"Fit progression across experiments   ({STAMP})", y=1.03)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def collect_sex_only():
    """Prevalence 15-49 by sex, including 016/017 which lack per-band counts."""
    rows = []
    specs = SEX_ONLY + [s for s in SPECS if s["cols"] == "popagesex"]
    for spec in specs:
        path = EXP / spec["exp"] / "outputs" / spec["file"]
        if not path.exists():
            continue
        df = pd.read_parquet(path)
        if spec["select"] is not None and spec["key"] in df.columns:
            df = df[df[spec["key"]] == spec["select"]]
        if spec.get("pset") and "pset" in df.columns:
            df = df[df.pset == spec["pset"]]
        if not len(df):
            continue
        for sex in ("f", "m"):
            col = f"popagesex.prevalence_{sex}_15_49"
            if col not in df.columns:
                # post-018 files: rebuild from per-band counts
                icols = [c for c in df.columns
                         if c.startswith(f"popagesex.n_infected_{sex}_")
                         and 15 <= int(c.split("_")[-2]) < 50]
                acols = [c.replace("n_infected", "n_alive") for c in icols]
                if not icols:
                    continue
                ser = df[icols].sum(axis=1) / df[acols].sum(axis=1).replace(0, np.nan)
                tmp = df[["timevec"]].assign(v=ser)
            else:
                tmp = df[["timevec"]].assign(v=df[col])
            m = tmp.groupby("timevec")["v"].mean()
            for yr, v in m.items():
                rows.append(dict(label=spec["label"], sex=sex,
                                 year=float(yr), model=v))
    return pd.DataFrame(rows)


def main():
    tg = targets()
    print(f"{len(tg)} PHIA target strata "
          f"({sorted(tg.year.unique())}, sexes {sorted(tg.sex.unique())})\n")

    collected = []
    for spec in SPECS:
        d = collect(spec, tg)
        if d is None:
            continue
        fig_dir = EXP / spec["exp"] / "figures"
        fig_dir.mkdir(exist_ok=True)
        out = fig_dir / "prevalence_fit_vs_phia.png"
        plot_one(d, spec, out)
        print(f"  {spec['exp']}: {len(d)} strata -> {out.relative_to(REPO)}")
        collected.append(d)

    if not collected:
        print("nothing to plot")
        return
    all_d = pd.concat(collected, ignore_index=True)
    sc = scorecard(all_d)
    sc.to_csv(ANALYSIS / "fit_scorecard.csv", index=False)
    print("\n=== Fit scorecard ===")
    print(sc.round(4).to_string(index=False))

    lfl = sc[sc.label.str.startswith(("018", "019", "021"))]
    if len(lfl) > 1:
        print("\nLike-for-like only (same parameters, N and seed count):")
        print(f"  MAE  {' -> '.join(f'{v:.4f}' for v in lfl.mae)}")
        print(f"  bias {' -> '.join(f'{v:+.4f}' for v in lfl.bias)}")

    sex_only = collect_sex_only()
    plot_progression(all_d, sc, sex_only, ANALYSIS / "fit_progression.png")
    print(f"\n-> {(ANALYSIS / 'fit_progression.png').relative_to(REPO)}")


if __name__ == "__main__":
    main()
