"""Retrofit the standard prevalence-fit figure onto past experiments.

Writes, for every experiment whose trajectories survive on disk:

  experiments/NNN_*/figures/prevalence_fit_vs_phia.png   the standard figure
  analysis/fit_progression.png                           cross-experiment view
  analysis/fit_scorecard.csv                             MAE / bias, both resolutions

Run: `python plot_fit_progression.py`

The figure comes from `standard_figures.plot_prevalence_fit` -- the same
function every experiment's `run.py` calls -- so a retrofit and a freshly
produced figure are identical by construction. This module only decides *which*
configuration of each experiment to plot.

Coverage and its limits
-----------------------
`experiments/*/outputs/` is gitignored, so this works only where outputs survive
locally. 001-007 have none and cannot be recovered; re-running them would not be
faithful anyway, since they predate the stack in use now.

016 and 017 get the 15-49 panel only: their local `PopByAgeSex` stored `n_alive`
per band but only aggregate infected, so per-band prevalence is unrecoverable.

What is and is not comparable
-----------------------------
009 and 014 are **prior ensembles** -- 50 draws over wide parameter space,
including regions with no epidemic -- drawn as 5-95% envelopes. 016 onward are
**fixed-parameter arms** at one high-transmission point, drawn as mean +/- SD.
A single "progress" axis across both would partly measure the change of object,
so the progression chart marks which is which, and flags the three arms that are
genuinely like-for-like (018, 019, 021 -- same parameters, N and seed count).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from standard_figures import (load_targets, plot_prevalence_fit, fit_by_stratum,
                              fit_15_49)

REPO = Path(__file__).resolve().parent
EXP = REPO / "experiments"
ANALYSIS = REPO / "analysis"
ANALYSIS.mkdir(exist_ok=True)
STAMP = "retrofit: plot_fit_progression.py"

# Each entry names the experiment's *adopted or headline* configuration, so the
# series compares like with like as far as the data allows.
SPECS = [
    dict(exp="009_coverage_check", file="ensemble.parquet", kind="ensemble",
         key="par_idx", select=None, label="009 first coverage check"),
    dict(exp="014_prior_expansion", file="ensemble.parquet", kind="ensemble",
         key="par_idx", select=None, label="014 prior expansion"),
    dict(exp="016_double_counted_mortality", file="results.parquet",
         key="arm", select="hiv_deleted", pset="high_transmission",
         label="016 HIV-deleted mortality"),
    dict(exp="017_version_bump", file="results.parquet",
         key="arm", select="C_1_5_11_hiv_deleted", pset="high_transmission",
         label="017 stisim 1.5.11"),
    dict(exp="018_adopt_and_size", file="adopt.parquet",
         key="tag", select="high_transmission",
         label="018 PrEP removed (model-v1.1)", lfl=True),
    dict(exp="019_age_dependent_survival", file="results.parquet",
         key="arm", select="A_flat_13",
         label="019 arm A (adopted nothing)", lfl=True),
    dict(exp="020_model_sizing", file="size.parquet",
         key="tag", select="n20000_high_transmission",
         label="020 N=20,000"),
    dict(exp="021_vls_input", file="results.parquet",
         key="arm", select="B_vls_phia",
         label="021 VLS input (model-v1.2)", lfl=True),
]


def emit_standard(spec, tg):
    """Write the standard figure for one experiment, via the shared function."""
    path = EXP / spec["exp"] / "outputs" / spec["file"]
    if not path.exists():
        print(f"  skip {spec['exp']}: {spec['file']} not on disk")
        return None
    df = pd.read_parquet(path)
    if spec.get("select") is not None and spec["key"] in df.columns:
        df = df[df[spec["key"]] == spec["select"]]
    if spec.get("pset") and "pset" in df.columns:
        df = df[df.pset == spec["pset"]]
    if not len(df):
        print(f"  skip {spec['exp']}: '{spec.get('select')}' not found")
        return None

    kind = spec.get("kind", "arm")
    out = EXP / spec["exp"] / "figures" / "prevalence_fit_vs_phia.png"
    sc = plot_prevalence_fit(df, spec["label"], out, kind=kind, tg=tg,
                             stamp=STAMP)
    st = fit_by_stratum(df, kind=kind, tg=tg)
    if len(st):
        st = st.assign(label=spec["label"], exp=spec["exp"], kind=kind)
    agg = fit_15_49(df, tg=tg).assign(label=spec["label"], exp=spec["exp"],
                                      kind=kind)
    print(f"  {spec['exp']}: {sc.get('n_strata', 0)} strata, "
          f"{sc.get('n_15_49', 0)} at 15-49 -> {out.relative_to(REPO)}")
    return st, agg, dict(label=spec["label"], exp=spec["exp"], kind=kind,
                         lfl=bool(spec.get("lfl")), **sc)


def plot_progression(sc, agg, path):
    fig, axes = plt.subplots(1, 3, figsize=(18.5, 5))
    x = np.arange(len(sc))
    short = [l.split()[0] for l in sc.label]
    colors = ["#f58518" if k == "ensemble" else ("#e45756" if lfl else "#4c78a8")
              for k, lfl in zip(sc.kind, sc.lfl)]

    def bars(ax, col, title, ylabel):
        have = sc[col].notna().to_numpy() if col in sc else np.zeros(len(sc), bool)
        ax.bar(x[have], sc[col].to_numpy()[have],
               color=[c for c, h in zip(colors, have) if h])
        ax.set_xticks(x)
        ax.set_xticklabels(short, fontsize=9)
        ax.set(title=title, ylabel=ylabel)
        for xi, v in zip(x[have], sc[col].to_numpy()[have]):
            ax.annotate(f"{v:.3f}", (xi, v), ha="center", va="bottom",
                        fontsize=8)
        ax.grid(alpha=0.3, axis="y")

    bars(axes[0], "mae", "MAE vs PHIA, age-stratified\n"
                         "orange = prior ensemble, red = like-for-like arms",
         "MAE in prevalence")
    bars(axes[1], "mae_15_49", "MAE vs PHIA, aggregated to 15-49\n"
                               "this resolution flatters the model",
         "MAE in prevalence 15-49")

    ax = axes[2]
    for label, g in agg.groupby("label", sort=False):
        for sex, ls in (("f", "-"), ("m", "--")):
            s = g[g.sex == sex].sort_values("year")
            if len(s):
                ax.plot(s.year, s.model, ls, marker="o", ms=3,
                        label=f"{label.split()[0]} {sex}")
    p = agg.drop_duplicates(["year", "sex"])
    for sex, mk in (("f", "ks"), ("m", "k^")):
        s = p[p.sex == sex].sort_values("year")
        ax.plot(s.year, s.phia, mk, ms=9, mfc="none", label=f"PHIA {sex}")
    ax.set(title="Prevalence 15-49 by sex, every experiment",
           xlabel="year", ylabel="prevalence 15-49")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, ncol=3)

    fig.suptitle(f"Prevalence fit across experiments   ({STAMP})", y=1.03)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    tg = load_targets()
    print(f"{len(tg)} PHIA target strata, years {sorted(tg.year.unique())}\n")

    strata, aggs, cards = [], [], []
    for spec in SPECS:
        got = emit_standard(spec, tg)
        if got is None:
            continue
        st, agg, card = got
        if len(st):
            strata.append(st)
        aggs.append(agg)
        cards.append(card)

    if not cards:
        print("nothing to plot")
        return
    sc = pd.DataFrame(cards)
    sc.to_csv(ANALYSIS / "fit_scorecard.csv", index=False)
    cols = [c for c in ("label", "kind", "n_strata", "mae", "bias",
                        "n_within_ci", "mae_15_49", "bias_15_49") if c in sc]
    print("\n=== Fit scorecard ===")
    print(sc[cols].round(4).to_string(index=False))

    lfl = sc[sc.lfl]
    if len(lfl) > 1:
        print("\nLike-for-like arms only (same parameters, N and seed count):")
        print("  age-stratified MAE  "
              + " -> ".join(f"{v:.4f}" for v in lfl.mae))
        print("  15-49 MAE           "
              + " -> ".join(f"{v:.4f}" for v in lfl.mae_15_49))

    agg = pd.concat(aggs, ignore_index=True)
    plot_progression(sc, agg, ANALYSIS / "fit_progression.png")
    print(f"\n-> {(ANALYSIS / 'fit_progression.png').relative_to(REPO)}")


if __name__ == "__main__":
    main()
