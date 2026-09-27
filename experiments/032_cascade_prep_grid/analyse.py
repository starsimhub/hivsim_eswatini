"""Exp 032 -- the cascade x PrEP surface.

Every cell is placed on the x axis at the viral suppression it ACHIEVED
(VLS among all PLHIV, 15+, 2030), not the one it was asked for. That is the
whole point: 031's `art_95` asked for 0.95 and delivered 0.818 in men 25-34,
and no amount of relabelling fixes an arm whose name is not true of it.

Outputs
  outputs/grid.csv          one row per cell: x, y, and both z's, with seed spread
  outputs/per_seed.csv      per-cell per-seed cumulative infections
  outputs/substitution.csv  how much PrEP buys the same as a step up the cascade
  figures/surface.png       the two panels -- infections averted, % incidence drop
  figures/surface_heatmap.png   the same in (x, y, z) form

Uncertainty is seed noise only. One parameter point; never a credible interval.
"""

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))

from cascade_analysis import aggregate, conditionals, INK, MUTED  # noqa: E402

OUT, FIG = HERE / "outputs", HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)

WINDOW = (2026, 2040)
X_YEAR = 2030          # the year the x axis is measured at -- state it everywhere
BASE_CELL = ("C0_baseline", "P0_none")
REF_031_BASELINE = 66178      # 031's baseline cum infections, for the stop check

CASC_LAB = {"C0_baseline": "baseline", "C1_linkage": "+ link all diagnosed",
            "C2_suppression": "+ VLS to 0.99", "C3_testing_x2": "+ testing x2",
            "C4_testing_x3": "+ testing x3", "C5_testing_x4": "+ testing x4"}
PREP_LAB = {"P0_none": "no PrEP", "P1_fsw": "FSW 60%",
            "P2_agyw_risk": "+ higher-risk AGYW", "P3_agyw_all": "+ all AGYW",
            "P4_women_25_34": "+ women 25-34"}
# Sequential ramp for the nested PrEP ladder -- ordered, not categorical, so a
# single-hue progression is the honest encoding.
PREP_COL = {"P0_none": "#2c3e50", "P1_fsw": "#1a6fa8", "P2_agyw_risk": "#3d94c7",
            "P3_agyw_all": "#7bb8dd", "P4_women_25_34": "#b8d9ec"}


def load():
    files = sorted(glob.glob(str(OUT / "sims" / "*.parquet")))
    if not files:
        sys.exit("no outputs/sims/*.parquet -- run run.py first")
    d = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    print(f"loaded {len(files)} runs: {d.cascade.nunique()} cascade x "
          f"{d.prep.nunique()} PrEP x {d.seed.nunique()} seeds")
    return d


def _band_cols(d, stem):
    out = []
    for c in d.columns:
        if not c.startswith(f"popagesex.{stem}_"):
            continue
        parts = c.split("_")
        if len(parts) >= 3 and parts[-3] in ("f", "m"):
            out.append(c)
    return out


def add_derived(d):
    d = d.copy()
    d["new_inf"] = d[_band_cols(d, "new_infections")].sum(axis=1)
    inf = _band_cols(d, "n_infected")
    alv = [c.replace("n_infected", "n_alive") for c in inf]
    d["n_susc"] = d[alv].sum(axis=1) - d[inf].sum(axis=1)
    d["incidence"] = 100 * d.new_inf / d.n_susc
    return d


def achieved_x(d):
    """VLS among all PLHIV, 15+, at X_YEAR -- the x coordinate of every cell.

    Pooled across seeds numerator-and-denominator before dividing (mean of
    ratios would bias small strata), which is what `conditionals` does.
    """
    rows = []
    for (c, p), g in d.groupby(["cascade", "prep"]):
        cond = conditionals(aggregate(g, {"15+": (15, 200)},
                                      group_cols=("seed",)),
                            keys=("year", "bin"))
        s = cond[cond.year == X_YEAR].set_index("measure").model
        rows.append(dict(cascade=c, prep=p,
                         aware=s["aware_of_all_plhiv"],
                         on_art_given_aware=s["on_art_given_aware"],
                         vls_given_art=s["vls_given_art"],
                         vls_of_plhiv=s["vls_of_all_plhiv"]))
    return pd.DataFrame(rows)


def build_grid(d):
    lo, hi = WINDOW
    w = d[(d.timevec >= lo) & (d.timevec <= hi)]
    per_seed = (w.groupby(["cascade", "prep", "seed"])
                 .agg(cum_inf=("new_inf", "sum"),
                      py_prep=("hiv.n_on_prep", "sum"))
                 .reset_index())
    # Incidence at the end of the horizon, per cell per seed.
    end = (d[d.timevec == hi].groupby(["cascade", "prep", "seed"])
             .incidence.mean().reset_index(name="inc_end"))
    per_seed = per_seed.merge(end, on=["cascade", "prep", "seed"])
    per_seed.to_csv(OUT / "per_seed.csv", index=False)

    bc, bp = BASE_CELL
    base = per_seed[(per_seed.cascade == bc) & (per_seed.prep == bp)] \
        .set_index("seed")

    rows = []
    for (c, p), g in per_seed.groupby(["cascade", "prep"]):
        g = g.set_index("seed")
        # Paired by seed against the single global baseline cell. Pairing is
        # load-bearing: seed spread (SD ~4,400) exceeds several cell effects.
        av = (base.cum_inf - g.cum_inf).dropna()
        dinc = (100 * (base.inc_end - g.inc_end) / base.inc_end).dropna()
        rows.append(dict(
            cascade=c, prep=p,
            cum_inf=g.cum_inf.mean(),
            averted=av.mean(), averted_sd=av.std(ddof=1),
            averted_lo=av.min(), averted_hi=av.max(),
            seeds_positive=int((av > 0).sum()), n_seeds=len(av),
            inc_end=g.inc_end.mean(),
            inc_drop_pct=dinc.mean(), inc_drop_sd=dinc.std(ddof=1),
            py_prep=g.py_prep.mean()))
    grid = pd.DataFrame(rows).merge(achieved_x(d), on=["cascade", "prep"])

    # The PrEP INCREMENT at each cascade rung, paired by seed against the
    # no-PrEP cell on the SAME rung.
    #
    # Dividing PrEP person-years by `averted` (which is measured against the
    # global baseline) would credit PrEP with the cascade's infections too, and
    # reads backwards: it makes FSW PrEP look MORE efficient at high suppression
    # (1.4 py/infection) than at baseline (3.1), when the truth is the reverse.
    # Efficiency must use the increment the PrEP arm itself is responsible for.
    inc, frac = [], []
    for _, r in grid.iterrows():
        cell = per_seed[(per_seed.cascade == r.cascade)
                        & (per_seed.prep == r.prep)].set_index("seed")
        ref = per_seed[(per_seed.cascade == r.cascade)
                       & (per_seed.prep == "P0_none")].set_index("seed")
        dd = (ref.cum_inf - cell.cum_inf).dropna()
        inc.append(dd.mean())
        # Share of the burden REMAINING at that cascade rung that PrEP removes
        # -- the quantity the abstract reports.
        frac.append(100 * dd.mean() / ref.cum_inf.mean())
    grid["prep_increment"] = inc
    grid["prep_pct_of_residual"] = frac
    grid["py_per_prep_averted"] = np.where(
        grid.prep_increment > 0, grid.py_prep / grid.prep_increment, np.nan)
    grid = grid.sort_values(["prep", "cascade"])
    grid.to_csv(OUT / "grid.csv", index=False)
    return grid, per_seed


def substitution(grid):
    """What does adding PrEP buy, against what a step up the cascade buys?

    The decision-relevant quantity the surface exists to produce: at a given
    achieved suppression, moving UP the PrEP ladder averts X; moving RIGHT along
    the cascade ladder averts Y. Reported as both, per step.
    """
    rows = []
    cascs, preps = list(CASC_LAB), list(PREP_LAB)
    for i, c in enumerate(cascs):
        for j, p in enumerate(preps):
            cur = grid[(grid.cascade == c) & (grid.prep == p)]
            if not len(cur):
                continue
            cur = cur.iloc[0]
            up = grid[(grid.cascade == c) & (grid.prep == preps[j + 1])] \
                if j + 1 < len(preps) else None
            right = grid[(grid.cascade == cascs[i + 1]) & (grid.prep == p)] \
                if i + 1 < len(cascs) else None
            rows.append(dict(
                cascade=c, prep=p, vls_of_plhiv=cur.vls_of_plhiv,
                averted=cur.averted,
                gain_from_more_prep=(up.iloc[0].averted - cur.averted
                                     if up is not None and len(up) else np.nan),
                gain_from_more_cascade=(right.iloc[0].averted - cur.averted
                                        if right is not None and len(right)
                                        else np.nan)))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "substitution.csv", index=False)
    return t


def fig_surface(grid):
    """The headline: two panels, x = achieved suppression, series = PrEP ladder.

    Individual cells are drawn as points, NOT smoothed into a fitted surface.
    Two cells at the same x can differ in z because WHO is suppressed matters
    (031 obs 10), and smoothing would hide exactly that.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
    panels = [("averted", "averted_sd",
               f"Cumulative infections averted, {WINDOW[0]}-{WINDOW[1]}"),
              ("inc_drop_pct", "inc_drop_sd",
               f"Reduction in HIV incidence at {WINDOW[1]} (%)")]
    for ax, (zcol, zsd, ylab) in zip(axes, panels):
        for p in PREP_LAB:
            g = grid[grid.prep == p].sort_values("vls_of_plhiv")
            if not len(g):
                continue
            ax.plot(g.vls_of_plhiv, g[zcol], "-o", ms=6, lw=1.8,
                    color=PREP_COL[p], label=PREP_LAB[p], zorder=3)
            ax.fill_between(g.vls_of_plhiv, g[zcol] - g[zsd], g[zcol] + g[zsd],
                            color=PREP_COL[p], alpha=0.13, lw=0, zorder=2)
        ax.axhline(0, color=INK, lw=0.9)
        ax.set_xlabel(f"viral suppression achieved among all PLHIV, 15+, "
                      f"{X_YEAR}")
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.28)
    axes[0].legend(fontsize=8.5, title="PrEP coverage", title_fontsize=9,
                   frameon=False, loc="upper left")
    fig.suptitle("What prevention adds at every level of treatment scale-up "
                 "Eswatini can actually reach", fontsize=13)
    fig.text(0.5, -0.03,
             "Each point is one scenario, placed at the suppression it ACHIEVED "
             "rather than the one it was asked for. Moving RIGHT is a more "
             "complete cascade\n(link everyone diagnosed, then raise suppression "
             "to 0.99, then scale testing); moving UP a line is wider PrEP. "
             "Bands are +/-1 SD over 10 seeds.",
             ha="center", fontsize=8.2, color=MUTED)
    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    fig.savefig(FIG / "surface.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_heatmap(grid):
    """The same data as (x, y, z), matching the grid's own geometry."""
    cascs, preps = list(CASC_LAB), list(PREP_LAB)
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    for ax, (zcol, title, fmt) in zip(axes, [
            ("averted", f"Infections averted, {WINDOW[0]}-{WINDOW[1]}", "{:,.0f}"),
            ("inc_drop_pct", f"Incidence reduction at {WINDOW[1]} (%)", "{:.0f}%")]):
        M = np.full((len(preps), len(cascs)), np.nan)
        for i, p in enumerate(preps):
            for j, c in enumerate(cascs):
                g = grid[(grid.cascade == c) & (grid.prep == p)]
                if len(g):
                    M[i, j] = g.iloc[0][zcol]
        im = ax.imshow(M, cmap="YlGnBu", origin="lower", aspect="auto")
        for i in range(len(preps)):
            for j in range(len(cascs)):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, fmt.format(M[i, j]), ha="center",
                            va="center", fontsize=7.5,
                            color="white" if M[i, j] > np.nanmax(M) * 0.6
                            else INK)
        xt = [f"{CASC_LAB[c]}\n{grid[(grid.cascade==c)&(grid.prep=='P0_none')].vls_of_plhiv.iloc[0]:.3f}"
              if len(grid[(grid.cascade==c)&(grid.prep=='P0_none')]) else CASC_LAB[c]
              for c in cascs]
        ax.set_xticks(range(len(cascs))); ax.set_xticklabels(xt, fontsize=7.5)
        ax.set_yticks(range(len(preps)))
        ax.set_yticklabels([PREP_LAB[p] for p in preps], fontsize=8.5)
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.035)
    fig.suptitle("Cascade rung (with the suppression it achieves) x PrEP rung",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG / "surface_heatmap.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    d = add_derived(load())
    pd.set_option("display.width", 220)
    grid, per_seed = build_grid(d)

    bc, bp = BASE_CELL
    b = per_seed[(per_seed.cascade == bc) & (per_seed.prep == bp)].cum_inf
    print("\n=== STOP CHECK: does C0/P0 reproduce 031's baseline? ===")
    print(f"  032 C0/P0: {b.mean():,.0f} (seed range {b.min():,.0f}-"
          f"{b.max():,.0f})   031: {REF_031_BASELINE:,.0f}   "
          f"diff {b.mean()-REF_031_BASELINE:+,.0f}")
    if not (b.min() <= REF_031_BASELINE <= b.max()):
        print("  -> OUTSIDE the seed range. The grid harness differs from 031's; "
              "resolve before reporting.")

    print(f"\n=== The x axis: suppression each cascade rung ACHIEVES "
          f"(no PrEP, {X_YEAR}) ===")
    x = grid[grid.prep == "P0_none"].sort_values("vls_of_plhiv")
    print(x[["cascade", "aware", "on_art_given_aware", "vls_given_art",
             "vls_of_plhiv"]].round(3).to_string(index=False))

    print(f"\n=== The grid: infections averted vs baseline, "
          f"{WINDOW[0]}-{WINDOW[1]} ===")
    piv = grid.pivot_table(index="prep", columns="cascade", values="averted")
    print(piv.reindex(index=list(PREP_LAB), columns=list(CASC_LAB))
             .round(0).to_string())

    print(f"\n=== Incidence reduction at {WINDOW[1]} (%) ===")
    piv2 = grid.pivot_table(index="prep", columns="cascade",
                            values="inc_drop_pct")
    print(piv2.reindex(index=list(PREP_LAB), columns=list(CASC_LAB))
              .round(1).to_string())

    print("\n=== What PrEP adds AT each cascade rung (increment over no-PrEP "
          "on the same rung) ===")
    p1 = grid.pivot_table(index="prep", columns="cascade",
                          values="prep_increment")
    print(p1.reindex(index=list(PREP_LAB), columns=list(CASC_LAB))
            .round(0).to_string())
    print("\n    as % of the burden REMAINING at that rung:")
    p2 = grid.pivot_table(index="prep", columns="cascade",
                          values="prep_pct_of_residual")
    print(p2.reindex(index=list(PREP_LAB), columns=list(CASC_LAB))
            .round(1).to_string())

    sub = substitution(grid)
    print("\n=== Substitution: what one more PrEP rung buys vs one more "
          "cascade rung ===")
    s = sub[["cascade", "prep", "vls_of_plhiv", "averted",
             "gain_from_more_prep", "gain_from_more_cascade"]].copy()
    s["vls_of_plhiv"] = s.vls_of_plhiv.round(3)
    print(s.round({"averted": 0, "gain_from_more_prep": 0,
                   "gain_from_more_cascade": 0}).to_string(index=False))

    print("\n=== PrEP efficiency, using the increment PrEP is responsible for "
          "(lower is better) ===")
    eff = (grid[grid.py_prep > 0]
           .sort_values(["prep", "vls_of_plhiv"]))
    e = eff[["cascade", "prep", "vls_of_plhiv", "prep_increment", "py_prep",
             "py_per_prep_averted"]].copy()
    e["vls_of_plhiv"] = e.vls_of_plhiv.round(3)
    print(e.round({"prep_increment": 0, "py_prep": 0,
                   "py_per_prep_averted": 1}).to_string(index=False))

    fig_surface(grid)
    fig_heatmap(grid)
    print("\nwrote outputs/{grid,per_seed,substitution}.csv")
    print("      figures/{surface,surface_heatmap}.png")


if __name__ == "__main__":
    main()
