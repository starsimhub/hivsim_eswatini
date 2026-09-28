"""Exp 033 -- how far can the cascade actually go, and what does PrEP add?

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
BASE_CELL = ("S0_status_quo", "P0_none")
REF_031_BASELINE = 66178      # 031/032's baseline, for the stop check

CASC_LAB = {"S0_status_quo": "status quo",
            "S1_testing_only": "testing x3 only",
            "S2_unaids_95": "95-95-95 every group",
            "S3_best_in_class": "best-in-class",
            "S4_bound": "BOUND (not a scenario)"}
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
    fig, axes = plt.subplots(1, 3, figsize=(18.5, 5.4))
    panels = [("averted", "averted_sd",
               f"Cumulative infections averted, {WINDOW[0]}-{WINDOW[1]}"),
              ("inc_drop_pct", "inc_drop_sd",
               f"Reduction in HIV incidence at {WINDOW[1]} (%)"),
              # The panel that answers "does prevention still matter when
              # suppression is high": PrEP's share of the burden REMAINING at
              # each rung, so the denominator shrinks as the cascade improves.
              ("prep_pct_of_residual", None,
               "PrEP's share of the burden REMAINING (%)")]
    for ax, (zcol, zsd, ylab) in zip(axes, panels):
        for p in PREP_LAB:
            g = grid[grid.prep == p].sort_values("vls_of_plhiv")
            if not len(g) or (zcol == "prep_pct_of_residual" and p == "P0_none"):
                continue
            ax.plot(g.vls_of_plhiv, g[zcol], "-o", ms=6, lw=1.8,
                    color=PREP_COL[p], label=PREP_LAB[p], zorder=3)
            if zsd is not None:
                ax.fill_between(g.vls_of_plhiv, g[zcol] - g[zsd],
                                g[zcol] + g[zsd], color=PREP_COL[p],
                                alpha=0.13, lw=0, zorder=2)
        ax.axhline(0, color=INK, lw=0.9)
        ax.set_xlabel(f"viral suppression achieved among all PLHIV, 15+, "
                      f"{X_YEAR}")
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.28)
    axes[2].set_ylim(bottom=0)
    axes[0].legend(fontsize=8.5, title="PrEP coverage", title_fontsize=9,
                   frameon=False, loc="upper left")
    fig.suptitle("What prevention adds at every level of treatment scale-up "
                 "Eswatini can actually reach", fontsize=13)
    fig.text(0.5, -0.045,
             "Each point is one scenario, placed at the suppression it ACHIEVED "
             "rather than the one it was asked for. Moving RIGHT is a more "
             "complete cascade (link everyone diagnosed,\nthen raise suppression "
             "to 0.99, then scale testing); moving UP a line is wider PrEP. The "
             "cascade axis SATURATES near 0.95 -- about 5% of people with HIV "
             "stay unsuppressed whatever is done.\nBands are +/-1 SD across 10 "
             "seeds on each line's own level; comparisons BETWEEN lines share "
             "seeds and are far tighter than the bands suggest.",
             ha="center", fontsize=8.2, color=MUTED)
    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    fig.savefig(FIG / "surface.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_heatmap(grid):
    """The same data as (x, y, z), matching the grid's own geometry.

    Reference row and column excluded -- see fig_xyz for why.
    """
    ref = grid
    grid = grid[grid.cascade != "S0_status_quo"]
    cascs = [c for c in CASC_LAB if c != "S0_status_quo"]
    preps = list(PREP_LAB)          # no-PrEP row kept as the reference
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
        # Achieved suppression comes from the no-PrEP cell of each rung, which
        # is in `ref` -- it was filtered out of the plotted data above.
        xt = []
        for c in cascs:
            r = ref[(ref.cascade == c) & (ref.prep == "P0_none")]
            xt.append(f"{CASC_LAB[c]}\n{r.vls_of_plhiv.iloc[0]:.3f}"
                      if len(r) else CASC_LAB[c])
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


def fig_xyz(grid):
    """The (x, y, z) scatter: suppression achieved, PrEP programme size, impact.

    Adam's design. Unlike the line plots, this puts BOTH interventions on real
    axes and lets the outcome be read as colour, so the shape of the trade-off
    is visible rather than inferred: how far right you can get on treatment,
    how far up you must go on prevention, and what each buys.

    y is the PrEP rung, EVENLY SPACED, with each tick labelled by the programme
    size it represents. A linear person-years axis was tried first and failed:
    FSW delivery is 31k person-years against 885k for the broadest arm, so the
    two smallest rungs collapsed onto each other at the bottom of the plot and
    their labels overlapped. Spacing the rungs evenly and putting the magnitude
    in the tick label keeps both the ordering and the size legible.

    x stays on the true achieved-suppression scale -- that is the whole point of
    the design, and rungs that land close together SHOULD look close together.
    """
    panels = [("averted", "Cumulative infections averted, "
               f"{WINDOW[0]}-{WINDOW[1]}", "{:,.0f}"),
              ("inc_drop_pct", f"Reduction in HIV incidence at {WINDOW[1]} (%)",
               "{:.0f}")]
    import matplotlib.patheffects as pe

    # The no-PrEP row IS kept: it is the counterfactual each cascade rung is
    # read against, and without it PrEP's contribution is the gap between a line
    # and nothing. The baseline-cascade COLUMN is dropped, which also removes
    # the degenerate C0/P0 cell whose "averted" is 0 by construction (each seed
    # differenced against itself). Its values remain in grid.csv.
    ref = grid                                   # keep for the x positions
    grid = grid[grid.cascade != "S0_status_quo"]

    preps = [p for p in PREP_LAB if (grid.prep == p).any()]
    ypos = {p: i for i, p in enumerate(preps)}
    pyk = grid.groupby("prep").py_prep.mean().to_dict()
    ylabels = [PREP_LAB[p] if pyk.get(p, 0) < 1 else
               f"{PREP_LAB[p]}\n({pyk[p]/1000:,.0f}k person-yrs)" for p in preps]

    fig, axes = plt.subplots(1, 2, figsize=(17.5, 6.4))
    for ax, (zcol, title, fmt) in zip(axes, panels):
        g = grid.copy()
        g["y"] = g.prep.map(ypos)
        z = g[zcol]
        # Faint connectors along both ladders, so the grid structure reads
        # behind the points.
        for p in preps:
            s = g[g.prep == p].sort_values("vls_of_plhiv")
            ax.plot(s.vls_of_plhiv, s.y, "-", lw=0.9, color="#cdd5db", zorder=1)
        for c in CASC_LAB:
            s = g[g.cascade == c].sort_values("y")
            if len(s) > 1:
                ax.plot(s.vls_of_plhiv, s.y, "-", lw=0.9, color="#cdd5db",
                        zorder=1)
        rng = max(z.max() - z.min(), 1e-9)
        sizes = 70 + 340 * (z - z.min()) / rng
        sc = ax.scatter(g.vls_of_plhiv, g.y, c=z, s=sizes, cmap="viridis",
                        edgecolor="white", linewidth=0.9, zorder=3)
        # Labels sit ABOVE each point with a white outline, not inside it --
        # inside-the-marker text collided wherever two rungs landed at a similar
        # suppression, which is exactly where the interesting cells are.
        # Several cascade rungs land at almost the same achieved suppression
        # (the testing rungs differ by <0.002), so a fixed label offset makes
        # their numbers overprint. Stagger the offset for any point that sits
        # within xtol of the previous one in the same row.
        xtol = 0.004
        for p in preps:
            s = g[g.prep == p].sort_values("vls_of_plhiv")
            offsets, last_x, k = [], None, 0
            for xv in s.vls_of_plhiv:
                if last_x is not None and abs(xv - last_x) < xtol:
                    k += 1
                else:
                    k = 0
                offsets.append([(0, 13), (0, -19), (0, 25)][k % 3])
                last_x = xv
            for (_, r), off in zip(s.iterrows(), offsets):
                ax.annotate(fmt.format(r[zcol]), (r.vls_of_plhiv, r.y),
                            xytext=off, textcoords="offset points",
                            fontsize=6.8, ha="center",
                            va="bottom" if off[1] > 0 else "top", color=INK,
                            zorder=4,
                            path_effects=[pe.withStroke(linewidth=2.4,
                                                        foreground="white")])
        ax.set_xlabel(f"viral suppression achieved among all PLHIV, 15+, "
                      f"{X_YEAR}  →  more complete cascade")
        ax.set_yticks(range(len(preps)))
        ax.set_yticklabels(ylabels, fontsize=8)
        ax.set_ylim(-0.6, len(preps) - 0.35)
        ax.margins(x=0.09)
        ax.set_title(title, fontsize=10.5)
        ax.grid(alpha=0.25, zorder=0)
        fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    axes[0].set_ylabel("PrEP programme  →  wider coverage")
    fig.suptitle("Treatment and prevention on their own axes: what each "
                 "combination averts", fontsize=13)
    fig.text(0.5, -0.02,
             "Each point is one scenario. RIGHT = a more complete cascade, "
             "placed at the suppression it ACHIEVED rather than the one it was "
             "asked for. UP = a larger PrEP programme,\nexpanding from female "
             "sex workers outward. Colour and size = the outcome. Grey lines "
             "trace the two ladders.",
             ha="center", fontsize=8.2, color=MUTED)
    fig.tight_layout(rect=[0, 0.01, 1, 0.94])
    fig.savefig(FIG / "xyz_scatter.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_prep_by_counterfactual(grid, per_seed):
    """PrEP's contribution under each cascade counterfactual -- the headline.

    Adam's reframing, and it is the better instrument. Reporting cells by TOTAL
    infections averted mixes the cascade's effect with PrEP's, and readers
    misattribute the sum. The increment is also far tighter statistically: it is
    paired by seed WITHIN a cascade rung, so the cascade effect cancels (SDs of
    ~1,100-2,400 against ~3,500+ for total averted). And it does not depend on
    where the cascade ceiling sits -- the question that has already been got
    wrong twice.

    Two panels because they answer different questions: the absolute increment
    is what a budget buys, the share of remaining burden is the fair comparison
    across counterfactuals whose denominators differ.
    """
    import matplotlib.patheffects as pe

    cascs = [c for c in CASC_LAB if (grid.cascade == c).any()]
    preps = [p for p in PREP_LAB if p != "P0_none"]
    xpos = {c: i for i, c in enumerate(cascs)}

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.8))
    panels = [("prep_increment", "Infections averted by PrEP, "
               f"{WINDOW[0]}-{WINDOW[1]}", "{:,.0f}"),
              ("prep_pct_of_residual",
               "PrEP's share of the burden remaining (%)", "{:.1f}")]
    for ax, (zcol, ylab, fmt) in zip(axes, panels):
        # The bound is not a scenario; shade it out rather than dropping it, so
        # the plausible rungs can be read against it without being confused for
        # it.
        if "S4_bound" in xpos:
            ax.axvspan(xpos["S4_bound"] - 0.5, xpos["S4_bound"] + 0.5,
                       color="#eceff1", zorder=0)
        for p in preps:
            g = grid[grid.prep == p].copy()
            g["x"] = g.cascade.map(xpos)
            g = g.sort_values("x")
            ax.plot(g.x, g[zcol], "-o", ms=7, lw=2, color=PREP_COL[p],
                    label=PREP_LAB[p], zorder=3)
            for _, r in g.iterrows():
                ax.annotate(fmt.format(r[zcol]), (r.x, r[zcol]),
                            xytext=(0, 9), textcoords="offset points",
                            fontsize=7, ha="center", color=INK, zorder=4,
                            path_effects=[pe.withStroke(linewidth=2.4,
                                                        foreground="white")])
        # Name each counterfactual by what it is AND what suppression it reaches
        labs = []
        for c in cascs:
            r = grid[(grid.cascade == c) & (grid.prep == "P0_none")]
            v = f"\n{r.vls_of_plhiv.iloc[0]:.3f}" if len(r) else ""
            labs.append(CASC_LAB[c].replace(" every", "\nevery") + v)
        ax.set_xticks(range(len(cascs)))
        ax.set_xticklabels(labs, fontsize=8)
        ax.set_ylabel(ylab)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.28)
        ax.margins(x=0.08)
    axes[0].legend(fontsize=8.5, title="PrEP programme", title_fontsize=9,
                   frameon=False)
    fig.suptitle("What prevention adds under each treatment-cascade "
                 "counterfactual", fontsize=13)
    fig.text(0.5, -0.035,
             "Each point is PrEP's increment over the SAME cascade with no "
             "PrEP, paired by seed. x labels give the scenario and the viral "
             "suppression among all PLHIV it reaches.\nThe shaded rung is a "
             "BOUND, not a scenario: it requires every age-sex group to exceed "
             "the best-performing group Eswatini has ever measured.",
             ha="center", fontsize=8.2, color=MUTED)
    fig.tight_layout(rect=[0, 0.01, 1, 0.94])
    fig.savefig(FIG / "prep_by_counterfactual.png", dpi=140,
                bbox_inches="tight")
    plt.close(fig)


def attribution(per_seed, out=True):
    """Split each scenario's infections averted into cascade- and PrEP-attributable.

    The two levers are partially substitutable, so their joint effect is
    SUB-ADDITIVE: cascade-alone plus PrEP-alone exceeds the two together, by an
    overlap of up to 10,183 infections here. That overlap has to be allocated,
    and the choice is not innocent.

    `sequential` gives the whole overlap to the cascade, purely because it is
    named first -- at the bound that reads 88/12 rather than 75/25.

    `shapley` averages over both orderings, which for two players is
    cascade = A + I/2, prep = B + I/2. It is symmetric, sums exactly to the
    joint effect, and does not depend on an arbitrary ordering. Preferred;
    sequential is kept as a sensitivity.
    """
    def cell(c, p):
        return per_seed[(per_seed.cascade == c)
                        & (per_seed.prep == p)].set_index("seed").cum_inf
    bc, bp = BASE_CELL
    base = cell(bc, bp)
    rows = []
    for c in CASC_LAB:
        for p in PREP_LAB:
            if not len(cell(c, p)):
                continue
            A = (base - cell(c, bp)).mean()      # cascade alone
            B = (base - cell(bc, p)).mean()      # PrEP alone
            J = (base - cell(c, p)).mean()       # both
            I = J - A - B
            rows.append(dict(
                cascade=c, prep=p, joint_averted=J, cascade_solo=A,
                prep_solo=B, interaction=I,
                seq_cascade=A, seq_prep=J - A,
                shap_cascade=A + I / 2, shap_prep=B + I / 2,
                shap_cascade_pct=100 * (A + I / 2) / J if J else np.nan,
                shap_prep_pct=100 * (B + I / 2) / J if J else np.nan,
                seq_cascade_pct=100 * A / J if J else np.nan,
                seq_prep_pct=100 * (J - A) / J if J else np.nan))
    t = pd.DataFrame(rows)
    if out:
        t.to_csv(OUT / "attribution.csv", index=False)
    return t


def fig_attribution(attr, baseline_total):
    """Every bar is the WHOLE epidemic, split by what did and did not prevent it.

    Adam's design. Expressing the split as a share of infections AVERTED makes
    every scenario look complete -- the segments always sum to 100% however
    small the effect. Against the full baseline burden instead, the grey block
    is the fraction neither lever prevents, which is the quantity the paper is
    actually about.

    Segments use Shapley attribution, so the split does not depend on which
    lever is counted first.
    """
    show_prep = ["P0_none", "P1_fsw", "P4_women_25_34"]
    cascs = [c for c in CASC_LAB if c != "S0_status_quo"]
    width, gap = 0.26, 0.02
    fig, ax = plt.subplots(figsize=(14, 6.8))
    if "S4_bound" in cascs:
        i = cascs.index("S4_bound")
        ax.axvspan(i - 0.5, i + 0.5, color="#eceff1", zorder=0)
    for j, p in enumerate(show_prep):
        xs = [i + (j - 1) * (width + gap) for i in range(len(cascs))]
        for x, c in zip(xs, cascs):
            r = attr[(attr.cascade == c) & (attr.prep == p)]
            cv = r.shap_cascade.iloc[0] if len(r) else 0.0
            pv = r.shap_prep.iloc[0] if len(r) else 0.0
            rem = baseline_total - cv - pv
            first = (j == 0 and x == xs[0])
            ax.bar(x, cv, width, color="#2c6fbb", zorder=2,
                   label="averted by the cascade" if first else None)
            ax.bar(x, pv, width, bottom=cv, color="#e8a33d", zorder=2,
                   label="averted by PrEP" if first else None)
            ax.bar(x, rem, width, bottom=cv + pv, color="#d5dade", zorder=2,
                   label="not averted" if first else None)
            for val, bot in ((cv, 0), (pv, cv), (rem, cv + pv)):
                if val > baseline_total * 0.035:
                    ax.text(x, bot + val / 2, f"{100*val/baseline_total:.0f}%",
                            ha="center", va="center", fontsize=7.4, zorder=3,
                            color="white" if val is not rem else INK)
        for x in xs:
            ax.text(x, -baseline_total * 0.035,
                    PREP_LAB[p].replace("+ ", "+\n"), ha="center", va="top",
                    fontsize=6.6, color=MUTED)
    ax.axhline(baseline_total, color=INK, lw=1.4, ls="--", zorder=4)
    ax.text(len(cascs) - 0.42, baseline_total * 1.012,
            f"baseline: {baseline_total:,.0f} infections with no cascade "
            f"improvement and no PrEP", ha="right", fontsize=8, color=INK)
    labs = [CASC_LAB[c].replace(" every", "\nevery") for c in cascs]
    ax.set_xticks(range(len(cascs)))
    ax.set_xticklabels(labs, fontsize=9)
    ax.tick_params(axis="x", pad=34)
    ax.set_ylabel(f"cumulative HIV infections, {WINDOW[0]}-{WINDOW[1]}")
    ax.set_ylim(0, baseline_total * 1.08)
    ax.grid(axis="y", alpha=0.28, zorder=0)
    # Legend above the axes: inside the plot it sat on top of the segment
    # labels in the leftmost group, where the bars are shortest.
    ax.legend(fontsize=9, frameon=False, ncol=3, loc="lower center",
              bbox_to_anchor=(0.5, 1.02))
    ax.set_title("Of every infection Eswatini would otherwise see, how many "
                 "does each lever prevent?", fontsize=12.5, pad=34)
    fig.text(0.5, -0.06,
             "Each bar is the full status-quo burden. Percentages are shares of "
             "ALL baseline infections, not of those averted, so the grey block "
             "is what neither lever prevents.\nThe two levers overlap, so the "
             "split uses Shapley attribution, which shares that overlap evenly "
             "and does not depend on which is counted first. The shaded rung is "
             "a BOUND, not a scenario.",
             ha="center", fontsize=8.2, color=MUTED)
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(FIG / "attribution.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_viremia(d):
    """Prevalence of detectable viral load among adults, by cascade scenario.

    Everything is computed on the 15+ population: the cascade analyzer's bands
    start at 15, so mixing its treated counts with an all-ages PLHIV
    denominator would inflate the viremic fraction by the paediatric
    infections it never sees.

    No PrEP in any curve -- PrEP changes who becomes infected, not how many of
    the infected are suppressed, so holding it at zero isolates the cascade.
    """
    sub = d[d.prep == "P0_none"].copy()
    inf = [c for c in sub.columns if c.startswith("cascadeage.n_infected_")]
    eff = [c for c in sub.columns
           if c.startswith("cascadeage.n_effective_art_")]
    alive = []
    for c in sub.columns:
        if not c.startswith("popagesex.n_alive_"):
            continue
        parts = c.split("_")
        if len(parts) >= 3 and parts[-3] in ("f", "m") and int(parts[-2]) >= 15:
            alive.append(c)
    sub["viremic"] = sub[inf].sum(axis=1) - sub[eff].sum(axis=1)
    sub["adults"] = sub[alive].sum(axis=1)
    sub["viremia_prev"] = 100 * sub.viremic / sub.adults

    cascs = [c for c in CASC_LAB if (sub.cascade == c).any()]
    cmap = plt.get_cmap("viridis")
    cols = {c: cmap(i / max(len(cascs) - 1, 1)) for i, c in enumerate(cascs)}

    fig, ax = plt.subplots(figsize=(9.6, 6.0))
    rows = []
    for c in cascs:
        g = sub[sub.cascade == c].groupby("timevec").viremia_prev
        m, sd = g.mean(), g.std(ddof=1).fillna(0)
        w = (m.index >= 2020) & (m.index <= 2040)
        style = dict(lw=2.4, color=cols[c])
        if c == "S4_bound":
            style.update(ls="--", lw=1.8)
        ax.plot(m.index[w], m.values[w],
                label=CASC_LAB[c].replace("\n", " "), zorder=3, **style)
        ax.fill_between(m.index[w], (m - sd).values[w], (m + sd).values[w],
                        color=cols[c], alpha=0.13, lw=0, zorder=2)
        rows.append(dict(cascade=c, y2026=m.get(2026), y2030=m.get(2030),
                         y2040=m.get(2040)))
    ax.axvline(2026, ls=":", color=MUTED, lw=1)
    ax.text(2026.15, ax.get_ylim()[1] * 0.97, "scenarios begin", fontsize=7.6,
            color=MUTED, va="top")
    ax.set_xlim(2020, 2040)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("year")
    ax.set_ylabel("adults (15+) with detectable viral load (%)")
    ax.set_title("Prevalence of detectable HIV viral load in the adult "
                 "population, by cascade scenario", fontsize=12)
    ax.legend(fontsize=8.5, frameon=False, title="ART cascade (no PrEP)",
              title_fontsize=9)
    ax.grid(alpha=0.28)
    fig.text(0.5, -0.04,
             "Share of all adults 15+ who are living with HIV and not virally "
             "suppressed -- the population reservoir available to transmit. "
             "Bands are +/-1 SD across 10 seeds.\nThe dashed curve is a BOUND, "
             "not a scenario: it requires every age-sex group to exceed the "
             "best-performing group Eswatini has ever measured.",
             ha="center", fontsize=8.2, color=MUTED)
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(FIG / "viremia_prevalence.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "viremia_prevalence.csv", index=False)
    return t


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

    # Where does the cascade actually stop, and who is left? 032 stopped at an
    # ART target of 0.95, called the resulting 0.948 a ceiling, and was wrong:
    # 78% of that residual was people DIAGNOSED BUT NOT ON ART, i.e. the target
    # itself. This table is what caught that, so it is a standing output now.
    x = grid[grid.prep == "P0_none"].copy()
    x["unaware"] = 1 - x.aware
    x["aware_not_on_art"] = x.aware * (1 - x.on_art_given_aware)
    x["on_art_not_suppr"] = (x.aware * x.on_art_given_aware
                             * (1 - x.vls_given_art))
    x["unsuppressed"] = 1 - x.vls_of_plhiv
    x = x.sort_values("vls_of_plhiv")
    x[["cascade", "aware", "on_art_given_aware", "vls_given_art",
       "vls_of_plhiv", "unaware", "aware_not_on_art", "on_art_not_suppr",
       "unsuppressed"]].to_csv(OUT / "residual_decomposition.csv", index=False)
    print("\n=== Who is NOT virally suppressed, by rung (share of all PLHIV) ===")
    print(x[["cascade", "vls_of_plhiv", "unaware", "aware_not_on_art",
             "on_art_not_suppr", "unsuppressed"]].round(4).to_string(index=False))

    attr = attribution(per_seed)
    print("\n=== Infections averted vs baseline, split by lever (Shapley) ===")
    print("    The levers overlap, so cascade-alone + PrEP-alone EXCEEDS the")
    print("    two together. Shapley shares that overlap evenly; 'seq' gives it")
    print("    all to the cascade purely because it is named first.")
    show = attr[attr.prep != "P0_none"].copy()
    print(show[["cascade", "prep", "joint_averted", "cascade_solo", "prep_solo",
                "interaction", "shap_cascade_pct", "shap_prep_pct",
                "seq_prep_pct"]]
          .round({"joint_averted": 0, "cascade_solo": 0, "prep_solo": 0,
                  "interaction": 0, "shap_cascade_pct": 0, "shap_prep_pct": 0,
                  "seq_prep_pct": 0}).to_string(index=False))
    bc, bp = BASE_CELL
    baseline_total = per_seed[(per_seed.cascade == bc)
                              & (per_seed.prep == bp)].cum_inf.mean()
    attr["cascade_pct_of_all"] = 100 * attr.shap_cascade / baseline_total
    attr["prep_pct_of_all"] = 100 * attr.shap_prep / baseline_total
    attr["not_averted_pct"] = 100 * (baseline_total - attr.shap_cascade
                                     - attr.shap_prep) / baseline_total
    attr.to_csv(OUT / "attribution.csv", index=False)
    print(f"\n=== Share of ALL {baseline_total:,.0f} baseline infections ===")
    print(attr[attr.cascade != "S0_status_quo"]
          [["cascade", "prep", "cascade_pct_of_all", "prep_pct_of_all",
            "not_averted_pct"]].round(1).to_string(index=False))
    fig_attribution(attr, baseline_total)

    vir = fig_viremia(d)
    print("\n=== Adults 15+ with detectable viral load (%), no PrEP ===")
    print(vir.round(3).to_string(index=False))

    fig_prep_by_counterfactual(grid, per_seed)
    fig_surface(grid)
    fig_heatmap(grid)
    fig_xyz(grid)
    print("\nwrote outputs/{grid,per_seed,substitution,"
          "residual_decomposition}.csv")
    print("      figures/{surface,surface_heatmap,xyz_scatter}.png")


if __name__ == "__main__":
    main()
