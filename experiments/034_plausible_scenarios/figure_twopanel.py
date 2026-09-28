"""Exp 034 -- the two-panel abstract figure.

Panel A: infections averted, grouped by PrEP strategy, with the two ART
         cascade scenarios side by side inside each group and the burden
         split into what each lever prevented (Shapley) and what neither did.
Panel B: HIV incidence 2020-2040 for the baseline, the age-gaps-filled
         cascade, and that cascade with each PrEP strategy added.

Panel A is grouped by PrEP rather than by cascade -- Adam's reorganisation --
so the question it answers is "for this PrEP programme, what does filling the
cascade gaps add?", read within a group, rather than across the figure.

Run where the per-run parquet live (the VM): panel B needs the full incidence
time series, which the summary tables do not carry.

Usage:  python experiments/034_plausible_scenarios/figure_twopanel.py
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
OUT, FIG, SIMS = HERE / "outputs", HERE / "figures", HERE / "outputs" / "sims"
FIG.mkdir(parents=True, exist_ok=True)

INK, MUTED = "#222222", "#6b6b6b"
AGE, SEX = "15+ (all)", "Both sexes"

CASCADES = [("status_quo", "S0_status_quo", "Status quo ART (95-95-95)"),
            ("unaids_95", "S2_unaids_95", "Age-gaps filled ART")]
PREPS = [("none", "P0_none", "No PrEP"),
         ("fsw", "P1_fsw", "FSW 60%"),
         ("agyw_risk", "P2_agyw_risk", "+ higher-risk AGYW"),
         ("agyw_all", "P3_agyw_all", "+ all AGYW"),
         ("women_25_34", "P4_women_25_34", "+ women 25-34")]

CASC_COL = {"status_quo": "#7fa9d0", "unaids_95": "#2c6fbb"}
PREP_COL = {"none": "#f0b862", "fsw": "#e8a33d", "agyw_risk": "#d98f27",
            "agyw_all": "#c47c18", "women_25_34": "#a86710"}
GREY = "#d5dade"


def adult_cols(d, stem):
    out = []
    for c in d.columns:
        if not c.startswith(f"popagesex.{stem}_"):
            continue
        p = c.split("_")
        if len(p) >= 3 and p[-3] in ("f", "m") and int(p[-2]) >= 15:
            out.append(c)
    return out


def incidence(casc_key, prep_key):
    """Per-year incidence, 15+, mean and SD across seeds."""
    f = glob.glob(str(SIMS / f"{casc_key}__{prep_key}__*.parquet"))
    if not f:
        return None
    d = pd.concat([pd.read_parquet(x) for x in f], ignore_index=True)
    ni = adult_cols(d, "new_infections")
    inf = adult_cols(d, "n_infected")
    alv = adult_cols(d, "n_alive")
    d["inc"] = 100 * d[ni].sum(axis=1) / (d[alv].sum(axis=1) - d[inf].sum(axis=1))
    g = d.groupby("timevec").inc
    return g.mean(), g.std(ddof=1).fillna(0)


def main():
    csv = OUT / "scenario_table.csv"
    if not csv.exists():
        sys.exit("need outputs/scenario_table.csv -- run scenario_table.py")
    t = pd.read_csv(csv)
    t = t[(t.age_group == AGE) & (t.sex == SEX)]

    def cum(c, p):
        r = t[(t.cascade_name == c) & (t.prep_name == p)]
        return float(r.cum_infections.iloc[0]) if len(r) else np.nan

    base = cum("status_quo", "none")

    fig, axes = plt.subplots(1, 2, figsize=(19, 7.2),
                             gridspec_kw=dict(width_ratios=[1.15, 1]))

    # --- Panel A ----------------------------------------------------------
    ax = axes[0]
    width, gap = 0.34, 0.03
    for i, (pk, _, plab) in enumerate(PREPS):
        for j, (ck, _, clab) in enumerate(CASCADES):
            x = i + (j - 0.5) * (width + gap)
            A = base - cum(ck, "none")          # cascade alone
            B = base - cum("status_quo", pk)    # PrEP alone
            J = base - cum(ck, pk)              # both
            I = J - A - B
            cv, pv = A + I / 2, B + I / 2
            rem = base - cv - pv
            first = (i == 0 and j == 0)
            ax.bar(x, cv, width, color=CASC_COL["unaids_95"], zorder=2,
                   label="averted by the ART cascade" if first else None)
            ax.bar(x, pv, width, bottom=cv, color=PREP_COL["fsw"], zorder=2,
                   label="averted by PrEP" if first else None)
            ax.bar(x, rem, width, bottom=cv + pv, color=GREY, zorder=2,
                   label="not averted" if first else None)
            for val, bot, light in ((cv, 0, True), (pv, cv, True),
                                    (rem, cv + pv, False)):
                if val > base * 0.05:
                    ax.text(x, bot + val / 2, f"{100*val/base:.0f}%",
                            ha="center", va="center", fontsize=7.6, zorder=3,
                            color="white" if light else INK)
            # Bar identity is the CASCADE here, rotated because the two labels
            # are long and the bars are narrow.
            ax.text(x, -base * 0.015, clab, ha="center", va="top",
                    rotation=90, fontsize=7.2, color=MUTED)
    ax.axhline(base, color=INK, lw=1.3, ls="--", zorder=4)
    ax.text(len(PREPS) - 0.55, base * 1.015,
            f"baseline: {base:,.0f} infections", ha="right", fontsize=8,
            color=INK)
    ax.set_xticks(range(len(PREPS)))
    ax.set_xticklabels([p[2] for p in PREPS], fontsize=9)
    ax.tick_params(axis="x", pad=104)
    ax.set_ylabel("cumulative HIV infections in adults 15+, 2026-2040")
    ax.set_ylim(0, base * 1.06)
    ax.grid(axis="y", alpha=0.28, zorder=0)
    ax.legend(fontsize=8.5, frameon=False, loc="lower left")
    ax.set_title("A. Who prevents the infections", fontsize=11.5, loc="left")

    # --- Panel B ----------------------------------------------------------
    ax = axes[1]
    curves = [("S0_status_quo", "P0_none", "Status quo ART (95-95-95)",
               "#2c3e50", "-")]
    curves += [("S2_unaids_95", "P0_none", "Age-gaps filled ART",
                CASC_COL["unaids_95"], "-")]
    curves += [("S2_unaids_95", pk_s, f"Age-gaps filled ART + {lab}",
                PREP_COL[pk], "-")
               for pk, pk_s, lab in PREPS if pk != "none"]
    for casc, prep, lab, col, ls in curves:
        r = incidence(casc, prep)
        if r is None:
            print(f"  ! missing {casc} {prep}")
            continue
        m, sd = r
        w = (m.index >= 2020) & (m.index <= 2040)
        ax.plot(m.index[w], m.values[w], ls, lw=2.2, color=col, label=lab,
                zorder=3)
        ax.fill_between(m.index[w], (m - sd).values[w], (m + sd).values[w],
                        color=col, alpha=0.10, lw=0, zorder=2)
    ax.axvline(2026, ls=":", color=MUTED, lw=1)
    ax.set_xlim(2020, 2040)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("year")
    ax.set_ylabel("HIV incidence, adults 15+ (per 100 person-years)")
    ax.grid(alpha=0.28)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("B. Incidence over time", fontsize=11.5, loc="left")
    ax.text(2026.2, ax.get_ylim()[1] * 0.02, "scenarios begin", fontsize=7.4,
            color=MUTED)

    fig.suptitle("Long-acting PrEP and ART cascade gaps in Eswatini, "
                 "2026-2040", fontsize=13)
    fig.text(0.5, -0.03,
             "Panel A: each bar is the full status-quo burden; percentages are "
             "shares of ALL baseline infections, so grey is what neither lever "
             "prevents. The two levers overlap, so the split uses Shapley "
             "attribution,\nwhich shares that overlap evenly. Panel B: bands "
             "are +/-1 SD across 50 seeds. 'Status quo ART (95-95-95)' is "
             "Eswatini's current national attainment; 'age-gaps filled' "
             "reaches 95-95-95 within every age and sex group.",
             ha="center", fontsize=8.2, color=MUTED)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    dest = FIG / "abstract_two_panel.png"
    fig.savefig(dest, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
