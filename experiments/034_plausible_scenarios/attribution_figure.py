"""Rebuild the cascade/PrEP attribution figure from the summary table.

Reads the scenario table rather than the raw per-run parquet, so it works
without the VM. Everything the Shapley split needs is already in there:

    A = cascade alone   (that cascade, no PrEP,  vs status quo)
    B = PrEP alone      (status quo, that PrEP,  vs status quo)
    J = both together   (that cascade + PrEP,    vs status quo)

The levers overlap, so J < A + B. The overlap I = J - A - B is shared evenly:
cascade gets A + I/2, PrEP gets B + I/2, and the two sum exactly to J. Giving
the whole overlap to whichever lever is named first would shift the split by
more than ten percentage points, which is why it is not done that way.

Bars run the full status-quo burden, so the grey block is the share neither
lever prevents -- the quantity the analysis is actually about. Percentages are
of ALL baseline infections, not of those averted.

Usage (repo root):
  python experiments/034_plausible_scenarios/attribution_figure.py
  python experiments/034_plausible_scenarios/attribution_figure.py --age "15-49"
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT, FIG = HERE / "outputs", HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)

INK, MUTED = "#222222", "#6b6b6b"
CASC_LAB = {"status_quo": "status quo ART",
            "testing_only": "testing x3 only",
            "unaids_95": "95-95-95\nevery group",
            "99_96_98": "ART 96% + VLS 98%"}
PREP_LAB = {"none": "no PrEP", "fsw": "FSW\n60%",
            "agyw_risk": "+ higher-\nrisk AGYW", "agyw_all": "+ all\nAGYW",
            "women_25_34": "+ women\n25-34"}

# Default is the abstract's set: two ART cascade scenarios against the four
# LA-PrEP strategies, plus the no-PrEP reference. The status-quo/no-PrEP bar
# is the baseline itself, so it renders as 100% not averted -- kept because it
# anchors the bar height rather than leaving the reader to infer it.
DEF_CASC = ["status_quo", "unaids_95"]
DEF_PREP = ["none", "fsw", "agyw_risk", "agyw_all", "women_25_34"]


def load(age, sex):
    """Prefer the CSV; fall back to the workbook if it is locked or stale."""
    # Check the CSV actually CONTAINS the requested stratum, not merely that
    # it has the right column names -- an older table has prep_name and
    # age_group but no both-sexes or 15-49 rows, and would silently return
    # nothing.
    csv = OUT / "scenario_table.csv"
    if csv.exists():
        try:
            t = pd.read_csv(csv)
            if "prep_name" in t.columns:
                sub = t[(t.age_group == age) & (t.sex == sex)]
                if len(sub):
                    return sub
                print(f"  ! {csv.name} has no {sex} / {age} rows "
                      f"(stale) -- falling back to the workbook")
        except (PermissionError, OSError):
            print(f"  ! {csv.name} is locked -- falling back to the workbook")
    xl = next((p for p in (OUT / "eswatini_scenario_table_v2.xlsx",
                           OUT / "eswatini_scenario_table.xlsx")
               if p.exists()), None)
    if xl is None:
        sys.exit("no scenario_table.csv and no workbook in outputs/")
    d = pd.read_excel(xl, sheet_name="Data")
    ren = {"ART cascade scenario": "cascade_pretty",
           "LA-PrEP scenario": "prep_pretty", "Sex": "sex",
           "Age group": "age_group", "Infections, status quo":
           "baseline_infections", "Infections, scenario": "cum_infections"}
    d = d.rename(columns=ren)
    cmap = {"1. Status quo": "status_quo", "2. Testing x3 only": "testing_only",
            "3. 95-95-95 every group": "unaids_95",
            "4. ART 96% + VLS 98%": "99_96_98"}
    pmap = {"No PrEP": "none", "FSW 60%": "fsw",
            "+ higher-risk AGYW": "agyw_risk", "+ all AGYW": "agyw_all",
            "+ women 25-34 (broad)": "women_25_34"}
    d["cascade_name"] = d.cascade_pretty.map(cmap)
    d["prep_name"] = d.prep_pretty.map(pmap)
    d["aware"] = d["Aware of status"]
    d["art_given_aware"] = d["On ART | aware"]
    d["vls_given_art"] = d["Suppressed | on ART"]
    return d[(d.age_group == age) & (d.sex == sex)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--age", default="15+ (all)")
    ap.add_argument("--sex", default="Both sexes")
    ap.add_argument("--cascades", nargs="*", default=DEF_CASC)
    ap.add_argument("--preps", nargs="*", default=DEF_PREP)
    ap.add_argument("--out", default="attribution.png")
    a = ap.parse_args()
    t = load(a.age, a.sex)
    if t.empty:
        sys.exit(f"no rows for age={a.age!r} sex={a.sex!r}")

    def cum(c, p):
        r = t[(t.cascade_name == c) & (t.prep_name == p)]
        return float(r.cum_infections.iloc[0]) if len(r) else np.nan

    base = cum("status_quo", "none")

    def triplet(c):
        r = t[(t.cascade_name == c) & (t.prep_name == "none")]
        if not len(r):
            return ""
        r = r.iloc[0]
        return (f"{100*r.aware:.0f}-{100*r.art_given_aware:.0f}"
                f"-{100*r.vls_given_art:.0f}")

    cascs = [c for c in a.cascades if not np.isnan(cum(c, "none"))]
    preps = [p for p in a.preps if not np.isnan(cum(cascs[0], p))]
    n = len(preps)
    width = min(0.8 / n, 0.26)
    gap = 0.01
    fig, ax = plt.subplots(figsize=(4.6 * len(cascs) + 4, 6.8))
    for j, p in enumerate(preps):
        xs = [i + (j - (n - 1) / 2) * (width + gap) for i in range(len(cascs))]
        for x, c in zip(xs, cascs):
            A = base - cum(c, "none")               # cascade alone
            B = base - cum("status_quo", p)         # PrEP alone
            J = base - cum(c, p)                    # both
            I = J - A - B
            cv, pv = A + I / 2, B + I / 2
            rem = base - cv - pv
            first = (j == 0 and x == xs[0])
            ax.bar(x, cv, width, color="#2c6fbb", zorder=2,
                   label="averted by the cascade" if first else None)
            ax.bar(x, pv, width, bottom=cv, color="#e8a33d", zorder=2,
                   label="averted by PrEP" if first else None)
            ax.bar(x, rem, width, bottom=cv + pv, color="#d5dade", zorder=2,
                   label="not averted" if first else None)
            for val, bot, light in ((cv, 0, True), (pv, cv, True),
                                    (rem, cv + pv, False)):
                if val > base * 0.045:
                    ax.text(x, bot + val / 2, f"{100*val/base:.0f}%",
                            ha="center", va="center", fontsize=7.4, zorder=3,
                            color="white" if light else INK)
        for x in xs:
            ax.text(x, -base * 0.03, PREP_LAB[p], ha="center", va="top",
                    fontsize=6.4, color=MUTED)

    ax.axhline(base, color=INK, lw=1.4, ls="--", zorder=4)
    ax.text(len(cascs) - 0.42, base * 1.012,
            f"baseline: {base:,.0f} infections with no cascade improvement "
            f"and no PrEP", ha="right", fontsize=8, color=INK)
    ax.set_xticks(range(len(cascs)))
    ax.set_xticklabels([f"{CASC_LAB[c]}\n{triplet(c)}" for c in cascs],
                       fontsize=9)
    ax.tick_params(axis="x", pad=34)
    ax.set_ylabel("cumulative HIV infections, 2026-2040")
    ax.set_ylim(0, base * 1.08)
    ax.grid(axis="y", alpha=0.28, zorder=0)
    ax.legend(fontsize=9, frameon=False, ncol=3, loc="lower center",
              bbox_to_anchor=(0.5, 1.02))
    ax.set_title("Of every infection Eswatini would otherwise see, how many "
                 f"does each lever prevent?  ({a.sex}, {a.age})",
                 fontsize=12.5, pad=34)
    fig.text(0.5, -0.06,
             "Each bar is the full status-quo burden. Percentages are shares "
             "of ALL baseline infections, not of those averted, so the grey "
             "block is what neither lever prevents.\nThe two levers overlap, "
             "so the split uses Shapley attribution, which shares that "
             "overlap evenly and does not depend on which is counted first.",
             ha="center", fontsize=8.2, color=MUTED)
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    dest = FIG / a.out
    fig.savefig(dest, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest}  (baseline {base:,.0f}, {a.sex}, {a.age})")
    for c in cascs:
        for p in [x for x in preps if x != "none"]:
            A, B = base - cum(c, "none"), base - cum("status_quo", p)
            J = base - cum(c, p)
            I = J - A - B
            print(f"  {c:14} +{p:12} cascade {100*(A+I/2)/base:4.1f}%  "
                  f"PrEP {100*(B+I/2)/base:4.1f}%  "
                  f"neither {100*(base-A-B-I)/base:4.1f}%")


if __name__ == "__main__":
    main()
