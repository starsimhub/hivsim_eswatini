"""Build the exp 034 scenario workbook from scenario_table.csv.

Four sheets: Notes (provenance and caveats), Scenarios (verbose definitions),
Data (all rows, with the derived columns as live formulas), Summary (the 15+
view, for reading).

The derived columns -- infections averted, % averted, % incidence difference --
are written as VALUES, not as Excel formulas. That is a deliberate departure
from the usual preference for live formulas.

openpyxl writes formulas with no cached value, so they read as blank to
pandas and to most previewers until a spreadsheet engine evaluates them. The
recalc step that would do that needs LibreOffice, which is available on
neither this Windows laptop (its helper requires socket.AF_UNIX) nor the VM.
Unverified formulas in a file nobody can recalculate here are worse than
values: a bad reference would ship looking fine and surface only when a
co-author opened it.

The values come from the same pandas computation as scenario_table.csv, and
each derivation is stated on the Notes sheet so the arithmetic stays checkable.

Usage (repo root):
  python experiments/034_plausible_scenarios/make_workbook.py
"""

import sys
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
SRC = OUT / "scenario_table.csv"
DEST = OUT / "eswatini_scenario_table.xlsx"

FONT = "Arial"
HDR_FILL = PatternFill("solid", fgColor="1F3864")
HDR_FONT = Font(name=FONT, bold=True, color="FFFFFF", size=10)
BODY = Font(name=FONT, size=10)
BOLD = Font(name=FONT, size=10, bold=True)
TITLE = Font(name=FONT, size=13, bold=True)
NOTE = Font(name=FONT, size=9, italic=True, color="595959")
THIN = Side(style="thin", color="BFBFBF")
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

CASC_ORDER = ["status_quo", "testing_only", "unaids_95", "99_96_98"]
CASC_PRETTY = {"status_quo": "1. Status quo",
               "testing_only": "2. Testing x3 only",
               "unaids_95": "3. 95-95-95 every group",
               "99_96_98": "4. ART 96% + VLS 98%"}
PREP_ORDER = ["none", "fsw", "agyw_risk", "agyw_all", "women_25_34",
              "broad_60", "broad_90"]
PREP_PRETTY = {"none": "No PrEP", "fsw": "FSW 60%",
               "agyw_risk": "+ higher-risk AGYW",
               "agyw_all": "+ all AGYW",
               "women_25_34": "+ women 25-34 (broad)",
               "broad_60": "All women 15-34 at 60%",
               "broad_90": "All women 15-34 at 90%"}
AGE_ORDER = ["15-24", "25-49", "50+", "15-49", "15+ (all)"]
SEX_ORDER = ["Both sexes", "Women", "Men"]

# (source column, header, number format, width)
COLS = [
    ("cascade_pretty",       "ART cascade scenario",            None,       26),
    ("prep_pretty",          "LA-PrEP scenario",                None,       22),
    ("sex",                  "Sex",                             None,       9),
    ("age_group",            "Age group",                       None,       11),
    ("aware",                "Aware of status",                 "0.0%",     14),
    ("art_given_aware",      "On ART | aware",                  "0.0%",     14),
    ("vls_given_art",        "Suppressed | on ART",             "0.0%",     18),
    ("vls_of_plhiv",         "Suppressed | PLHIV",              "0.0%",     18),
    ("plhiv_2030",           "PLHIV (2030)",                    "#,##0",    14),
    ("_prep_cov",            "PrEP coverage, % of HIV-negative adults 15-49 "
                             "(2030)",                          "0.0%",     26),
    ("incidence_2026",       "Incidence 2026 (scenario start)", "0.000",    22),
    ("_pct_decline",         "% incidence decline, 2026-2040",  "0.0%",     22),
    ("baseline_infections",  "Infections, status quo",          "#,##0",    20),
    ("cum_infections",       "Infections, scenario",            "#,##0",    20),
    ("_averted",             "Infections averted",              "#,##0",    18),
    ("_pct_averted",         "% infections averted",            "0.0%",     18),
    ("baseline_incidence_2030", "Incidence 2030, status quo",   "0.000",    22),
    ("incidence_2030",       "Incidence 2030, scenario",        "0.000",    22),
    ("_pct_inc_2030",        "% lower incidence, 2030",         "0.0%",     20),
    ("baseline_incidence_2040", "Incidence 2040, status quo",   "0.000",    22),
    ("incidence_2040",       "Incidence 2040, scenario",        "0.000",    22),
    ("_pct_inc_2040",        "% lower incidence, 2040",         "0.0%",     20),
]

NOTES = [
    ("Eswatini HIV scenario table", TITLE),
    ("", None),
    ("Source: experiments/034_plausible_scenarios, model-v1.6 "
     "(starsim 3.5.2 / stisim 1.5.11).", BODY),
    ("Agent-based model HIVsim, fitted to Eswatini prevalence by age, sex and "
     "year (PHIA/SHIMS), AIDS mortality, adult incidence, and the SHIMS3 2021 "
     "treatment cascade.", BODY),
    ("", None),
    ("WHAT EACH NUMBER IS", BOLD),
    ("Cascade percentages and PLHIV are read at 2030, when every scenario is "
     "fully scaled up.", BODY),
    ("Infections are cumulative over 2026-2040.", BODY),
    ("Incidence is new infections per 100 susceptible person-years.", BODY),
    ("Every comparison is against the status-quo, no-PrEP cell OF THE SAME "
     "stratum -- not against the overall baseline.", BODY),
    ("Incidence 2026 is the scenario start year and the reference for the "
     "2026-2040 decline. Scenarios ramp THROUGH 2026 rather than switching on "
     "at its start, so arms can differ slightly at that point -- that is the "
     "intervention acting in its first year, not a contaminated baseline.",
     BODY),
    ("% incidence decline 2026-2040 is measured WITHIN a scenario. The "
     "'% lower incidence' columns are measured AGAINST the status quo at the "
     "stated year. The two answer different questions.", BODY),
    ("Strata are Women, Men and Both sexes, by 15-24 / 25-49 / 50+ / 15-49 / "
     "15+. Cascade percentages CANNOT be summed across strata by hand -- they "
     "are ratios with different denominators -- so the combined rows are "
     "computed from pooled counts, not averaged.", BODY),
    ("PrEP coverage is expressed against HIV-NEGATIVE adults 15-49, since "
     "only they are eligible. It is a SCENARIO-level figure repeated across "
     "every stratum row: the model exports PrEP recipients as a single total "
     "with no age or sex breakdown, so a per-stratum coverage cannot be "
     "derived. A few sex workers above 49 sit outside the denominator, which "
     "inflates it very slightly.", BODY),
    ("Two summary sheets, both sexes combined: 15-49 (the conventional HIV "
     "reporting denominator) and all adults 15+ (the denominator the cascade "
     "and viremia findings are framed in). They give different answers -- 50+ "
     "carries a large and growing share of unsuppressed HIV -- so check which "
     "one a number came from before quoting it.", BODY),
    ("", None),
    ("UNCERTAINTY", BOLD),
    ("A single calibrated parameter set was used, so all variation is "
     "stochastic (between-seed) and NOT parameter uncertainty. Nothing here is "
     "a credible interval.", BODY),
    ("", None),
    ("KNOWN LIMITATIONS -- all three make the cascade scenarios OPTIMISTIC", BOLD),
    ("1. The model has no never-testing subgroup. SHIMS3 implies ~3% of adults "
     "over 50 remain unaware after decades; the model reaches ~1% and "
     "structurally cannot do worse.", BODY),
    ("2. stisim's stratified ART pathway ignores the linkage delay and the "
     "initiation probability, so linkage can be driven higher than any real "
     "programme achieves.", BODY),
    ("3. Suppression given ART carries no age gradient, reading 6-10 "
     "percentage points too high in 15-24s against SHIMS3. The youth "
     "suppression deficit a real adherence or long-acting-ART programme would "
     "target is therefore absent.", BODY),
    ("", None),
    ("TREAT WITH CAUTION", BOLD),
    ("Testing x3 only produces MORE infections in men and MORE AIDS deaths "
     "overall than status quo. This is a real result, not noise. See the "
     "Scenarios sheet.", BODY),
    ("", None),
    ("HOW THE DERIVED COLUMNS ARE CALCULATED", BOLD),
    ("Infections averted = (infections, status quo) - (infections, scenario), "
     "for the same sex and age group.", BODY),
    ("% infections averted = infections averted / (infections, status quo).",
     BODY),
    ("% lower incidence = ((incidence, status quo) - (incidence, scenario)) / "
     "(incidence, status quo), at the stated year.", BODY),
    ("These are stored as values rather than live formulas: openpyxl writes "
     "formulas with no cached result, and neither machine used here has a "
     "spreadsheet engine available to evaluate and verify them. Shipping "
     "unverified formulas would hide a bad cell reference until someone "
     "opened the file.", NOTE),
]


def build():
    if not SRC.exists():
        sys.exit(f"missing {SRC} -- run scenario_table.py first")
    t = pd.read_csv(SRC)
    t["cascade_pretty"] = t.cascade_name.map(CASC_PRETTY)
    t["prep_pretty"] = t.prep_name.map(PREP_PRETTY)
    t["_c"] = t.cascade_name.map({c: i for i, c in enumerate(CASC_ORDER)})
    t["_p"] = t.prep_name.map({p: i for i, p in enumerate(PREP_ORDER)})
    t["_a"] = t.age_group.map({a: i for i, a in enumerate(AGE_ORDER)})
    t["_s"] = t.sex.map({s: i for i, s in enumerate(SEX_ORDER)})
    t = t.sort_values(["_c", "_p", "_a", "_s"]).reset_index(drop=True)

    wb = Workbook()

    # --- Notes -------------------------------------------------------------
    ws = wb.active
    ws.title = "Notes"
    ws.column_dimensions["A"].width = 110
    for i, (text, font) in enumerate(NOTES, start=1):
        c = ws.cell(row=i, column=1, value=text)
        c.font = font or BODY
        c.alignment = Alignment(wrap_text=True, vertical="top")

    # --- Scenarios ---------------------------------------------------------
    ws = wb.create_sheet("Scenarios")
    hdr = ["ART cascade scenario", "LA-PrEP scenario", "Full description"]
    for j, h in enumerate(hdr, start=1):
        c = ws.cell(row=1, column=j, value=h)
        c.font, c.fill, c.border = HDR_FONT, HDR_FILL, BOX
        c.alignment = Alignment(vertical="center")
    for w, col in zip((26, 22, 140), "ABC"):
        ws.column_dimensions[col].width = w
    seen, r = set(), 2
    for row in t.itertuples():
        key = (row.cascade_pretty, row.prep_pretty)
        if key in seen:
            continue
        seen.add(key)
        for j, v in enumerate([row.cascade_pretty, row.prep_pretty,
                               row.description], start=1):
            c = ws.cell(row=r, column=j, value=v)
            c.font, c.border = BODY, BOX
            c.alignment = Alignment(wrap_text=True, vertical="top")
        r += 1
    r += 1
    ws.cell(row=r, column=1, value="Caveat on scenario 2 (testing only)").font = BOLD
    ws.cell(row=r + 1, column=1,
            value="In men this scenario produces MORE infections than status "
                  "quo. The ART coverage target is a fraction of all PLHIV and "
                  "fills from the diagnosed pool, so tripling testing enlarges "
                  "that pool without enlarging the number treated. Direction "
                  "was consistent across seeds; magnitude should be read "
                  "against the seed spread before being reported."
            ).alignment = Alignment(wrap_text=True, vertical="top")
    ws.cell(row=r + 1, column=1).font = NOTE

    # --- Data --------------------------------------------------------------
    ws = wb.create_sheet("Data")
    for j, (_, head, fmt, width) in enumerate(COLS, start=1):
        c = ws.cell(row=1, column=j, value=head)
        c.font, c.fill, c.border = HDR_FONT, HDR_FILL, BOX
        c.alignment = Alignment(wrap_text=True, vertical="center",
                                horizontal="center")
        ws.column_dimensions[get_column_letter(j)].width = width
    ws.row_dimensions[1].height = 34

    for i, row in enumerate(t.itertuples(), start=2):
        # Percentages are stored as FRACTIONS so the 0.0% format renders them
        # correctly; the CSV keeps them in percentage points.
        f = {
            "_averted": row.infections_averted,
            "_pct_averted": row.pct_infections_averted / 100.0,
            "_pct_decline": row.pct_decline_2026_2040 / 100.0,
            "_prep_cov": row.prep_coverage_hivneg_15_49 / 100.0,
            "_pct_inc_2030": row.pct_diff_incidence_2030 / 100.0,
            "_pct_inc_2040": row.pct_diff_incidence_2040 / 100.0,
        }
        for j, (src, _h, fmt, _w) in enumerate(COLS, start=1):
            v = f[src] if src in f else getattr(row, src)
            if isinstance(v, float) and pd.isna(v):
                v = None
            c = ws.cell(row=i, column=j, value=v)
            c.font, c.border = BODY, BOX
            if fmt:
                c.number_format = fmt
            if src in ("cascade_pretty", "prep_pretty"):
                c.alignment = Alignment(vertical="center")
    ws.freeze_panes = "E2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(COLS))}{len(t) + 1}"

    # --- Summary sheets ----------------------------------------------------
    # BOTH reporting strata get their own tab, rather than one standing in for
    # the other: 15-49 both sexes is the conventional HIV reporting
    # denominator, while 15+ is all adults and is the denominator the cascade
    # and viremia findings are framed in. They give different answers and the
    # abstract uses both.
    keep = [c for c in COLS if c[0] not in
            ("sex", "age_group", "plhiv_2030", "baseline_incidence_2030",
             "incidence_2030", "_pct_inc_2030")]
    summaries = [("Summary 15-49", "15-49"),
                 ("Summary all adults 15+", "15+ (all)")]
    counts = []
    for sheet_name, age in summaries:
        ws = wb.create_sheet(sheet_name)
        for j, (_, head, fmt, width) in enumerate(keep, start=1):
            c = ws.cell(row=1, column=j, value=head)
            c.font, c.fill, c.border = HDR_FONT, HDR_FILL, BOX
            c.alignment = Alignment(wrap_text=True, vertical="center",
                                    horizontal="center")
            ws.column_dimensions[get_column_letter(j)].width = width
        ws.row_dimensions[1].height = 34
        sub = t[(t.age_group == age)
                & (t.sex == "Both sexes")].reset_index(drop=True)
        for i, row in enumerate(sub.itertuples(), start=2):
            f = {
                "_averted": row.infections_averted,
                "_pct_averted": row.pct_infections_averted / 100.0,
                "_pct_decline": row.pct_decline_2026_2040 / 100.0,
                "_prep_cov": row.prep_coverage_hivneg_15_49 / 100.0,
                "_pct_inc_2040": row.pct_diff_incidence_2040 / 100.0,
            }
            for j, (src, _h, fmt, _w) in enumerate(keep, start=1):
                v = f[src] if src in f else getattr(row, src)
                if isinstance(v, float) and pd.isna(v):
                    v = None
                c = ws.cell(row=i, column=j, value=v)
                c.font, c.border = BODY, BOX
                if fmt:
                    c.number_format = fmt
        ws.freeze_panes = "C2"
        ws.auto_filter.ref = f"A1:{get_column_letter(len(keep))}{len(sub) + 1}"
        counts.append(f"{sheet_name}: {len(sub)} rows")

    wb.save(DEST)
    print(f"wrote {DEST}")
    print(f"  Data: {len(t)} rows | " + " | ".join(counts))


if __name__ == "__main__":
    build()
