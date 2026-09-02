"""Fix the extrapolated 45+ band in data/art_coverage.csv from SHIMS3.

Run `python art_coverage_construction.py` to regenerate.

What was wrong
--------------
`data/art_coverage.csv` has four age bands -- [15,25), [25,35), [35,45),
[45,100) -- but its source, `external_data/SWAZILAND_calibration_nationalARTprevalence.csv`,
has only three, stopping at [35,45). The 45+ band was filled by **copying the
[35,45) value**, in every year and for both sexes. Exp 020 flagged this when the
derived suppression-given-ART ratio exceeded 1.0 above age 45 -- more suppressed
people than treated ones, which is impossible -- and that was the tell.

What this changes, and only this
--------------------------------
SHIMS3 2021 Table 9.1.A reports ART coverage among all PLHIV for a **50+** band,
so the 2021 [45,100) cells can be replaced with a measured value instead of a
copy:

    women 2021 [45,100)   0.970 (copied)  ->  0.966 (SHIMS3 50+)
    men   2021 [45,100)   0.902 (copied)  ->  0.953 (SHIMS3 50+)

The men's cell moves 5.1 points; the women's barely moves, which is why the
problem went unnoticed -- female ART coverage really is nearly flat above 35.

**2011 and 2016 keep the copied value**, because neither SHIMS1 nor SHIMS2
publishes a 50+ ART-coverage band (SHIMS2 Table 10.3.B stops at 35-49 and then
jumps to 15+). Those cells remain an extrapolation and are marked as such in
PROVENANCE below rather than silently corrected. [45,100) also spans 45-49,
which belongs to SHIMS3's 35-49 band, so even the 2021 value is an
approximation -- it is weighted toward 50+ because most of the band is.

Also fixed: 2016 men 15-24
--------------------------
The file carried 0.360 there. Traced to SHIMS2 **Table 8.3.A**, the 20-24 row
(self-reported on ART = 34.6%, n = 29, parenthesised) -- see SHIMS2_MEN_15_24
below for the full argument. Corrected to 0.556, the value on the biomarker
basis the file's other cells use.

The Akullian et al. 2020 appendix was checked as a possible source and does not
supply age/sex ART coverage as data at all: that model *fitted* coverage through
`ART_Link_Max/Mid/Rate`, comparing emergent coverage against a 15-49 aggregate.

Idempotent: SRC and OUT are the same path, so re-running is safe -- every fix
assigns an absolute value rather than adjusting the existing one.

Sex convention
--------------
stisim's, which `parse_coverage` documents as **0 = female, 1 = male**. The
existing file follows it correctly -- all nine 2011/2016 age x sex cells match
their raw source once the source's opposite convention is applied.
"""

import pandas as pd

SRC = "data/art_coverage.csv"
OUT = "data/art_coverage.csv"

# SHIMS3 2021 Table 9.1.A, "On Treatment", percentage on ART among all PLHIV,
# 50+ band. Sex keyed by stisim's convention: 0 = female, 1 = male.
SHIMS3_50PLUS = {0: 0.966, 1: 0.953}
FIX_YEAR, FIX_BAND = 2021, "[45,100)"

# 2016 men 15-24: the file carried 0.36, which traces to SHIMS2 Table 8.3.A's
# **20-24** row (self-reported on ART = 34.6%, n = 29, parenthesised by the
# report as a small-denominator estimate). Applied to the whole [15,25) band it
# is wrong twice: the published Total 15-24 self-report figure is 46.8% -- men
# 15-19 are at 59.5%, well above 20-24, and pull the band up -- and every other
# cell in this file is on the *biomarker* basis, not self-report.
#
# Deriving each 2016 cell from SHIMS2 Table 10.3.B as
# diagnosed x on-ART-among-diagnosed gives:
#
#     men   15-24  0.556    file 0.36   <-- off by 20 points
#     men   25-34  0.583    file 0.60
#     men   35-49  0.762    file 0.76   <-- exact match: this is the file's basis
#     women 15-24  0.632    file 0.60
#     women 25-34  0.779    file 0.80
#     women 35-49  0.855    file 0.90
#
# Men 35-49 reproducing 0.762 to three decimals identifies the basis, and every
# cell sits within ~5 points of its biomarker value except men 15-24. So 0.36 is
# an anomaly, not a choice, and 0.556 is the value on the file's own basis.
#
# This matters out of proportion to one cell: young men are where the prevalence
# fit is worst (-65% vs PHIA), where 020's rare-event floor binds, and
# understating their ART coverage by 20 points leaves young HIV-positive men
# transmitting far more than they should -- inflating onward transmission to
# young women, whose prevalence the model also badly misses.
SHIMS2_MEN_15_24 = 0.556

PROVENANCE = {
    (2004, "*"): "zero by construction: ART programme had not started",
    (2011, "[15,25)"): "SWAZILAND_calibration_nationalARTprevalence.csv (sex-flipped to stisim convention)",
    (2011, "[25,35)"): "as above",
    (2011, "[35,45)"): "as above",
    (2011, "[45,100)"): "EXTRAPOLATED: copy of [35,45); no 50+ band published for 2011",
    (2016, "[15,25)"): "women: source file. men: 0.556, derived from SHIMS2 Table 10.3.B -- FIXED by this script, was 0.360 (Table 8.3.A 20-24 row)",
    (2016, "[25,35)"): "as above",
    (2016, "[35,45)"): "as above",
    (2016, "[45,100)"): "EXTRAPOLATED: copy of [35,45); no 50+ band published for 2016",
    (2021, "[15,25)"): "as above -- NB men 15-24 = 0.833 vs SHIMS3 Table 9.1.A 0.877",
    (2021, "[25,35)"): "as above; matches SHIMS3 Table 9.1.A to <0.6 pp",
    (2021, "[35,45)"): "as above; matches SHIMS3 Table 9.1.A to <0.2 pp",
    (2021, "[45,100)"): "SHIMS3 2021 Table 9.1.A, 50+ band -- FIXED by this script",
}


def build():
    df = pd.read_csv(SRC)
    before = df.copy()
    mask = (df.Year == FIX_YEAR) & (df.AgeBin == FIX_BAND)
    assert mask.sum() == 2, f"expected 2 cells to fix, found {mask.sum()}"
    df.loc[mask, "p_art"] = df.loc[mask, "Gender"].map(SHIMS3_50PLUS)

    # 2016 men 15-24 (stisim convention: Gender 1 = male)
    m2 = (df.Year == 2016) & (df.AgeBin == "[15,25)") & (df.Gender == 1)
    assert m2.sum() == 1, f"expected 1 cell, found {m2.sum()}"
    df.loc[m2, "p_art"] = SHIMS2_MEN_15_24

    assert df.p_art.between(0, 1).all(), "p_art must be a proportion"

    changed = df.loc[mask].merge(
        before.loc[mask, ["Year", "Gender", "AgeBin", "p_art"]],
        on=["Year", "Gender", "AgeBin"], suffixes=("_new", "_old"))
    return df, changed


if __name__ == "__main__":
    df, changed = build()
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT}\n")
    changed["sex"] = changed.Gender.map({0: "female", 1: "male"})
    print(changed[["Year", "sex", "AgeBin", "p_art_old", "p_art_new"]]
          .to_string(index=False))
    print("\nBands still extrapolated (no 50+ published for those rounds):")
    for (yr, band), note in PROVENANCE.items():
        if note.startswith("EXTRAPOLATED"):
            print(f"  {yr} {band}")
