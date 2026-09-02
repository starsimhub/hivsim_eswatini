"""Recode external_data/Swaziland_PPDV.csv onto one consistent sex convention.

Run `python ppdv_construction.py` to regenerate data/eswatini_ppdv.csv.

The problem
-----------
The raw file mixes two sex conventions *within the same row*:

  Prev, Prev_lb/ub, VLS, VLS_lb/ub, ppdv, ppdv_lb/ub   Gender 0 = FEMALE
  Inc, Inc_lb, Inc_ub                                  Gender 0 = MALE

Verified against primary sources, all three survey rounds:

  Prev 2016   G0 0.343 / G1 0.189   SHIMS2: women 34.3, men 18.9
  VLS  2016   G0 0.748 / G1 0.623   SHIMS2 Table 9.4.A: women 74.8, men 62.3
  VLS  2021   G0 0.886 / G1 0.824   SHIMS3 Table 8.2:   women 88.6, men 82.4
  Inc  2021   G0 0.20  / G1 1.45    SHIMS3 Table 5.1:   MEN 0.20, women 1.45
  Inc  2011   G0 1.65  / G1 3.14    SHIMS1:             MEN ~1.7, women ~3.1

and `ppdv` reproduces Prev x (1 - VLS) on the Prev/VLS convention
(2021: 0.316 x 0.114 = 0.0360, matching the file's 0.036024), which pins it to
that block rather than to Inc.

The fix
-------
Emit an explicit **`sex` column of 'm'/'f' strings**, so the ambiguity cannot
recur -- a string label cannot be silently inverted the way an integer can.
`Gender` is retained alongside it in the convention the other *target* files
use, **0 = Male, 1 = Female**, matching
`calibration_data/prevalence_by_age_sex.csv` and
`calibration_data/incidence_2021_VALIDATION_ONLY.csv`.

Note that is the opposite of stisim's internal convention (0 = female), which
`data/art_coverage.csv` correctly follows. The two conventions are irreducible:
model inputs follow stisim, targets follow the survey files. What matters is
that each file states which it uses -- see exp 021's config.yaml for the map.

Why the raw file is left untouched
----------------------------------
`external_data/` holds raw source data. Editing it in place would destroy the
provenance trail and make a future re-download indistinguishable from the
corrected version. This follows the existing pattern in this repo:
`mortality_construction.py` and `vls_construction.py` both read raw inputs and
emit processed files into `data/`.

ART_prev is carried through unmapped and flagged. It exists only for 2011, where
the two values are 0.34 and 0.33 -- too close to tell the sexes apart, and
matching neither the 0.32/0.23 in the external ART prevalence file. Its
convention is indeterminate and its meaning unclear, so it is not to be used
without checking.
"""

import numpy as np
import pandas as pd

RAW = "external_data/Swaziland_PPDV.csv"
OUT = "data/eswatini_ppdv.csv"

# Column blocks and the convention each uses in the RAW file.
FEMALE_IS_0 = ["Prev", "Prev_lb", "Prev_ub", "VLS", "VLS_lb", "VLS_ub",
               "ppdv", "ppdv_lb", "ppdv_ub"]
MALE_IS_0 = ["Inc", "Inc_lb", "Inc_ub"]
INDETERMINATE = ["ART_prev", "ART_prev_lb", "ART_prev_ub"]

# Primary-source anchors. The build fails loudly if the raw file ever changes
# such that these no longer hold, rather than silently producing a wrong recode.
CHECKS = [
    # (year, sex, column, expected, source)
    (2016, "f", "Prev", 0.343, "SHIMS2 prevalence 15-49, women"),
    (2016, "m", "Prev", 0.189, "SHIMS2 prevalence 15-49, men"),
    (2016, "f", "VLS", 0.747840, "SHIMS2 Table 9.4.A, women 74.8"),
    (2016, "m", "VLS", 0.623272, "SHIMS2 Table 9.4.A, men 62.3"),
    (2021, "f", "VLS", 0.886, "SHIMS3 Table 8.2, women 88.6"),
    (2021, "m", "VLS", 0.824, "SHIMS3 Table 8.2, men 82.4"),
    (2021, "m", "Inc", 0.20, "SHIMS3 Table 5.1, men 0.20"),
    (2021, "f", "Inc", 1.45, "SHIMS3 Table 5.1, women"),
    (2011, "m", "Inc", 1.65, "SHIMS1, men ~1.7"),
    (2011, "f", "Inc", 3.14, "SHIMS1, women ~3.1"),
]


def build():
    raw = pd.read_csv(RAW)
    out_rows = []
    for (year, agecat), g in raw.groupby(["Year", "AgeCat"], sort=False):
        for sex in ("m", "f"):
            row = dict(Year=int(year), AgeCat=agecat, sex=sex,
                       Gender=0 if sex == "m" else 1)   # target convention
            # Prev/VLS/ppdv block: raw Gender 0 is female
            src = g[g.Gender == (0 if sex == "f" else 1)]
            for c in FEMALE_IS_0:
                row[c] = float(src[c].iloc[0]) if len(src) else np.nan
            # Inc block: raw Gender 0 is male
            src = g[g.Gender == (0 if sex == "m" else 1)]
            for c in MALE_IS_0:
                row[c] = float(src[c].iloc[0]) if len(src) else np.nan
            # Indeterminate: carried on the Prev block's convention, flagged
            src = g[g.Gender == (0 if sex == "f" else 1)]
            for c in INDETERMINATE:
                row[c] = float(src[c].iloc[0]) if len(src) else np.nan
            out_rows.append(row)

    df = pd.DataFrame(out_rows)
    df["art_prev_convention"] = "INDETERMINATE - do not use without checking"

    for year, sex, col, expected, source in CHECKS:
        got = df[(df.Year == year) & (df.sex == sex)][col]
        assert len(got), f"missing {year} {sex}"
        assert abs(float(got.iloc[0]) - expected) < 1e-6, (
            f"{year} {sex} {col}: got {float(got.iloc[0])}, expected {expected} "
            f"({source}) -- the raw file's coding may have changed")

    # ppdv must still reproduce Prev x (1 - VLS) after the recode, which it can
    # only do if both came from the same original row.
    ok = df.ppdv.notna() & df.Prev.notna() & df.VLS.notna()
    resid = (df.loc[ok, "Prev"] * (1 - df.loc[ok, "VLS"]) - df.loc[ok, "ppdv"]).abs()
    assert (resid < 5e-4).all(), (
        f"ppdv no longer equals Prev x (1 - VLS) after recode; max residual "
        f"{resid.max():.2g} -- the block assignment is wrong")

    cols = (["Year", "AgeCat", "sex", "Gender", "Prev", "Prev_lb", "Prev_ub",
             "VLS", "VLS_lb", "VLS_ub", "ppdv", "ppdv_lb", "ppdv_ub",
             "Inc", "Inc_lb", "Inc_ub"] + INDETERMINATE
            + ["art_prev_convention"])
    return df[cols].sort_values(["Year", "sex"], ascending=[False, True])


if __name__ == "__main__":
    df = build()
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(df)} rows, all {len(CHECKS)} primary-source "
          f"checks passed\n")
    print(df[["Year", "sex", "Gender", "Prev", "VLS", "Inc"]].to_string(index=False))
