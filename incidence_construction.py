"""Build calibration_data/incidence_by_age_sex.csv from SHIMS1 and SHIMS2.

Run `python incidence_construction.py` to regenerate.

Why incidence is worth adding as a target
-----------------------------------------
Everything the wave-1 transmission parameters govern is a *rate*, and prevalence
is a *stock* that confounds incidence with mortality, ART and cohort history.
Seven experiments (016-022) failed to move the age-shape defect in prevalence
precisely because susceptibility errors, mixing errors and mortality errors all
produce similar prevalence deviations. They produce different incidence
profiles, so incidence is the observation that discriminates among them.

Adding earlier-round incidence also makes the 2021 hold-out worth more, not
less: validating on 2021 after fitting 2011 and 2016 is a genuine out-of-period
prediction, where validating on it having never fitted incidence at all only
shows that the model cannot do incidence, without showing why.

Sources
-------
**2011 -- SHIMS1.** Ages 18-49, by sex. Taken from
`external_data/Swaziland_PPDV.csv` (recoded by `ppdv_construction.py`), whose
2021 and 2016 incidence values reproduce SHIMS3 Table 5.1 and SHIMS2 Table
5.3.B exactly, which is what licenses trusting its 2011 row. SHIMS1 measured
incidence with a **longitudinal cohort** rather than a cross-sectional recency
assay, and the intervals show it -- roughly half the width of the 2016 ones. It
is the best-measured incidence observation available and should carry the most
weight.

**2016-17 -- SHIMS2 Table 5.3.B**, "Annual HIV incidence using limiting
antigen/viral load/antiretroviral testing algorithm", by sex and age.

Two algorithms are published. 5.3.A uses LAg + viral load; **5.3.B adds ARV
biomarker testing**, which excludes people on treatment who could otherwise be
misclassified as recent infections. 5.3.B is used here for two reasons: it is
the more refined algorithm, and its 15-49 row (m 0.85, f 1.73) reproduces the
project's own PPDV file exactly, so it is already the series in use.

**2021 -- SHIMS3 Table 5.1 is deliberately NOT included.** It is the
designated validation hold-out
(`calibration_data/incidence_2021_VALIDATION_ONLY.csv`).

Cautions the likelihood must respect
------------------------------------
1. **Several cells are uninformative by construction.** SHIMS2 men 25-34 is
   1.50 with a 95% CI of (0.00, 3.06) -- an interval spanning the entire
   plausible range. Fitting the point estimate would chase noise. The published
   bounds are carried in this file so the likelihood can down-weight
   accordingly rather than treating all cells alike.
2. **Incidence and prevalence from the same survey are not independent.** They
   share a sampling frame, so treating them as independent likelihood
   components overstates the information. Worth an explicit down-weight.
3. **The absolute level is more assay-sensitive than the shape.** LAg-avidity
   estimates depend on the assumed mean duration of recent infection, which
   scales all cells roughly uniformly. The female:male ratio and the age
   profile are therefore more robust than the levels, which is why `fm_ratio`
   is emitted as a derived target alongside them.

What the data says, in the survey's own words
---------------------------------------------
SHIMS2's chapter conclusion: "New HIV infections continue at high rates among
males aged 25-34 years and females aged 35-49 years." Both are bands where this
model's prevalence residual is largest -- men 25-34 at -12.7 pp and women 35-44
at +4.5 pp -- which is the reason to bring incidence in.
"""

import pandas as pd

OUT = "calibration_data/incidence_by_age_sex.csv"

# (year, sex, age_low, age_high): (incidence_pct, lb, ub, source)
# Percentage annual incidence, i.e. per 100 person-years.
OBS = {
    # SHIMS1 2011, longitudinal cohort, ages 18-49. Sex-aggregated only -- no
    # published age breakdown available, so this enters as two rows.
    (2011, "m", 18, 50): (1.65, 1.28, 2.11, "SHIMS1 2011 (cohort)"),
    (2011, "f", 18, 50): (3.14, 2.63, 3.74, "SHIMS1 2011 (cohort)"),

    # SHIMS2 2016-2017, Table 5.3.B (LAg/VL/ARV algorithm)
    (2016, "m", 15, 25): (0.52, 0.00, 1.14, "SHIMS2 Table 5.3.B"),
    (2016, "m", 25, 35): (1.50, 0.00, 3.06, "SHIMS2 Table 5.3.B"),
    (2016, "m", 35, 50): (0.68, 0.00, 1.85, "SHIMS2 Table 5.3.B"),
    (2016, "f", 15, 25): (1.67, 0.62, 2.71, "SHIMS2 Table 5.3.B"),
    (2016, "f", 25, 35): (1.54, 0.00, 3.09, "SHIMS2 Table 5.3.B"),
    (2016, "f", 35, 50): (2.09, 0.23, 3.92, "SHIMS2 Table 5.3.B"),
}

# Published 15-49 aggregates. NOT emitted as fitting rows -- that would
# double-count the age bands above -- but used to build the F:M ratio target and
# to check the bands are consistent with their own published total.
AGG_15_49 = {
    (2011, "m"): 1.65, (2011, "f"): 3.14,     # 18-49 for SHIMS1
    (2016, "m"): 0.85, (2016, "f"): 1.73,
}


def build():
    rows = []
    for (year, sex, lo, hi), (val, lb, ub, src) in OBS.items():
        # A CI reaching zero means the cell cannot distinguish any incidence
        # from none. Flagged so a likelihood can drop or heavily down-weight it
        # rather than fitting a point estimate that carries no information.
        rows.append(dict(
            year=year, sex=sex, age_low=lo, age_high=hi,
            incidence_pct=val, lb=lb, ub=ub,
            ci_width=ub - lb, uninformative=(lb <= 0.0),
            source=src))
    df = pd.DataFrame(rows).sort_values(["year", "sex", "age_low"])

    # Derived: female:male incidence ratio at 15-49, more robust to the recency
    # assay's mean-duration assumption than the absolute levels are.
    ratios = []
    for year in sorted({y for y, _ in AGG_15_49}):
        m, f = AGG_15_49[(year, "m")], AGG_15_49[(year, "f")]
        ratios.append(dict(year=year, quantity="fm_ratio_15_49",
                           value=f / m,
                           note=f"female {f} / male {m}, published aggregates"))
    ratio_df = pd.DataFrame(ratios)

    # Consistency: each round's age bands should straddle its published total.
    for year in (2016,):
        for sex in ("m", "f"):
            bands = df[(df.year == year) & (df.sex == sex)].incidence_pct
            total = AGG_15_49[(year, sex)]
            assert bands.min() <= total <= bands.max(), (
                f"{year} {sex}: published total {total} outside its own band "
                f"range [{bands.min()}, {bands.max()}] -- check transcription")
    return df, ratio_df


if __name__ == "__main__":
    df, ratio_df = build()
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(df)} target rows "
          f"({int(df.uninformative.sum())} flagged uninformative)\n")
    print(df.to_string(index=False))
    print("\nderived targets:")
    print(ratio_df.round(3).to_string(index=False))
