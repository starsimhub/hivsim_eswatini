"""Build calibration_data/incidence_by_age_sex.csv from SHIMS1 and SHIMS2.

Run `python incidence_construction.py` to regenerate.

Revision, 2026-09-04
--------------------
**Only the adult aggregates are fitted now.** The six age-banded rows are
retained with `fit=False`. See FIT_POLICY below for the false-recency argument
that motivated the change, and note that caution 3 below has been REVERSED as a
result -- it originally claimed the age profile was more robust than the levels,
which is the opposite of what the prevalence-amplification arithmetic implies.

Why incidence is worth adding as a target
-----------------------------------------
Everything the wave-1 transmission parameters govern is a *rate*, and prevalence
is a *stock* that confounds incidence with mortality, ART and cohort history.
Seven experiments (016-022) failed to move the age-shape defect in prevalence
precisely because susceptibility errors, mixing errors and mortality errors all
produce similar prevalence deviations. They produce different incidence
LEVELS, so incidence still discriminates among them -- but via the level and the
sex ratio, not via the age profile, which the 2026-09-04 revision drops.

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
3. **REVERSED 2026-09-04. The age profile is LESS trustworthy than the level,
   not more.** The original claim here was that because LAg-avidity estimates
   depend on an assumed mean duration of recent infection (MDRI) that scales all
   cells roughly uniformly, the ratios and the age profile were more robust than
   the levels. The MDRI part is right and is why `fm_ratio` is still emitted.
   The age-profile part was wrong: it accounted only for MDRI, and ignored
   false-recency bias, which does NOT scale uniformly -- it scales with
   P/(1-P), so it grows with age exactly as prevalence does. See FIT_POLICY.

What the data says, in the survey's own words
---------------------------------------------
SHIMS2's chapter conclusion: "New HIV infections continue at high rates among
males aged 25-34 years and females aged 35-49 years."

Read that with caution. Both are high-prevalence bands, which is precisely where
false-recency bias inflates a cross-sectional recency estimate most, so the
conclusion is partly a property of the assay rather than of the epidemic. The
2011 cohort estimate does not share this problem -- a longitudinal design has no
recency misclassification -- which is a further reason to weight it above 2016.
"""

import numpy as np
import pandas as pd

OUT = "calibration_data/incidence_by_age_sex.csv"
# The six age-banded rows are NOT calibration targets and are written here
# instead, under a name that cannot be mistaken for one. They stay available for
# the age-profile figure and for the false-recency argument.
OUT_REF = "data/shims2_incidence_by_age_REFERENCE_NOT_A_TARGET.csv"

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

# Published adult aggregates: (value, lb, ub). 2011 is 18-49 (SHIMS1's published
# range), 2016 is 15-49. CIs from data/eswatini_ppdv.csv, which reproduces
# SHIMS3 Table 5.1 and SHIMS2 Table 5.3.B exactly.
#
# **These are the fitting rows as of 2026-09-04** -- see FIT_POLICY below. They
# used to be excluded to avoid double-counting the bands; now it is the bands
# that are excluded, so the double-counting concern is satisfied from the other
# direction.
AGG_15_49 = {
    (2011, "m"): (1.65, 1.28, 2.11), (2011, "f"): (3.14, 2.63, 3.74),
    (2016, "m"): (0.85, 0.21, 1.49), (2016, "f"): (1.73, 0.96, 2.50),
}

# SUPERSEDED 2026-09-04 (same day). This was an aggregation sigma of f 0.126,
# m 0.194, added because SHIMS1 publishes 18-49 while PopByAgeSex works in
# 5-year bands, so 18-49 looked bracketable but not expressible.
#
# It IS expressible. Pro-rate the 15-19 band to 2/5 of its person-time and the
# model lands on the survey's own basis, needing only that incidence is roughly
# flat across 15-19 -- no assumption about 15-17 incidence, because the survey
# excludes them and so does the pro-rated estimate. Model 18-49 then comes out
# at f 3.11 and m 1.70 against targets of 3.14 and 1.65, i.e. z = -0.11 and
# +0.19, against -0.04 and -0.49 on the naive 15-49 basis.
#
# So the sigma widening is removed and the published CI stands alone. The
# `model_age_basis` column tells the consumer how to compute the model side.
# Note the direction of travel: putting the MODEL on the DATA's basis makes the
# targets sharper, where adjusting the data to the model's basis would have
# required inventing a 15-17 incidence from either the FRR-contaminated SHIMS2
# 15-24 cell or from the model itself, which would be circular.
BAND_PRORATE = 0.4       # 18-19 is 2/5 of the 15-19 band

FIT_POLICY = """\
Fitted rows are the adult aggregates only; the age-banded rows are retained for
provenance with fit=False.

Researcher decision, 2026-09-04. The published female age profile is essentially
flat (1.67, 1.54, 2.09 across 15-24, 25-34, 35-49), which is not a credible
epidemiological shape in this setting -- and it is the expected signature of
false-recency bias, which is amplified by prevalence. A recency assay credits
some fraction FRR of long-standing infections as recent, and because those are
counted against a SUSCEPTIBLE denominator the spurious incidence scales as
FRR x P/(1-P)/MDRI. At MDRI = 130 d that term is 5x larger in women 35-49 than
in women 15-24, and 18x larger in men. At FRR = 0.3% it accounts for ~40% of the
reported value in both older female bands and ~95% in men 35-49.

So the age SHAPE of these estimates is not trustworthy even where the CI looks
usable, and 4 of the 6 banded rows already had CIs reaching zero. Aggregating to
15-49 does not remove the bias but turns it into a population-weighted average
rather than something concentrated in the band we would most want to learn from.

Keeps incidence level as a target; abandons its age shape. The model's own age
profile (peak at 22.5 in women, 32.5 in men) is therefore unconstrained by data,
which is recorded in exp 024 rather than hidden.
"""


def _model_basis(year):
    """How the model side must be computed to match this row's published range.

    Consumed by run.py. Spelling it out in the file beats a convention nobody
    reads: run.py previously compared a model 15-50 estimate to an 18-49 target
    and nothing in the data said it was wrong.
    """
    return "18_50_prorate_first_band_0.4" if year == 2011 else "15_50"


def build():
    """Returns (targets, reference, ratios). Targets are aggregates ONLY."""
    tgt = []
    for (year, sex), (val, lb, ub) in AGG_15_49.items():
        lo, hi = (18, 50) if year == 2011 else (15, 50)
        tgt.append(dict(
            year=year, sex=sex, age_low=lo, age_high=hi,
            incidence_pct=val, lb=lb, ub=ub, ci_width=ub - lb,
            uninformative=(lb <= 0.0), sigma=(ub - lb) / 3.92,
            model_age_basis=_model_basis(year),
            source=("SHIMS1 2011 (longitudinal cohort), 18-49" if year == 2011
                    else "SHIMS2 Table 5.3.B, 15-49")))
    df = pd.DataFrame(tgt).sort_values(["year", "sex"])

    # Not targets. Kept so the age profile stays plottable and the FRR argument
    # auditable, under a filename that cannot be mistaken for a target sheet.
    ref = []
    for (year, sex, lo, hi), (val, lb, ub, src) in OBS.items():
        if hi - lo >= 25:
            continue                       # the aggregates, already emitted
        ref.append(dict(
            year=year, sex=sex, age_low=lo, age_high=hi,
            incidence_pct=val, lb=lb, ub=ub, ci_width=ub - lb,
            uninformative=(lb <= 0.0), source=src,
            excluded_because="age shape not trusted; see FIT_POLICY"))
    ref_df = pd.DataFrame(ref).sort_values(["year", "sex", "age_low"])

    ratios = []
    for year in sorted({y for y, _ in AGG_15_49}):
        m, f = AGG_15_49[(year, "m")][0], AGG_15_49[(year, "f")][0]
        ratios.append(dict(year=year, quantity="fm_ratio_15_49", value=f / m,
                           note=f"female {f} / male {m}, published aggregates"))
    ratio_df = pd.DataFrame(ratios)

    # The excluded bands should still straddle their own published total -- a
    # transcription check that survives the exclusion.
    for year in (2016,):
        for sex in ("m", "f"):
            bands = ref_df[(ref_df.year == year) & (ref_df.sex == sex)].incidence_pct
            total = AGG_15_49[(year, sex)][0]
            assert bands.min() <= total <= bands.max(), (
                f"{year} {sex}: published total {total} outside its own band "
                f"range [{bands.min()}, {bands.max()}] -- check transcription")
    assert len(df) == 4, f"expected 4 fitting rows, got {len(df)}"
    assert not (df.age_high - df.age_low < 25).any(), "an age band leaked in"
    return df, ref_df, ratio_df


if __name__ == "__main__":
    df, ref_df, ratio_df = build()
    print(FIT_POLICY)
    df.to_csv(OUT, index=False)
    ref_df.to_csv(OUT_REF, index=False)
    print(f"wrote {OUT}: {len(df)} FITTING rows (sex-specific aggregates only)")
    print(df.to_string(index=False))
    print(f"\nwrote {OUT_REF}: {len(ref_df)} rows, NOT calibration targets")
    print(ref_df[["year", "sex", "age_low", "age_high", "incidence_pct",
                  "lb", "ub", "uninformative"]].to_string(index=False))
    print("\nderived (not emitted as a target row):")
    print(ratio_df.round(3).to_string(index=False))
