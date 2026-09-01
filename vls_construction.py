"""Build data/eswatini_vls.csv from the SHIMS2 and SHIMS3 primary sources.

Every value here is transcribed by hand from a numbered table in a survey final
report, so the transcription is auditable against the source rather than
inherited from a derived file. Run `python vls_construction.py` to regenerate.

Sources
-------
**SHIMS3 2021** — `data/241123_SHIMS_ENG_RR3_Final-1.pdf` (in this repo)
  - Table 8.1  VLS by demographic characteristics, incl. HIV diagnosis and
               treatment status. Report p.69.
  - Table 8.2  VLS by age and sex, 5-year bands. Report p.71.

**SHIMS2 2016-2017** — SHIMS2 Final Report, PHIA Project / ICAP at Columbia,
  https://phia.icap.columbia.edu/wp-content/uploads/2019/05/SHIMS2_Final-Report_05.03.2019_forWEB.pdf
  - Table 9.3.A  VLS by demographic characteristics, incl. self-reported
                 diagnosis and ART status. Report p.71.
  - Table 9.4.A  VLS by age and sex, 5-year bands. Report p.72.

Two distinct quantities, and conflating them is the trap
--------------------------------------------------------
**`vls_among_plhiv`** — the headline PHIA indicator (Global AIDS Monitoring
1.3/1.4): of everyone living with HIV, what fraction is virally suppressed.
This is a *cascade outcome*: it rises when either ART coverage or suppression
among the treated rises. Useful for validation. **Not** the model input.

**`vls_given_art`** — of those on ART, what fraction is suppressed. This is the
quantity `sti.ART(vls_coverage=...)` wants: "fraction of newly-initiated agents
who achieve viral suppression (effective ART) rather than non-suppressive ART".
The model already takes ART coverage as an input from `data/art_coverage.csv`,
so feeding it the population figure would double-count the coverage ramp.

The numbers make the distinction stark. Among all PLHIV aged 15-49, suppression
went 62-75% (2016) to 82-89% (2021); among those *on ART* it went 91-92% to
96-97%. Almost the whole population-level gain is the coverage ramp the model
already has. The model's actual error from defaulting `vls_coverage` to 1.0 is
the missing 4-9% of treated patients who are not suppressed -- who transmit at
`nonsupp_art_efficacy` = 0.35 rather than `effective_art_efficacy` = 0.99.

Comparability caveat between rounds
-----------------------------------
SHIMS2 classified ART status by **self-report** ("Previously diagnosed, on
ART"). SHIMS3 classified it by self-report **adjusted by ARV biomarker
testing** ("Aware of HIV status and on ART"). The biomarker adjustment
reclassifies some self-reported-untreated people as treated, which mechanically
raises the measured denominator and can move `vls_given_art` in either
direction. The 91.9 -> 96.2 rise is therefore not a clean like-for-like trend.

1985-2011 is not covered by either survey
-----------------------------------------
SHIMS1 (2011) was primarily an incidence survey and does not provide a VLS
table comparable to these. `external_data/Swaziland_PPDV.csv` carries a 2011
VLS of 0.306 (F) / 0.297 (M) among PLHIV for ages 18-49, which is consistent
with ART coverage of ~0.33 multiplied by a ~0.92 suppression rate -- i.e. it
looks derived rather than measured. It is deliberately NOT included here.
Any pre-2016 `vls_given_art` used by a model run is an assumption, and the
experiment that makes it should say so rather than hiding it in this file.

NB `vls_coverage` defaults any stratum it is not given to 100% suppression, so
a partial table silently reintroduces the very default it is meant to replace.
"""

import pandas as pd

OUT = "data/eswatini_vls.csv"

# --- VLS among all PLHIV, by 5-year age band -------------------------------
# (value_pct, n_unweighted). Parenthesised estimates in the source -- based on
# 25-49 unweighted cases -- are flagged small_denominator=True. Asterisked
# estimates (<25 cases) are suppressed in the source and omitted here.
PLHIV = {
    # SHIMS2 2016-2017, Table 9.4.A
    (2016, "m"): {(15, 20): (40.2, 40, True), (20, 25): (25.8, 29, True),
                  (25, 30): (48.2, 81, False), (30, 35): (58.9, 139, False),
                  (35, 40): (65.3, 184, False), (40, 45): (74.3, 133, False),
                  (45, 50): (77.4, 123, False), (50, 55): (83.7, 84, False),
                  (55, 60): (86.8, 58, False), (60, 65): (92.3, 53, False),
                  (65, 100): (84.9, 48, True)},
    (2016, "f"): {(15, 20): (55.0, 72, False), (20, 25): (55.7, 197, False),
                  (25, 30): (71.5, 319, False), (30, 35): (75.3, 394, False),
                  (35, 40): (80.7, 319, False), (40, 45): (87.0, 224, False),
                  (45, 50): (80.4, 169, False), (50, 55): (83.2, 135, False),
                  (55, 60): (87.7, 91, False), (60, 65): (86.5, 67, False),
                  (65, 100): (86.6, 44, True)},
    # SHIMS3 2021, Table 8.2
    (2021, "m"): {(15, 20): (83.9, 28, True), (20, 25): (77.4, 28, True),
                  (25, 30): (60.2, 33, True), (30, 35): (63.8, 88, False),
                  (35, 40): (84.5, 130, False), (40, 45): (89.5, 137, False),
                  (45, 50): (92.2, 138, False), (50, 55): (96.8, 99, False),
                  (55, 60): (89.2, 65, False), (60, 65): (93.1, 66, False),
                  (65, 100): (95.6, 62, False)},
    (2021, "f"): {(15, 20): (74.4, 54, False), (20, 25): (76.7, 139, False),
                  (25, 30): (78.0, 251, False), (30, 35): (91.9, 324, False),
                  (35, 40): (92.0, 342, False), (40, 45): (94.2, 278, False),
                  (45, 50): (96.2, 187, False), (50, 55): (97.5, 157, False),
                  (55, 60): (95.0, 127, False), (60, 65): (92.7, 60, False),
                  (65, 100): (97.6, 92, False)},
}

# Aggregate rows carried through as published rather than recomputed, since the
# survey weights are not in the report. SHIMS2 Table 9.4.A; SHIMS3 Table 8.2.
PLHIV_TOTALS = {
    (2016, "m", (15, 50)): (62.3, 729), (2016, "f", (15, 50)): (74.8, 1694),
    (2016, "m", (15, 100)): (67.6, 972), (2016, "f", (15, 100)): (76.0, 2031),
    (2021, "m", (15, 50)): (82.4, 582), (2021, "f", (15, 50)): (88.6, 1575),
    (2021, "m", (15, 100)): (86.1, 874), (2021, "f", (15, 100)): (90.1, 2011),
}

# --- VLS conditional on being on ART -- the model input --------------------
# SHIMS2 Table 9.3.A row "Previously diagnosed, on ART" (self-report basis).
# SHIMS3 Table 8.1 row "Aware of HIV status and on ART" (ARV-biomarker
# adjusted). Both are among adults 15+; neither report stratifies this row by
# age, which is why the model input has no age dimension.
GIVEN_ART = {
    (2016, "m"): (91.3, 694), (2016, "f"): (92.2, 1587),
    (2021, "m"): (96.7, 775), (2021, "f"): (95.9, 1885),
}

# --- VLS among all PLHIV, published 10-year aggregates ---------------------
# These bands match data/art_coverage.csv's bins exactly, so the derived
# suppression-given-ART needs no re-weighting of 5-year bands.
# SHIMS2 Table 9.4.B; SHIMS3 Table 8.2 (its aggregate rows).
PLHIV_10YR = {
    (2016, "m"): {(15, 25): 32.9, (25, 35): 54.8, (35, 45): 69.2,
                  (45, 55): 80.0, (55, 65): 89.3},
    (2016, "f"): {(15, 25): 55.5, (25, 35): 73.5, (35, 45): 83.4,
                  (45, 55): 81.5, (55, 65): 87.3},
    (2021, "m"): {(15, 25): 80.5, (25, 35): 62.9, (35, 45): 87.2,
                  (45, 55): 94.2, (55, 65): 90.8},
    (2021, "f"): {(15, 25): 76.1, (25, 35): 85.7, (35, 45): 93.0,
                  (45, 55): 96.8, (55, 65): 94.3},
}

SURVEY = {2016: "SHIMS2 2016-2017", 2021: "SHIMS3 2021"}
SRC_PLHIV = {2016: "SHIMS2 Table 9.4.A", 2021: "SHIMS3 Table 8.2"}
SRC_ART = {2016: "SHIMS2 Table 9.3.A", 2021: "SHIMS3 Table 8.1"}


def build():
    rows = []
    for (year, sex), bands in PLHIV.items():
        for (lo, hi), (pct, n, small) in bands.items():
            rows.append(dict(
                year=year, survey=SURVEY[year], sex=sex, age_low=lo, age_high=hi,
                measure="vls_among_plhiv", value=pct / 100.0, n_unweighted=n,
                small_denominator=small, source_table=SRC_PLHIV[year]))

    for (year, sex, (lo, hi)), (pct, n) in PLHIV_TOTALS.items():
        rows.append(dict(
            year=year, survey=SURVEY[year], sex=sex, age_low=lo, age_high=hi,
            measure="vls_among_plhiv", value=pct / 100.0, n_unweighted=n,
            small_denominator=False, source_table=SRC_PLHIV[year]))

    for (year, sex), (pct, n) in GIVEN_ART.items():
        rows.append(dict(
            year=year, survey=SURVEY[year], sex=sex, age_low=15, age_high=100,
            measure="vls_given_art", value=pct / 100.0, n_unweighted=n,
            small_denominator=False, source_table=SRC_ART[year]))

    df = pd.DataFrame(rows).sort_values(
        ["measure", "year", "sex", "age_low"]).reset_index(drop=True)

    # Internal consistency: the published 15-49 aggregate must sit inside the
    # range of its own 5-year bands. Catches a transcription slip that a
    # spot-check by eye would miss.
    for year in (2016, 2021):
        for sex in ("m", "f"):
            bands = [v[0] / 100 for (lo, hi), v in PLHIV[(year, sex)].items()
                     if lo >= 15 and hi <= 50]
            tot = PLHIV_TOTALS[(year, sex, (15, 50))][0] / 100
            assert min(bands) <= tot <= max(bands), (
                f"{year} {sex}: published 15-49 total {tot} outside its own "
                f"band range [{min(bands)}, {max(bands)}] -- check transcription")
    return df


def derive_given_art_by_age(art_coverage_path="data/art_coverage.csv"):
    """Suppression-given-ART by age and sex, as VLS-among-PLHIV / ART coverage.

    The age gradient in *population* suppression could be a gradient in
    coverage, a gradient in suppression among the treated, or both. PHIA
    publishes the first two quantities but only reports suppression-given-ART by
    sex, so the age split has to be derived.

    Sex convention: `data/art_coverage.csv` follows stisim's, which
    `parse_coverage` documents as **0 = female, 1 = male** (utils.py). Note this
    is the *opposite* of `calibration_data/prevalence_by_age_sex.csv`, where
    0 = male. Four files in this project use two different conventions; the
    inversion is the bug that cost exp 018 its rare-event floor estimate. Here
    the mapping is applied once, explicitly, and then everything is keyed on
    'm'/'f' strings.

    The derivation carries its own validation: collapsed over age it must
    reproduce the directly measured `vls_given_art`. If it does not, either a
    transcription is wrong or the two surveys' denominators disagree.
    """
    art = pd.read_csv(art_coverage_path)
    art["sex"] = art.Gender.map({0: "f", 1: "m"})   # stisim convention
    art["age_low"] = art.AgeBin.str.extract(r"\[(\d+),").astype(int)

    rows = []
    for (year, sex), bands in PLHIV_10YR.items():
        for (lo, hi), pct in bands.items():
            a = art[(art.Year == year) & (art.sex == sex) & (art.age_low == lo)]
            if not len(a) or a.p_art.iloc[0] <= 0:
                continue
            p_art = float(a.p_art.iloc[0])
            ratio = (pct / 100.0) / p_art
            # data/art_coverage.csv has no measured 45+ band: it copies the
            # [35,45) value into [45,100), because the source survey table only
            # reports three bands. VLS keeps rising after 45 while the copied
            # coverage does not, so the ratio there exceeds 1 -- more suppressed
            # people than treated ones, which is impossible. Flagged rather than
            # silently clipped.
            extrapolated = lo >= 45
            rows.append(dict(year=year, sex=sex, age_low=lo, age_high=hi,
                             vls_among_plhiv=pct / 100.0, p_art=p_art,
                             vls_given_art=ratio,
                             art_coverage_extrapolated=extrapolated,
                             implausible=ratio > 1.0))
    return pd.DataFrame(rows)


def to_vls_coverage(df, years=None, fill_back_to=None):
    """Reshape to the DataFrame `sti.ART(vls_coverage=...)` expects.

    Columns Year / AgeBin / Gender / p_vls, from the `vls_given_art` rows.

    `fill_back_to` prepends a row at that year holding the earliest observed
    value, because the surveys start in 2016 and the model starts in 1985.
    That is an ASSUMPTION -- suppression among the treated is taken as flat
    before the first measurement -- and the caller should record it as one.
    Left as None by default so it cannot be applied silently.
    """
    art = df[df.measure == "vls_given_art"]
    rows = []
    for _, r in art.iterrows():
        rows.append(dict(Year=r.year, AgeBin="[15,100)", Gender=r.sex,
                         p_vls=r.value))
    if fill_back_to is not None:
        first = art[art.year == art.year.min()]
        for _, r in first.iterrows():
            rows.append(dict(Year=fill_back_to, AgeBin="[15,100)",
                             Gender=r.sex, p_vls=r.value))
    out = pd.DataFrame(rows).sort_values(["Year", "Gender"])
    if years is not None:
        out = out[out.Year.isin(years)]
    return out.reset_index(drop=True)


if __name__ == "__main__":
    df = build()
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(df)} rows")
    print(df[df.measure == "vls_given_art"].to_string(index=False))
    print("\nas vls_coverage input (back-filled to 1985 -- an assumption):")
    print(to_vls_coverage(df, fill_back_to=1985).to_string(index=False))
