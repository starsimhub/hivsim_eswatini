# Exp 027 — Cascade inputs and PrEP, then recalibration

> Daniel Citron, 2026-10-28. Plan A from `DTC_check/context_dtc.md`.
> Starts from 024's design row 868 (the current point fit).

## Question

Row 868 puts far too many people living with HIV in the "diagnosed" group
before 2020 (men 15–49: 94% in early 2017 against SHIMS2's 77%), because the
testing ramps in `interventions.py` were never fitted. Can fixed, sex-specific
general-population testing rates reproduce SHIMS2's first 90 by sex? Does
adding PrEP and a realistic ANC testing start leave the prevalence fit intact?
Then: can a recalibration of the 7 transmission parameters bring women's
20–40 prevalence down?

## Environment

All 027 runs: starsim 3.7.2 (`starsim` main @ `2a99fc85`, 2026-10-06) and
stisim 1.7.0 (`stisim` main @ `3730ddc`, 2026-10-01), Python 3.13 (conda env
`hivsim-dev`), editable checkouts beside this repo. NB row 868 was fitted on
stisim ~1.5.11 / starsim 3.5.2. Chosen parameters for each step are in
`params/` (copied out of the gitignored `outputs/`).

## Model changes (all opt-in; defaults reproduce 024–026 exactly)

- **Testing** (`interventions.get_testing_products(testing_pars)`):
  - general-population testing split into men's and women's streams, each the
    existing ramp times a fixed multiplier (`rel_test_m`, `rel_test_f`);
  - ANC testing starts in 1995, ramping 70% → 90% of pregnant women tested,
    midpoint 2005 (previously ~100% from 1985).
- **PrEP** (`prep_initiations.py`, data `data/prep_initiations.csv`):
  oral TDF (PrEPWatch, 75% efficacy) and LA-PrEP (MOH, 92%), both 1-year
  courses, driven by annual *initiation* counts spread over each year.
  Eligibility: HIV-negative, sexually active 15+, and FSW, FSW clients, more
  than one current partner, or no current stable partner.

## Decisions taken (context_dtc.md)

- SHIMS2 values are compared to the model at **2017.0** (fieldwork 2016–17).
- **2021 (SHIMS3) 90s and prevalence are held out** for checking only.
- Testing rates are fixed inputs, not calibrated parameters.

## Plan

1. Add PrEP and the ANC change. *(done: code above)*
2. `tune_testing.py`: with row 868 held fixed, grid over the testing
   multipliers; pick `rel_test_m`, `rel_test_f` so % diagnosed matches
   SHIMS2 (2017.0) by sex for ages 15–49 and 15+.
3. `check_inputs.py`: row 868 with old vs new inputs, 10 seeds — prevalence
   and incidence should barely move, since ART coverage is forced either way.
4. Re-run the calibration (025's design plus a women 25–34 target, SHIMS2
   targets at 2017.0) with the step-2 inputs fixed. Needs raccoon.
5. Re-check the 90s at the new point; repeat step 2 once if they drifted.

## Status 2026-10-28: steps 1–3 done, step 4 not started

**Extra change found during step 2.** `data/art_coverage.csv` and the VLS
input label SHIMS2 as 2016. Compared at 2017.0, the forced ART coverage had
already moved a fifth of the way to its 2021 value, which pushed the second
90 ~4 points above SHIMS2. `shims2_year=2017` (opt-in) moves both inputs to
2017, consistent with the decision above.

**Step 2 result** (`outputs/chosen_testing.json`; 10 seeds, at 2017.0):
`rel_test_m = 0.426`, `rel_test_f = 0.465`.

| | model | SHIMS2 |
|---|---|---|
| diagnosed, men 15–49 / 15+ | 77.4% / 80.7% | 77.3% / 80.1% |
| diagnosed, women 15–49 / 15+ | 90.5% / 91.1% | 90.2% / 90.7% |
| on ART of diagnosed, men 15–49 / 15+ | 89.6% / 87.9% | 88.4% / 90.2% |
| on ART of diagnosed, women 15–49 / 15+ | 89.7% / 90.0% | 87.5% / 88.1% |

All four first-90 values are inside SHIMS2's 95% CI; the second 90 is within
2.3 points without being fitted. No ART stratum is capped by the diagnosed
pool, though men 25–34 are close (96.6% of diagnosed on ART).
Held out (2021): men 15–49 diagnosed 86.6% vs SHIMS3 89.1%; women 94.1% vs 94.5%.

**Step 3 result** (`DTC_check/outputs/027/`): 15–49 prevalence by sex moves
≤ 0.34 points in 2007, 2011 and 2017; prevalence-fit MAE 0.042 → 0.040.
Criterion met.

**Open problems, before step 4:**
1. Diagnosis before ~2015 is still far above Kasturi's series (2011: men
   60% vs ~32%, women 81% vs ~36%). The ANC ramp (70% of pregnancies from
   1995) drives women's rise in 1995–2000; the general-population ramp,
   starting in 1990, drives the rest. Goal 3.1 needs a later testing start
   and/or a lower early ANC level — a decision for the researcher, since the
   ANC values were specified.
2. Within men, the age pattern is off: 15–24 diagnosed 70% vs 60%, 25–34
   62% vs 69%. Testing has no age dependence.
3. ART initiations run ~85,000/year after 2015, far above real new
   initiations: the stratified coverage correction stops and restarts people
   to hold each age band on target. Not a cascade-fit problem, but the
   counts panel shouldn't be read as programme initiations.

## Step 2b (2026-10-28): testing start years and ramps calibrated to SHIMS2

`calib_testing.py`: six testing inputs (general-population start year and
level by sex; ANC start year, starting share and ramp length), row 868
fixed, PrEP and `shims2_year=2017` on. **Targets: SHIMS2 first and second 90
at 2017.0, by sex, 15–49 and 15+** (σ from the Wilson CI). The third 90 is a
model input and is reported only. **Constraint:** ART must reach its forced
input — the share of forced ART not delivered (pooled over the ART-input
bins, worst of 2007.5/2011.5/2014.5) is penalised above 5%. 480 Latin-hypercube
points × 2 seeds; best 10 re-run on 8 fresh seeds.
Kasturi's series is a sanity check only (researcher, 2026-10-28).

Chosen (`outputs/calib_chosen.json`):

| | value |
|---|---|
| general-population testing starts | 2003.2 |
| level (× 0.5/yr in 2020), men / women | 0.60 / 0.57 |
| ANC testing starts | 1992.5 at 36%, reaching 90% by 1994.5 |
| forced ART not delivered (worst year) | 1.5% |

| early 2017 (target SHIMS2) | dx | on ART of dx | VLS of on ART (input) |
|---|---|---|---|
| men 15–49 | 78.8% (77.3) | 88.2% (88.4) | 90.3% (87.6) |
| men 15+ | 81.3% (80.1) | 87.2% (90.2) | 90.5% (90.5) |
| women 15–49 | 90.1% (90.2) | 89.9% (87.5) | 91.9% (91.4) |
| women 15+ | 90.6% (90.7) | 90.2% (88.1) | 91.8% (91.9) |

Sanity check, % diagnosed mid-2011: men 15–49 48% (Kasturi 32%), women 78%
(Kasturi 36%). Step 2's point (1990 start): men 60%, women 81%.

**SHIMS2 alone does not pin down the start years.** Step 2's point (testing
from 1990, ANC 1995 70%→90%) scores about as well on the same targets
(loss 21 on 4 seeds vs 27 on 8 fresh seeds; 2-seed noise in the search is
of order ±10). Fits differing by 10+ years in start date and by ~10 points
in 2011 % diagnosed are not distinguished by one survey year. The ANC result
(~90% of pregnancies from the mid-1990s) is close to the researcher's
original specification.

Knock-on (`DTC_check/outputs/027b/`, 10 seeds, vs row 868): prevalence 15–49
moves ≤ 0.5 points except men 2011 (−1.0); HIV deaths 2005–2011 up ~15%
(2008: 5,300 vs 4,400); prevalence-fit MAE 0.042 → 0.038.

Superseded: a first run also scored Kasturi's 2011 values (outputs
`calib_v1_kasturi_*`). It matched Kasturi but left 23–31% of the forced 2011
ART undeliverable — SHIMS1's ART coverage is about as large as Kasturi's
% aware. Kasturi was never meant as a target.

## Step 2c (2026-11-18): revised inputs, testing recalibrated

Revised instructions in `DTC_check/context_dtc.md` (2026-10-28 revision).
New inputs, all opt-in:

| input | 2c | was (2b) |
|---|---|---|
| ANC testing, per pregnancy | 0 in 1992 → 45% 2003 → 90% 2010, linear, flat after (`interventions.ANC_KNOTS_2C`); fixed | calibrated: 36% in 1992.5 → 90% 1994.5 |
| PMTCT | new `pmtct.py`: diagnosed pregnant women not on ART enrolled once per pregnancy, 0 in 2003 → 90% 2010; protects the unborn child (efficacy 0.96, as ART's `pmtct_efficacy`); mother's state untouched | none outside ART |
| VLS among on ART | logistic per sex through SHIMS2 (2017.0) and SHIMS3 (2021.0), floor 80%, ceiling 98% (`art_vls_coverage='phia_logistic'`) | flat to 2017, linear to 2021 |
| general-population testing | separate start year per sex (`gp_start_m/f`) | one start year |

"Max" for ANC and PMTCT is not given in the instructions; 90% is the
researcher's earlier ANC level. The VLS floor and ceiling are assumptions
(`vls_construction.to_vls_logistic`). Old note 6 (ANC 1995, 70→90%) is
superseded by goal 2.

**The ANC schedule sets women's early diagnosis.** With *no*
general-population testing for women, women 15–49 are already 34–37%
diagnosed in 2005 (Kasturi: 7%), linear or logistic ANC ramp alike: over
1992–2005 ANC makes ~5× more diagnoses than FSW and low-CD4 testing
combined, and they accumulate across a woman's pregnancies. So Kasturi's
level for women is unreachable under the specified ANC, and pulling men to
Kasturi's level would widen the men-vs-women gap the researcher objected to.

`calib_testing_2c.py`: four free inputs (`gp_start_m/f`, `rel_test_m/f`).
Loss = SHIMS2 first + second 90 at 2017.0 (as 2b) + 2b ART-delivery penalty
+ Kasturi guide, 15–49 at 2005.5/2008.5/2011.5/2014.5: men:women ratio of %
diagnosed (σ 0.10, the guide proper) and each sex's level (σ 0.20, a loose
pull toward later diagnosis). 420 LHS points × 2 seeds; best 10 on 8 fresh
seeds.

**Chosen** (`outputs/calib2c_chosen.json`): men start 1997.6, women 2002.8,
levels 0.52 / 0.53. Not the argmin: third by total loss (71.7 vs 68.4,
inside 8-seed noise) but the lowest Kasturi-guide loss. The argmin
(`calib2c_argmin.json`: men 2000.7, women 2008.9 at 1.13) puts women at
93.0% diagnosed in 2017 (outside SHIMS2's CI) and keeps men's on-ART-of-
diagnosed at 0.39 in 2005.

| 15–49 | 2005 | 2008 | 2011 | 2014 | 2017 (SHIMS2) |
|---|---|---|---|---|---|
| women dx (Kasturi) | 40% (7) | 61% (17) | 75% (36) | — | 89.4% (90.2) |
| men dx (Kasturi) | 27% (5) | 44% (14) | 58% (31) | — | 79.6% (77.3) |
| men:women ratio (Kasturi) | 0.68 (0.79) | 0.73 (0.82) | 0.78 (0.86) | 0.83 (0.91) | |
| on ART of dx, women / men | 0.16 / 0.23 | 0.33 / 0.43 | 0.49 / 0.57 | | 0.906 / 0.876 (0.875 / 0.884) |

The men-only spike in on-ART-of-diagnosed is gone (2b: men 0.77, women 0.17
in 2005). The cost is earlier diagnosis of men than 2b (2005: 27% vs 8%):
with women fixed high by ANC, balance can only come from raising men.

Dashboard (`DTC_check/outputs/027c/`, seeds 1–10): first 90 2017 women
91.0% (CI 88.7–91.5), men 80.9% (CI 74.1–80.2, just above); second 90
women 90.7% vs 87.5%, men 88.1% vs 88.4%. 15–49 prevalence moves ≤ 0.4
points from 2b, ≤ 0.6 from row 868; women 2016–2021 are +0.3–0.5 points
above 2b (lower early suppression → more transmission), the wrong direction
for goal 1. Prevalence MAE 0.039 (2b 0.038, row 868 0.042). Children 0–14
living with HIV in 2015: 14,800 vs 21,900 in 2b (PMTCT).

**Not addressed here:** goal 1 (women 20–40 prevalence) is step 4. Men's
diagnosis has no age dependence (open problem 2 above). SHIMS2 men 15–49
VLS (87.6%) is below the model (91.6%) because the VLS input is 15+ only.

## Step 2d (2026-11-19): age-structured testing and treatment

Instructions: `DTC_check/context_dtc.md`, "Next Iteration of calibration
2026-11-19". New inputs, all opt-in (`DTC_check` config `"027d"`):

| input | 2d | was (2c) |
|---|---|---|
| general-population testing | level per sex × cohort 0–24, 25–34, 35–49, 50+ (`testing_pars['rel_test_age']`, `interventions.AGE_COHORTS`); under-15s share the 0–24 rate | one level per sex |
| forced ART | same cohorts, from the SHIMS 90-90-90 file (dx × on ART of dx); 0–14 take 15–24 (`art_coverage='cohorts'`, `data/art_coverage_cohorts.csv`, `art_coverage_cohort_construction.py`) | [15,25) [25,35) [35,45) 45+, nothing under 15 |
| ANC testing | logistic per pregnancy: 0 before 1992.5, midpoint 2003, 95% of the 90% max by 2010 (`anc_logistic=True`); 1998: 10%, 2003: 45%, 2010: 86% | linear knots 1992 / 2003 / 2010 |
| VMMC | SHIMS2 Table 15.5.A medical + non-medical totals in the SHIMS2 row, moved to 2017 (`vmmc_shims2=True`) | EMOD-derived 2016 row |
| PMTCT, VLS | unchanged from 2c | |

Assumptions to note: SHIMS2 publishes no 50+ band, so the 2017 ART input for
50+ is backed out of the 15+ and 15–49 rows using unweighted n (men 88.5%
diagnosed × 95.0% on ART; women 93.2% × 91.0%). SHIMS1 (2011) ART is
remapped onto the new cohorts ([35,45) → 35–49, 45+ → 50+), not re-tabulated.

`calib_testing_2d.py`: 10 free inputs (start year per sex; level per sex ×
cohort). Loss as 2c but by age: SHIMS2 first + second 90 at 2017.0 for
15–24, 25–34, 35–49 and 15+ by sex; ART-delivery penalty on the cohort ART;
Kasturi guide (15–49 men:women ratio σ 0.10, levels σ 0.20). optuna TPE,
360 trials × 2 seeds, 2c's point enqueued; best 10 on 8 fresh seeds.

**Chosen** (`outputs/calib2d_chosen.json`), the argmin on fresh seeds:
men start 1995.4, women 2011.0; levels (0–24 / 25–34 / 35–49 / 50+)
men 0.18 / 0.57 / 0.31 / 2.28, women 0.99 / 1.25 / 1.38 / 0.68.
Loss 39.2 (SHIMS 25.7, Kasturi 13.5) vs 2c's point 60.1 (39.7, 20.4) on the
same seeds; forced ART not delivered: 2.5%.

| early 2017 | dx model [SHIMS2] | on ART of dx model [SHIMS2] |
|---|---|---|
| women 15–24 / 25–34 / 35–49 / 15+ | 76.4 / 91.5 / 96.7 / 91.1 [75.5 / 91.1 / 94.6 / 90.7] | 82.7 / 84.8 / 88.4 / 86.5 [83.7 / 85.5 / 90.4 / 88.1] |
| men 15–24 / 25–34 / 35–49 / 15+ | 68.3 / 65.0 / 82.7 / 80.7 [60.3 / 69.2 / 84.8 / 80.1] | 88.5 / 89.3 / 92.0 / 89.6 [92.2 / 84.3 / 89.9 / 90.2] |

| 15–49 % diagnosed (Kasturi) | 2005 | 2008 | 2011 | 2014 |
|---|---|---|---|---|
| women | 34 (7) | 51 (17) | 63 (36) | 80 (61) |
| men | 28 (5) | 43 (14) | 56 (31) | 68 (56) |
| men:women ratio | 0.82 (0.80) | 0.84 (0.83) | 0.89 (0.87) | 0.85 (0.91) |

The men:women balance now follows Kasturi; levels before ~2014 remain above
it. Women's general-population testing doesn't start until 2011, so women's
pre-2011 diagnosis is all ANC, FSW and low-CD4 testing: under the specified
ANC schedule (45% of pregnancies tested in 2003) the early level for women
can't come down further, and men are held near women by the ratio guide.

Dashboard (`DTC_check/outputs/027d/`, seeds 1–10): prevalence MAE 0.038
(2c 0.039, row 868 0.042); women 20–39 prevalence essentially unchanged
(2017: 46.5% vs 46.8% in 2c) — goal 1 still needs step 4. Children 0–14
living with HIV, 2017: 13,000 (2c 11,000, row 868 16,000) — children now
receive forced ART, so more survive. HIV deaths 2017: 1,300 (2c 2,100).

**Open, from 2d:**
1. Men 50+ are ~100% diagnosed in 2017 (backed-out SHIMS2 88.5%): no 2017
   target constrains 50+ directly, and level 2.28 is what fills the 50+ ART
   input. The backed-out 50+ could be added as a soft target.
2. VLS among those on ART is one curve per sex for all ages; SHIMS2 has men
   15–24 at 57.6% and women 15–24 at 81.2% vs the model's ~91%.
3. Diagnosis before 2014 is still 15–25 points above Kasturi for both sexes;
   the remaining lever is the ANC schedule.

## Success criteria (set before looking)

- Step 2: first 90 within SHIMS2's 95% CI for each sex, at 15–49 and 15+.
- Second 90 within 3 points of SHIMS2 without being fitted, and no stratum
  where ART falls short of its forced target for lack of diagnosed people.
- Step 3: 15–49 prevalence by sex moves less than 1 percentage point in
  2007, 2011 and 2017.
