# SUMMARY — experiment 008 (calibration targets)

**Status:** complete.
**Date:** 2026-05-07.

## What we have

89 target rows in [`outputs/calibration_targets.csv`](outputs/calibration_targets.csv),
spanning two quantities:

| Quantity | Rows | Source | Uncertainty |
|---|---|---|---|
| HIV prevalence (age × sex × year) | 54 | PHIA 2007, 2011, 2016 (`calibration_data/prevalence_by_age_sex.csv`) | **Sourced** — `lb`/`ub` from PHIA design |
| Annual AIDS deaths | 35 | UNAIDS Spectrum point estimates 1990–2024 (`data/eswatini_hiv_calib.csv → hiv.new_deaths`) | **Placeholder** — ±15% relative, applied uniformly |

Plot: [`outputs/calibration_targets.png`](outputs/calibration_targets.png).
Per-row audit: [`outputs/calibration_targets_audit.csv`](outputs/calibration_targets_audit.csv).

## Findings

### Prevalence (in good shape)

- 3 PHIA waves cover 1990s through 2016 — adequate to constrain epidemic
  peak, post-peak plateau, and post-ART-rollout decline.
- Male and female age-prevalence curves both show the expected forward
  shift in modal age between 2011 and 2016 (likely ART-driven survival
  effect).
- **Caveat: PHIA 2011 only covers ages 15–49.** PHIA 2007 and 2016
  both extend to 60–65. This means 2011 contributes 16 rows vs. 20 each
  from 2007 and 2016. Not a problem — just means the older-age trajectory
  is anchored only at 2007 and 2016.

### AIDS deaths (uncertainty needs work)

- Source file is annual point estimates with no bounds. **±15% is a
  placeholder consistent with typical Spectrum uncertainty**, not a
  sourced value.
- **Follow-up:** download the UNAIDS AIDSinfo CSV for Eswatini (it
  publishes lower/upper bounds explicitly), or extract from Spectrum
  output if you have access. Replace the placeholder with sourced bounds
  before the likelihood for deaths is finalised in `likelihood-design`.
- Until then, the placeholder lets us proceed with experiment 009 (prior
  predictive check) — coverage check is robust to small uncertainty
  miscalibration.

### Age distribution — not loaded; needs your decision

The user spec mentioned "age distribution from UNAIDS" as a possible
target, with a question mark. The available file (`data/eswatini_age_1985.csv`)
is the **1985 baseline used as the demographic-module initial condition,**
not a multi-year fitting target. Three options:

1. **Drop as fitting target.** Treat the population age structure as a
   qualitative validation check — at end of run, plot model age pyramid
   vs UN WPP for ~2020 and visually verify it's sensible. Recommended
   for the 1-month timeline.
2. **Add as fitting target.** Source year-stratified UN WPP age
   distribution for ~2010 and ~2020, add to the targets CSV. Adds
   complexity; demographic structure is not the primary thing the
   calibration parameters can move (most are HIV-specific).
3. **Defer.** Decide after seeing experiment 009 coverage check —
   if the model's demographic structure is obviously wrong, add as
   target then.

Default if not chosen: option 1 (drop, qualitative validation only).

### Validation hold-out

- `calibration_data/incidence_2021_VALIDATION_ONLY.csv` exists and is
  **not loaded** by `run.py`. Confirmed segregation.
- Also note `calibration_data/prevalence_2021_VALIDATION_ONLY.csv`
  exists alongside — likely a 2021 PHIA result. Decide whether 2021
  prevalence is also held out for validation, or whether to add it to
  the fitting set (more recent data, post-PEPFAR scaleup).

## Open follow-ups

- [ ] Source UNAIDS AIDSinfo bounds for `aids_deaths`; replace ±15%
      placeholder.
- [ ] Decide age-distribution target stance (drop / add / defer).
- [ ] Decide 2021 prevalence — fitting target or validation hold-out?
- [ ] Decide whether to include `hiv.n_infected` (PLHIV count) as a
      target — currently only its derivative form (prevalence_15_49) is
      implied through PHIA. UNAIDS publishes PLHIV with bounds.

## Next experiment

**009 — Prior predictive coverage check.** With this target set, run
50–200 sims drawn from the prior on the calibration parameters and check
whether the observed data falls inside the simulated ensemble. This is
the model-data sanity check before any calibration. See
`coverage-check` skill for the full protocol. This is also where we'll
need raccoon (it'll run faster on 120 cores than locally).
