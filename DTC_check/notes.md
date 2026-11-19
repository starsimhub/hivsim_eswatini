# Calibration review notes — HIVSim Eswatini

*2026-10-27. Reviewing Adam Akullian's calibrated model: design row 868 from
[`experiments/024_hm_wave1`](../experiments/024_hm_wave1/SUMMARY.md), run
through [`build_and_check_calibration.py`](build_and_check_calibration.py).*

**Bottom line.** The testing and ART cascade was never a calibration target.
Issues 1–3 below all trace back to the model's HIV testing rates. These are
fixed values typed directly into the code, not estimated from data and not
adjusted by the calibration.
Issue 4 is a genuine calibration weakness.

---

## 1. What the calibration was fitted to

Wave 1 fitted **7 sexual-behaviour and transmission parameters only**:

- `beta_m2f`, `rel_beta_f2m`
- `s_f_young`
- `age_gap_shift`, `age_gap_sd_mult`
- `prop_f0`, `prop_m0`

These were scored against **48 targets**, built in
[`run.py` `build_observations()`](../experiments/024_hm_wave1/run.py#L122-L184).
History matching has no weights, only standard deviations (σ), so σ sets how
tightly each target binds.

### Group A — headline targets (8)

| Target | n | Source | σ |
|---|---|---|---|
| Prevalence 15–49 by sex, 2007 / 2011 / 2016 | 6 | `calibration_data/prevalence_by_age_sex.csv` | survey CI + 0.02 allowance for model error |
| Mean of those six | 1 | same | √(0.005² + 0.02²) |
| Peak AIDS deaths: 11,000 in 2004 | 1 | UNAIDS Spectrum, `data/eswatini_hiv_calib.csv` (`hiv.new_deaths`) | 2,000, deliberately loose |

> **Only the mean of the six prevalence values (`prev_15_49_all_mean`) was
> actually used to cut the parameter space.** All other targets were scored
> but not used to cut, and were used only to pick the best single draw.

The deaths σ is loose on purpose: the model peaks at about 7,000 deaths against
UNAIDS' 11,000.

### Group B — shape of the epidemic (16)

| Target | n | Source | σ |
|---|---|---|---|
| Female:male prevalence ratio, 15–49, per year | 3 | derived from the prevalence file | 15% relative |
| Female 15–24 : 35–44 prevalence ratio, per year | 3 | derived from the prevalence file | 15% relative |
| Incidence by sex (and by age band in 2016) | 8 | SHIMS1 2011 cohort (ages 18–49); SHIMS2 2016 Table 5.3.B | published CI, floored at 0.25 |
| Female:male incidence ratio, 2011 and 2016 | 2 | SHIMS1/2 (1.90, 2.04) | 0.35 |

### Group C — age-specific prevalence (24)

| Target | n | Source | σ |
|---|---|---|---|
| Prevalence by sex × 10-year age band, 2007 / 2011 / 2016 | 24 | prevalence file | survey CI + 0.02 |

The bands are 15–24, 25–34, 35–44 and 45+. The 2011 survey stops at 45–49, so
its top band is capped there.

### Notes on the sources

- **Prevalence surveys.** The repo labels all three years "PHIA":
  - **2011** is SHIMS1, and covers ages 15–49 only.
  - **2016** is SHIMS2 (Eswatini's PHIA), and covers ages 15–64.
  - **2007** is, *in my reading*, the 2006–07 Swaziland DHS, since PHIA didn't
    exist then. The repo doesn't say. It covers ages 15–64.
- **The incidence targets have changed since wave 1.**
  `calibration_data/incidence_by_age_sex.csv` now holds only the 4 sex-level
  totals. The age bands were dropped on 2026-09-04 because of false-recency
  bias (see the repo's `CLAUDE.md`). **The model's incidence age profile is now
  unconstrained by data.**
- **Held out for validation, never fitted:** 2021 prevalence and incidence
  (SHIMS3).
- **Explicitly excluded from fitting** (per `CLAUDE.md`): ART coverage, VMMC and
  everything in the testing/treatment cascade.

---

## 2. The four issues

### Issue 1 — Diagnosis far too high before 2020

**Cause.** Diagnosis is not fitted. It comes from hard-coded linear testing
ramps in [`interventions.py:18-23`](../interventions.py#L18-L23):

| Testing stream | Testing probability | Notes |
|---|---|---|
| General population | 0 (1990) → 0.50 (2020) → 0.60 (2040) | linear |
| Female sex workers | 0 → 0.75 → 0.85 | linear |
| CD4 < 200 | 0 → 0.85 → 0.95 | linear |
| Antenatal care | 0.9 | flat, first trimester |

These are assumptions, not data. With these ramps, diagnosis among men is
already 16% by 1995 and **94% in 2017, against SHIMS2's 77%**.

**Possible fix.**
- Make the testing ramp free calibration parameters: the ramp's start year,
  shape and level.
- Add first-90 targets: SHIMS2 2017 (men 77%, women 90%) and SHIMS3 2021
  (men 89%, women 95%).
- Use Kasturi's modelled series as a soft prior on the ramp's shape before 2016.

### Issues 2 and 3 — Men's on-ART-among-diagnosed too low in 2017, and the apparent "swap"

**Cause: the same as issue 1. This is not a data swap.**

- ART is **forced to match SHIMS coverage among *all* people living with HIV**,
  by age and sex, from `data/art_coverage.csv`.
- That table is set for 2004, 2011, 2016 and 2021 and interpolated in between.
- So the model's on-ART-among-diagnosed is simply the forced coverage divided
  by the model's own diagnosed fraction:

  *% on ART among diagnosed = ART coverage among PLHIV ÷ % diagnosed*

**Men, ages 15–49, 2017:**

| | % diagnosed | % on ART among diagnosed | Product (on ART among PLHIV) |
|---|---|---|---|
| Model | 0.938 | 0.788 | **0.74** |
| SHIMS2 | 0.773 | 0.884 | **0.68** |

The product is about right, because it is forced. The diagnosed fraction is
too high, so the second 90 is pushed down, and the two values look transposed.

Women diagnose less excessively (96% against 90%), so their second 90 lands on
the data (0.880 against 0.875).

**Possible fix.** Fixing the testing ramp (issue 1) should largely fix issues 2
and 3 as well.

### Issue 4 — Female prevalence too high at ages 20–40

**Cause: a genuine calibration weakness.**

- Only the overall 15–49 mean actually narrowed the parameter space. The
  age-band targets were used only to pick row 868, which was accepted as any
  draw within 3σ on all 48 targets.
- The 0.02 allowance for model error is about 90% of each prevalence σ, so a
  gap of ±2 percentage points was tolerated.
- Row 868's worst residual is **women 25–34 in 2016, at +3.0σ**.
- Its young-female susceptibility multiplier, `s_f_young`, is 2.22 on a
  0.8–3.0 range.

**Possible fix.** Adam's wave-2 plan
([`experiments/025_hm_wave2`](../experiments/025_hm_wave2/), not yet run)
already does two relevant things:
- cuts the allowance for model error to 0.01;
- fits the female:male and young:old prevalence ratios directly.

That is the natural place to address this issue.

---

## 3. A data-file problem

In `calibration_data/ppdv.csv`, **the incidence column has the sexes swapped**
relative to the prevalence and VLS columns.

For example, the 2021 `Gender = 1` row has men's prevalence (0.156) but women's
incidence (1.45). It is not a calibration target, but don't take incidence from
that file.
