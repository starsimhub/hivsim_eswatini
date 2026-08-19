# Experiment 005 — VMMC integration

## Goal

Add Voluntary Medical Male Circumcision (VMMC) to the model with age-stratified
coverage targets, fitted to PHIA prevalence rather than treated as a hazard.

## Why prevalence-target semantics

`stisim.VMMC` originally interpreted `p_vmmc` as a per-step hazard ("circumcise
this fraction of remaining uncircumcised men this timestep"), which doesn't map
cleanly onto PHIA-style cross-sectional prevalence data. ART has an analogous
mechanism (`art_coverage_correction()`) that *does* hit a prevalence target.
Mirrored that for VMMC: a new `prevalence_target=True` flag (default on) in the
patched stisim VMMC class, with a `_apply_prevalence_target()` method that tops
up circumcisions per (age_bin, sex) stratum to match the target proportion of
all alive males in that stratum. Never reduces (VMMC is irreversible).

Patch lives in `star_sim/stisim/stisim/interventions/hiv_interventions.py`
(local mod, candidate for upstream PR).

## Data triangulation

### Available local sources

1. **`/EMOD_eswatini/eSwatini2/InputFiles/Templates/campaign_PEPFAR_DISC_BASELINE.json`**
   — your previous EMOD calibration. Contains VMMC coverage targets at
   single-year ages 10/15/20/25/30/35/40/45/50/55, with `Times = [2004.5, 2013.5,
   2018.5, 2050]`. The commented-out `.Times = [2007, 2016, 2021, 2050]` reveals
   that the active `Times` were time-shifted ~2.5 years backward from the actual
   PHIA survey years.

2. **`external_data/mmc1.pdf`** — supplementary appendix to your 2020 Lancet HIV
   90-90-90 paper. Contains HIV prevalence by age/sex from SDHS 2007 + SHIMS 2011
   + SHIMS2 2016, but **no MMC tables**. Useful for HIV calibration only.

3. **`data/241123_SHIMS_ENG_RR3_Final-1.pdf`** — official SHIMS3 2021 final report.
   Table 12.5 has age-stratified male circumcision (medical + nonmedical) by
   5-year bin from 15-19 through 65+. **Authoritative source for 2021 values** —
   used to replace the EMOD JSON extrapolations.

4. **`data/pone.0298387.s005.docx`** — PLOS ONE supplement: cross-country MMC
   meta-analysis. Reference context, not a primary Eswatini source.

### Authoritative sources (need to verify against, when convenient)

- **SDHS 2006-07** — pre-program baseline; very low MMC (~5-10% in adult men).
- **SHIMS2 2016-17** — first major post-VMMC survey. ~30% adult MMC prevalence
  with strong age gradient (younger = higher due to school/youth-targeted
  programs).
- **SHIMS3 2021** — primary 2021 source. May report aggregate MMC for 15-49
  rather than full 5-year-bin stratification. Per recollection, ~50-60% overall.
- **PEPFAR/MOH operational data** — annual VMMC procedure counts by age band.
  Useful for cross-validating, but flow not stock.

## Recommended values (`data/vmmc_coverage.csv`)

2021 values are **SHIMS3 official Table 12.5** (medical + nonmedical, since the
model treats all circumcision as protective). 2007 and 2016 values are EMOD
JSON, which appear to be SDHS 2007 + SHIMS2 2016 derived.

| Age bin | 2007 | 2016 | 2021 | 2021 source |
|---|---|---|---|---|
| [10,15) | 0.04 | 0.44 | 0.85 | EMOD extrapolation (no SHIMS data <15) |
| [15,20) | 0.04 | 0.38 | 0.738 | SHIMS3 Table 12.5 |
| [20,25) | 0.07 | 0.31 | 0.604 | SHIMS3 Table 12.5 |
| [25,30) | 0.08 | 0.28 | 0.439 | SHIMS3 Table 12.5 |
| [30,35) | 0.10 | 0.25 | 0.345 | SHIMS3 Table 12.5 |
| [35,40) | 0.20 | 0.27 | 0.298 | SHIMS3 Table 12.5 |
| [40,45) | 0.12 | 0.15 | 0.325 | SHIMS3 Table 12.5 |
| [45,50) | 0.12 | 0.13 | 0.267 | SHIMS3 Table 12.5 |
| [50,55) | 0.12 | 0.16 | 0.252 | SHIMS3 Table 12.5 |
| [55,60) | 0.12 | 0.14 | 0.241 | SHIMS3 Table 12.5 |
| [60,65) | 0.12 | 0.14 | 0.166 | SHIMS3 Table 12.5 |

### SHIMS3 official totals (cross-checks)

- 15-24: 67.6% (medical 66.9 + nonmedical 0.7)
- **15-49: 48.3%** (medical 47.2 + nonmedical 1.1) — vs my CSV equal-weight
  average of [15,49) bins ≈ 41%. SHIMS3 weights by population (more young men
  due to age structure), so the official 48.3% is higher than the unweighted
  average — consistent.
- 50+: 19.5% (medical 16.3 + nonmedical 3.2)
- 15+ overall: 43.3%

### EMOD JSON had been close but slightly off for older bins

EMOD's extrapolated 2021 values were within ±5pp of the SHIMS3 actuals for ages
15-39, but underestimated 45-59 (~0.18-0.20 vs 0.24-0.27 actual) and the cohort
fix I applied ([40,45)=0.30) was almost exactly right (SHIMS3 actual 0.325).
Now using SHIMS3 directly, no projection needed.

### Where I retained EMOD values despite uncertainty

- **2007 [35,40) = 0.20** is anomalously high vs [30,35) = 0.10 and [40,45) = 0.12.
  Could reflect a genuine cohort with elevated traditional MC, or could be a
  data-entry oddity in the original SDHS analysis. Kept as-is — flag for review.
- **2007 [40,55] = 0.12** is suspiciously flat — looks like a placeholder
  (perhaps "no data, assume baseline"). Real SDHS values likely vary modestly.
  Kept at 0.12 pending verification against an actual SDHS report.

### Cohort-consistency check (population-weighted average for 15-49 in 2021)

Equal-weight across 7 bins [15,20)…[45,50): avg ≈ 41%. This is plausibly aligned
with SHIMS3 (recollection: 50–60%), though slightly low. Possible explanations:
(1) my recollection of SHIMS3 is inflated, (2) the 2016 EMOD source values are
slightly low vs SHIMS2 actuals, (3) the 2021 projections are conservative. Worth
cross-checking against the actual SHIMS3 report when available.

### Pre-program baseline (1990 anchor — added in followup)

`sc.smoothinterp` holds boundary values constant outside the data range, so the
earliest data point in the CSV becomes the "1985 baseline" by extrapolation.
With only 2007/2016/2021 anchors, this propagated the SDHS 2007 [35,40)=0.20
value back to 1985, giving an unrealistically high 20% MC baseline for that
cohort in the 1980s.

Fix: added a **1990 anchor row** with low values calibrated to SHIMS3 nonmedical
MMC (which captures lifetime traditional MC):

| Age bin | 1990 baseline |
|---|---|
| [10,15) | 0.005 |
| [15,20) | 0.01 |
| [20,25)–[25,30) | 0.02 |
| [30,35) | 0.03 |
| [35,40)–[60,65) | 0.04 |

Pre-2007 trajectory now linearly ramps from 1990 baseline to SDHS 2007 values.
Pre-1990 (rare in sim — start year is 1985) is held constant at 1990 baseline.

### SDHS 2007 validation

User provided actual SDHS 2007 circumcision-by-age data, which validates the
EMOD JSON 2007 values (and confirms the [35,40)=0.20 is real, not placeholder):

| Age | SDHS 2007 % | EMOD/CSV |
|---|---|---|
| 15-19 | 4.2 | 0.04 ✓ |
| 20-24 | 6.5 | 0.07 ✓ |
| 25-29 | 7.9 | 0.08 ✓ |
| 30-34 | 9.9 | 0.10 ✓ |
| 35-39 | **19.7** | **0.20 ✓** |
| 40-44 | 12.5 | 0.12 ✓ |
| 45-49 | 11.9 | 0.12 ✓ |

The 35-39 bump is real (cohort born 1968-1972 had elevated traditional MC; SDHS
shows 56.9% of those circumcisions happened pre-age-13). 50+ values are EMOD
extrapolations since SDHS 2007 didn't survey men 50+.

## Scope of experiment 005

- Apply stisim VMMC patch (prevalence-target semantics)
- Build vmmc_coverage.csv as above
- Enable VMMC in `interventions.py`
- Add a `VMMCPrevByAgeSex` analyzer to capture circumcision prevalence by age
  bin over time
- Replace the dashboard's panel F placeholder with a real VMMC coverage panel
  (lines = sim, points = data targets)
- Run 10-seed dashboards, compare male incidence (15–49) — VMMC should reduce
  it modestly from the 004 baseline since ~30%+ of adult men become circumcised
  with 60% acquisition reduction

## Out of scope (deferred)

- **Bellan acute-phase params** — was queued for 005 but split into 006
- **Re-calibration with new condom + VMMC** — VM-based work
- **Eff_circ as a calibration parameter** — use stisim default 0.6 (RCT-derived)
