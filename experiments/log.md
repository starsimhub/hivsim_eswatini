# Model Experiment Log

Running record of changes to the Eswatini HIV model, their rationale, and observed effects.

---

## 001 — Stratified ART coverage by age, sex, and year
- **Date**: 2026-04-07
- **Commit**: `7b0811d`
- **Branch**: adam/model-dev

### What changed
- ART input switched from `data/n_art.csv` (national total counts + future_coverage=90%) to `data/art_coverage.csv` (proportion of PLHIV on ART by age bin, gender, and year from PHIA surveys)
- Code: `interventions.py` — `make_interventions()` now uses `sti.ART(coverage=art_data)`

### Why
- The unstratified input applied the same ART coverage to all age/sex groups, but PHIA surveys show large differences (e.g., in 2016: 90% for females 35-45 vs 36% for males 15-25)
- Stratified input allows the model to reproduce age/sex-specific treatment patterns

### Results
- **Fig 1 (HIV calibration)**: Panel C (People on ART) shows lower ART numbers post-2015. Other panels broadly similar.
- **ART by age/sex plot**: Model captures general scale-up trajectory but does not fully reproduce the male-female gap, especially in 25-35 age group. See experiment 002 for root cause.
- **Incidence/prevalence**: Modest changes. Tighter uncertainty bands post-2015.

### Figures
- `figures/archive/2026-04-06_before_stratified_art/` (baseline)
- `figures/archive/2026-04-07_after_stratified_art/` (after change)

### Decision
Keep change. Proceed to fix underlying allocation bug (experiment 002).

---

## 002 — Fix ART allocation to respect age/sex strata
- **Date**: 2026-04-14
- **Commit**: (pending)
- **Branch**: adam/model-dev

### What changed
- Fixed `_get_n_to_treat_stratified()` and `art_coverage_correction()` in stisim source to allocate ART **within each (age bin x sex) stratum** rather than computing a single global total
- Previously, stisim summed stratum targets into one number and then allocated globally by CD4 priority, ignoring which stratum agents belonged to

### Why
- The model was not reproducing its own ART coverage inputs: older/sicker agents were over-treated at the expense of younger groups
- The male-female ART gap was smoothed out because prioritization ignored sex
- Root cause identified by reading `stisim/interventions/hiv_interventions.py` lines 434-459 and 533-556

### Results
- **Dramatic improvement** in ART coverage by age/sex fit
- Before fix: male and female ART curves nearly overlapped (global CD4 allocation washed out sex differential)
- After fix: clear separation between male and female curves matching PHIA data
  - Ages 15-25: Female ~60% vs Male ~36% in 2016 (matches PHIA)
  - Ages 25-35: Female ~80% vs Male ~60% in 2016 (matches PHIA)
  - Ages 35-45: Both reach 90%+ by 2021, gap narrows (matches PHIA)
- Uncertainty bands much tighter — per-stratum allocation is more deterministic
- Model now reproduces its own stratified ART inputs

### Figures
- `experiments/002_art_strata_fix/before_fix.png` — pre-fix (experiment 001 baseline)
- `experiments/002_art_strata_fix/after_fix.png` — post-fix
- `figures/art_coverage_by_age_sex.png` — current (post-fix)

### Upstream
- Bug filed as GitHub issue on starsimhub/stisim (see `experiments/002_art_strata_fix/github_issue.md`)
- Local fix applied to `star_sim/stisim/stisim/interventions/hiv_interventions.py`

### Decision
Keep fix. Major improvement in reproducing stratified ART inputs. Submit PR to stisim upstream.

---

## 003 — Lower sexual debut age to DHS eSwatini values
- **Date**: 2026-04-14
- **Commit**: (pending)
- **Branch**: adam/model-dev

### What changed
- Added `debut_pars_f=[17.5, 2.5]` and `debut_pars_m=[18.5, 2.5]` to `StructuredSexual` in `run_sims.py`
- Previously using stisim package defaults: F=20yr, M=21yr (lognormal)
- New values based on DHS eSwatini (median ~17-18yr F, ~18-19yr M)

### Why
- Package defaults were generic, not country-specific
- EMOD calibrated to ~16.3yr F, ~17.5yr M
- Late debut (20yr) delayed network entry by 3-4 years, reducing HIV exposure for young adults
- Model was under-predicting prevalence in young women (15-25)

### Results
- **Incidence**: Higher overall, especially early epidemic. Peak ~5% vs ~4% before. Female incidence now clearly above male (more realistic sex differential). Model now overshoots PHIA 2011 incidence.
- **Prevalence by age/sex**: Young women (15-25) prevalence now higher and closer to SDHS/PHIA survey data. Overall prevalence somewhat elevated across age groups.
- **Implication**: Model is running "hotter" — needs re-calibration to let beta_m2f and other free parameters adjust to compensate for earlier debut.

### Figures
- `experiments/003_debut_age/before_*.png` — pre-change (from experiment 002)
- `experiments/003_debut_age/after_*.png` — post-change

### Decision
Keep change — direction is correct and produces more realistic age patterns. Re-calibration needed before further parameter tuning.
