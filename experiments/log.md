# Model Experiment Log

Running record of changes to the Eswatini HIV model, their rationale, and observed effects.

## Structure

- One folder per experiment: `NNN_topic/`
- Figures live in `NNN_topic/figures/` (use `before/` + `after/` subfolders if comparing pre/post state). At minimum, every closed-out experiment ends with `dashboard_fit.png` and `dashboard_network.png` summarizing the post-experiment model state.
- This `log.md` is the single source of truth for what changed, why, and the decision
- The repo's top-level `figures/` directory is **scratch only** — overwritten freely by plot scripts. Promote keepers into `experiments/NNN_topic/figures/` when closing out an experiment.
- `plot_dashboard.py` and `plot_beta_sweep.py` accept `--outdir`, so final runs can write directly into the experiment folder, e.g.:
    ```
    python plot_dashboard.py --label final --outdir experiments/004_beta_m2f/figures
    python plot_beta_sweep.py --outdir experiments/004_beta_m2f/figures
    ```

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
- **Incidence**: Higher overall, especially early epidemic. Peak ~5% vs ~4% before. Female incidence now clearly above male (more realistic sex differential). Model still within PHIA 2011/2016 CIs but running hotter than before.
- **Prevalence by age/sex**: Young women (15-25) prevalence now higher and closer to SDHS/PHIA survey data. Overall prevalence somewhat elevated across age groups.
- **Implication**: Model is running somewhat hotter — re-calibration would let beta_m2f and other free parameters re-optimize, but current fit is still reasonable.

### Figures
- `experiments/003_debut_age/before_*.png` — pre-change (from experiment 002)
- `experiments/003_debut_age/after_*.png` — post-change

### Decision
Keep change — direction is correct and produces more realistic age patterns. Re-calibration needed before further parameter tuning.

---

## 003b — Debut age follow-up: tighten SD and sensitivity sweep
- **Date**: 2026-04-22
- **Branch**: master

### What changed
- Tightened debut SD from 2.5 → 1 in `run_sims.py` (F=17.5 ± 1, M=18.5 ± 1)
- Added new analysis scripts:
  - `plot_debut_check.py` — boxplot + KDE density comparison (old vs new params)
  - `plot_debut_sensitivity.py` — sweep female debut 13–20 (male = female + 1), 10 seeds each
  - `plot_debut_vs_targets.py` — 10-seed run at current settings vs PHIA/UNAIDS calibration targets
- Fixed a latent bug in `plot_debut_check.py`: starsim deep-copies analyzers during sim init, so reading data from the original analyzer reference after `sim.run()` gives empty state. Retrieve via `sim.analyzers[...]` instead.

### Why
- SD=2.5 produced a very wide debut distribution (IQR ~3.4 yr) that overlapped heavily with the "old" setting — obscuring the effect of the mean change
- Wanted to understand the dose-response of debut age on incidence
- Wanted to check whether the new debut (F=17.5, M=18.5) now hits sex-stratified calibration targets

### Results
- **Tighter SD (1.0)** produces visibly separated old/new distributions; 2021 IQR ~1.4 yr for both sexes
- **Sensitivity sweep** (8 debut ages × 10 seeds): clear monotonic gradient — each year earlier in debut age increases peak incidence by roughly 0.3–0.5 per 100 PY during the 1990s–2000s epidemic peak. Gradient is comparable for both sexes. Curves converge after ~2015 as ART scale-up dominates.
- **Target comparison** (10 seeds, current settings):
  - **Male incidence**: tracks 2016 (~0.5 sim vs 0.85 target) and 2021 (~0.25 sim vs 0.2 target) reasonably well
  - **Female incidence**: sim underestimates — 0.7 vs 1.7 target (2016); 0.4 vs 1.4 target (2021). The ~2× female-to-male ratio in targets is not reproduced
  - **Prevalence by age**: young women (15–25) still low vs targets; 30–50 age range tracks well; 50–65 fits within bands

### Implication
Debut age change alone is insufficient to close the F:M incidence gap. The remaining gap points to `beta_m2f` (male-to-female transmission probability) as the next parameter to revisit.

### Figures
- `experiments/003_debut_age/debut_age_check.png` — boxplot of debut ages at 2000/2010/2016/2021
- `experiments/003_debut_age/debut_age_density_2010.png` — KDE densities at year 2010
- `experiments/003_debut_age/incidence_debut_comparison.png` — 1-seed incidence old vs new
- `experiments/003_debut_age/incidence_debut_sensitivity.png` — 10-seed sweep across debut ages
- `experiments/003_debut_age/incidence_vs_targets.png` — sim mean/band vs PHIA targets
- `experiments/003_debut_age/prevalence_by_age_sex_vs_targets.png` — age-stratified prevalence vs targets
- Archive mirror: `figures/archive/2026-04-22_debut_age_sensitivity/`

### Decision
Keep SD=1. Proceed to experiment 004: revisit `beta_m2f` to close the sex-stratified incidence gap.

---

## 004 — Close the F:M incidence and absolute-level gaps

- **Date**: 2026-05-04
- **Branch**: master
- **Motivation**: After 003b, female incidence under-predicts by ~2× and the F:M ratio is too narrow (~1.4× sim vs ~2× PHIA). Several candidate levers: `rel_beta_f2m`, `beta_m2f`, condom use, age-stratified susceptibility.

### Three sub-experiments

**(a) `rel_beta_f2m` sweep** — `plot_beta_sweep.py` over [0.125, 0.25, 0.5], 3 seeds each.
- Lower `rel_beta_f2m` widens F:M ratio but *also drops absolute incidence* (less F→M transmission means fewer infected men, fewer onward infections to women)
- 0.5 (default): F:M ~1.4×; 0.25: F:M ~2× (matches PHIA); 0.125: F:M ~3× (overshoots)
- **Conclusion**: `rel_beta_f2m=0.25` reproduces the F:M ratio but does not close the absolute-incidence gap on its own.

**(b) Condom-use sensitivity sweep** — `plot_condom_sweep.py` scaling stisim defaults by [1.0, 0.7, 0.5, 0.3], holding `rel_beta_f2m=0.25`. 3 seeds × 4 scenarios.
- Stisim defaults (post-2005 plateau ~0.7 mixed, 0.9 MM/HH) appear to over-protect when interpreted as act-level usage (DHS surveys ask "ever used at last sex" — consistency is lower than ever-use)
- **0.5× scaling** brings 2016 sim to F~1.4 / M~0.8 vs PHIA F=1.7 / M=0.85 — closest fit
- 0.3× overshoots female slightly
- **Decision**: Adopt 0.5× scaling for non-LL pairings. LL kept at marital baseline (~1%). FSW–Client also scaled (could argue for holding 0.95 — left for future calibration).

**(c) Adopt `rel_sus_age` for young women** — cherry-picked from unmerged stisim PR `fix/395-age-dep-rel-beta`.
- New API: `rel_sus_age=[(15,25,'f',1.7), (25,50,'f',1.0), (15,50,'m',1.0)]`
- 1.7× susceptibility for women 15–24 (per-agent multiplier on `rel_sus`)
- Motivated by mucosal/cervical immaturity literature; partly compensates for the model's tendency to under-predict young-women prevalence
- Combined with `rel_beta_f2m=0.25`, young women carry ~6.8× per-act acquisition risk vs same-age men in this configuration

### Stisim provenance

- Cherry-picked commits onto local `main` (preserving the local ART allocation fix from experiment 002):
  - `827d720` — Add age-dependent M→F susceptibility: `rel_beta_m2f_by_age`
  - `d18413e` — Rename `rel_beta_m2f_by_age` → `rel_sus_age`; add sex stratification
- Skipped the v1.5.4 version-bump commit. Did NOT pull from upstream `main` because subsequent upstream changes (PMTCT uplift, coverage refactor) would have wiped the local ART fix.

### Final configuration (committed in `run_sims.py`)

```python
hiv = sti.HIV(
    beta_m2f=0.01,
    rel_beta_f2m=0.25,                          # was 0.5 (default)
    rel_sus_age=[(15,25,'f',1.7), (25,50,'f',1.0), (15,50,'m',1.0)],  # new
    eff_condom=0.85,
    ...
)
condom_data = csv * 0.5 (non-LL rows only)      # was csv unchanged
```

### Results (10-seed dashboards)

- **Female incidence (15–49)**: 2016 sim ~1.5 vs PHIA 1.7 — within band; 2021 sim ~0.8 vs PHIA 1.4 — below
- **Male incidence (15–49)**: 2016 sim ~0.5 vs PHIA 0.85 — below; 2021 sim ~0.3 vs PHIA 0.2 — slightly above
- **F:M incidence ratio**: ~2× across time, matches PHIA
- **Prevalence by age & sex**: Now tracks all four survey years (2007/2011/2016/2021) tightly across age bins. Largest residual gap: female prevalence 50+ slightly low.
- **ART coverage by age & sex**: Matches PHIA across age groups for both sexes.

### Caveats

- Original calibration found `beta_m2f` and `eff_condom` against the *unscaled* stisim default `condom_data`. Halving the condom data without re-calibrating means those two parameters are now miscalibrated — `beta_m2f=0.01` is likely too low, `eff_condom=0.85` is conditional on the old condom_data. **Re-calibration with a `condom_scale` free parameter is the right next move** before further structural changes.
- `rel_sus_age=1.7` for young women is the dev's example value — not literature-calibrated. Should also be a candidate free parameter.

### Figures (in `experiments/004_beta_m2f/figures/`)

- `beta_f2m_sweep.png` — `rel_beta_f2m` sweep (sub-experiment a)
- `condom_use_sweep.png` — condom scaling sweep (sub-experiment b)
- `dashboard_fit_004_final.png` — 10-seed fit dashboard at adopted configuration
- `dashboard_network_004_final.png` — 10-seed network dashboard at adopted configuration

### Decision

Adopt all three changes as the new baseline (committed to `run_sims.py`). Defer re-calibration to a separate cycle on a VM. Proceed to experiment 005.

---

## 005 — VMMC integration with prevalence-target semantics

- **Date**: 2026-05-06
- **Branch**: master (HIVsim) / local main (stisim)
- **Motivation**: Add Voluntary Medical Male Circumcision to the model with PHIA-aligned prevalence targets, since the previous panel F was a placeholder. Bellan acute params shifted to experiment 006.

### What changed

**Stisim patch** (local mod to `star_sim/stisim/stisim/interventions/hiv_interventions.py`):
- Added `prevalence_target=True` flag to `VMMC` (default on for stratified data)
- New `_apply_prevalence_target()` method: for each (age_bin, sex) stratum, top up circumcisions to match `p_vmmc` interpreted as target prevalence (proportion of all alive males in stratum). Never reduces (VMMC is irreversible).
- Mirrors ART's `art_coverage_correction()` design pattern
- Legacy hazard interpretation preserved if `prevalence_target=False` or coverage is non-stratified
- Candidate for upstream PR (alongside the experiment-002 ART fix)

**Eswatini repo**:
- `data/vmmc_coverage.csv` — 5-year age bins ([10,15) through [60,65)), male only, at 2007/2016/2021
- `interventions.py` — VMMC enabled with the new CSV
- `plot_dashboard.py` — added `VMMCPrevByAge` analyzer; replaced panel F placeholder with real VMMC coverage plot (4 broader bins: [15,25), [25,35), [35,45), [45,65))

### Data triangulation

See `experiments/005_vmmc/notes.md` for full detail. Summary:
- **2007 values**: from EMOD JSON, **validated against actual SDHS 2007 data** (15-19 through 45-49 match within rounding; the [35,39)=0.197 ≈ 0.20 anomaly is a real cohort effect, not a placeholder). 50+ values remain EMOD extrapolations since SDHS 2007 didn't survey men 50+.
- **2016 values**: from EMOD JSON (SHIMS2 — strongly age-stratified, well-attested)
- **2021 values**: replaced with **SHIMS3 official Table 12.5** (medical + nonmedical circumcision prevalence by 5-year bin), found in `data/241123_SHIMS_ENG_RR3_Final-1.pdf`. Cross-checks: SHIMS3 reports 47.2% medical + 1.1% nonmedical = 48.3% total for 15-49.
- The EMOD 2021 values were "EXTRAPOLATED" per the JSON comment. They turned out
  to be within ±5pp of SHIMS3 actuals for ages 15-39, but underestimated 45-59
  (~0.18-0.20 EMOD vs 0.24-0.27 SHIMS3 actual). My initial cohort-progression
  fix to [40,45)=0.30 was nearly right (SHIMS3 actual 0.325).

### Stisim provenance (cumulative)

Local main now has 3 commits ahead of upstream `main` (cherry-picked) plus 1 local-only mod:
- `827d720` — Add age-dependent M→F susceptibility (cherry-pick from PR #395)
- `d18413e` — Rename rel_beta_m2f_by_age → rel_sus_age (cherry-pick)
- **NEW (uncommitted)**: VMMC prevalence-target patch in `hiv_interventions.py`
- **Still uncommitted**: ART allocation fix from experiment 002 (also in `hiv_interventions.py`)

Both local mods are candidates for upstream PRs. Should be committed to a feature branch in stisim before VM workflow.

### Results (10-seed dashboards)

**VMMC panel (F)**:
- 15-25 bin: sim tracks the 2016/2021 targets cleanly (rapid scale-up captured)
- 25-35 bin: tracks well; modest overshoot in 2021 due to cohort aging from [15,25)
- 35-45 bin: tracks 2016 well; overshoots 2021 target (cohort effect from younger bin)
- 45-65 bin: similar pattern; overshoots due to inflow from [35,45)

The cohort-aging overshoots are *biologically correct* — VMMC is irreversible, so high prevalence in younger bins propagates upward over time. The PHIA point estimates are cross-sectional snapshots and don't fully reflect this stock-vs-flow distinction.

**Male incidence (panel A)**: drops further from 004 (sim 2016: ~0.3 vs PHIA target 0.85). Expected — VMMC adds a 60% acquisition reduction for ~30%+ of adult males by 2016. Closing the gap requires re-calibration (especially `beta_m2f` upward).

**Female incidence (panel A)**: roughly unchanged from 004 (~1.2 vs target 1.7). VMMC affects men's susceptibility, indirect effect on women via fewer infected men.

**Prevalence by age & sex (B, C)**: unchanged fit quality from 004 — strong tracking across all 4 survey years.

**ART coverage (D, E)**: unchanged from 004.

### Caveat on display

The network dashboard's panel C (condom use by partnership) still shows the *unscaled* CSV values, while the sim runs with the 0.5× scaled version (per experiment 004). Display-only issue — the sim itself uses the scaled values. Fix: have the plot apply the same scaling, or have run_sims.py write the scaled CSV to disk. Deferred.

### Caveats on the data

- 2007 ages 40+ are flat at 0.12 in the EMOD JSON — likely placeholder, real SDHS may differ
- 2007 [35,40) at 0.20 is anomalously high vs surrounding bins; could reflect traditional MC cohort effect or data error
- 2021 values are projections, not direct SHIMS3 measurements. Population-weighted average for [15,49) under our values ≈ 41%; SHIMS3 published overall MMC may be higher (~50–60% per recollection). If SHIMS3 reports a higher aggregate, scale all 2021 values up proportionally.

### Figures (in `experiments/005_vmmc/figures/`)

- `dashboard_fit_005_final.png` — 10-seed fit dashboard with new VMMC panel
- `dashboard_network_005_final.png` — 10-seed network dashboard (unchanged from 004 except cohort effects)

### Decision

Adopt VMMC integration as the new baseline. Male incidence under-prediction is a known consequence of adding the VMMC effect without re-calibrating `beta_m2f` — the right place to address it is the VM-based recalibration cycle, with `condom_scale`, `eff_circ`, `rel_sus_age` and the existing free parameters jointly tuned. Proceed to experiment 006: cherry-pick Bellan acute params (`a5e9ec1`).

### Followup (same-day): pre-program baseline correction

The initial 005 commit had only 2007/2016/2021 data points in the CSV, which caused `sc.smoothinterp` to propagate the SDHS 2007 [35,40)=0.20 value back to the sim start (1985). Pre-2000 [35,40) MMC was therefore an unrealistic 20%.

Fix: added a **1990 anchor row** to `data/vmmc_coverage.csv` with values 0.005–0.04 calibrated to SHIMS3 nonmedical-MMC rates. Pre-2007 trajectory now ramps linearly from low traditional-MC baseline up to the validated SDHS 2007 values. User confirmed SDHS 2007 measurements directly, validating the EMOD JSON for 15-49 within rounding (including the 35-39 bump).
