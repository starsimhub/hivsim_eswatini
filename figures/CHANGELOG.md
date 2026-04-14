# Figures Changelog

## 2026-04-07 -- Stratified ART coverage by age/sex (single sim)
- **Change**: ART input switched from `data/n_art.csv` (national total counts) to `data/art_coverage.csv` (proportion on ART by age, sex, and year from PHIA surveys)
- **Code modified**: `interventions.py` -- `make_interventions()` now reads `art_coverage.csv` and passes it to `sti.ART(coverage=art_data)` instead of using `n_art.csv` + `future_coverage`
- **Figures regenerated**: Fig 1 only (single sim, seed=1). Fig 2 requires `raw_data/` (not available locally). Fig 3 (network structure) unaffected.
- **Observations**:
  - Panel C (People on ART): lower ART numbers, especially post-2015 -- model now targets within-stratum proportions rather than national totals
  - Panel E (Deaths): somewhat higher, consistent with less aggressive ART scale-up
  - Panels A, B, D, F: broadly similar to baseline
- **Caveat**: Baseline was multi-sim (200 runs, median + CI bands); new result is single sim (seed=1). Stochastic variation contributes to differences. Full re-calibration needed for proper comparison.
- **Archive**: `archive/2026-04-06_before_stratified_art/` (baseline) vs `archive/2026-04-07_after_stratified_art/` (new)

## 2026-04-06 -- Baseline (pre-stratified ART)
- **ART input**: `data/n_art.csv` (unstratified national counts, 1990-2024) + `future_coverage={'year': 2022, 'prop': 0.90}`
- **Source**: Multi-sim calibration run (200 best parameter sets from 1000 Optuna trials)
- **Archive**: `archive/2026-04-06_before_stratified_art/`
