# Changelog

All notable changes to this project will be documented in this file.
The format loosely follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project uses semantic versioning.

## [Unreleased]

### Added
- `pyproject.toml` with pinned dependencies (numpy, pandas, sciris, matplotlib, starsim, stisim).
- `tests/` directory with a smoke test and analyzer-invariant tests (pytest).
- `CONTRIBUTING.md` with dev setup and test instructions.
- `data/README.md` documenting input-data provenance (UNAIDS, PHIA, census, ART).
- Shared `run_one` / `stack` / `summarize` helpers in `utils.py` for multi-seed plotting.

### Changed
- All hardcoded relative paths in run/plot scripts now anchor to `Path(__file__).parent`,
  so scripts work regardless of working directory.
- `plot_calibration.py` no longer executes simulations at import time; the work is
  inside `main()` and only runs under `__main__`.
- `NetworkSnapshot._capture_snapshot` vectorizes its per-age-bin loop with `np.bincount`.
- PHIA incidence reference values in `plot_incidence.py` are now in a named, documented
  constant (`PHIA_INCIDENCE`) instead of magic numbers.

### Removed
- Commented-out VMMC intervention scaffolding in `interventions.py`.

## [0.1.0]

Initial commit: HIV transmission model for Eswatini with structured sexual network,
HIV testing (FSW/general/low-CD4), ART, PrEP, calibration to UNAIDS/PHIA, and
plotting scripts for prevalence, incidence, ART coverage, and network structure.
