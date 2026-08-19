# hivsim_eswatini — calibration project

## Intake (2026-05-07)

**Researcher.** Adam Akullian (Gates Foundation), solo, ~1 month timeline.

**Model.** Starsim/STIsim HIV model for Eswatini. ~10s per simulation.
Working state — recent fixes complete; ready for a real calibration. The
seven prior `experiments/001-007` were exploratory model-development, not
calibration.

**Question.** Calibration to support downstream decision analyses, most
prominently: *relative value of scaling up long-acting PrEP vs. improving
the treatment cascade to raise population viral suppression — and where
prevention still matters when suppression is high.* Eventual sensitivity
analyses will identify data gaps.

**Aim now.** Better-than-point-fit. Researcher is not committing to a
full posterior yet but wants to see whether characterising plausible
parameter space pays off. **Leading method (deferred to
`method-selection`):** history matching + trajectory selection.

**Calibration targets.**
- PHIA prevalence by age, sex, year (`calibration_data/prevalence_by_age_sex.csv` and PHIA files in `external_data/`)
- UNAIDS deaths (`data/eswatini_deaths.csv`)
- UNAIDS age distribution (`data/eswatini_age_1985.csv`)

**Validation hold-out (NOT a fitting target):**
- 2021 incidence (`calibration_data/incidence_2021_VALIDATION_ONLY.csv`)

**Explicitly excluded from fitting now:** ART coverage, VMMC, anything
else not listed above.

**Compute.** Local laptop now (10 sims at a time). Azure VM allocated,
setup in progress via `idm-azure` skill.

**Environment.**
- **Local (Windows laptop):** system Python 3.14 at `C:\Python314\python.exe`,
  with editable installs of starsim and stisim from
  `~/Dropbox/star_sim/{starsim,stisim}`. This is the existing working
  setup that has run experiments 001–007. Use it directly for any local
  Python work (`python experiments/NNN_*/run.py`). Corporate AV blocked
  the official `uv` Windows installer; revisit later via winget if
  needed.
- **Raccoon (Azure VM):** `uv` is installed at `~/.local/bin/uv` and
  managed via `pyproject.toml` at the repo root. On VM: `uv sync`
  resolves starsim/stisim from the `[tool.uv.sources]` editable paths
  *relative to the VM's checkout*, so the VM needs its own checkouts
  of starsim and stisim alongside `hivsim_eswatini` (or we adjust
  `pyproject.toml` to use PyPI versions on the VM). To revisit before
  experiment 009.

**Do not place a venv inside this directory** — it lives in
OneDrive/Dropbox and would cause sync churn.

**Repo state.** This is its own GitHub repo (`starsimhub/hivsim_eswatini`),
on branch `adam/model-dev`. No new branch for calibration — calibration
work continues the existing `experiments/NNN_*` numbering on the same
branch (current: `008_calibration_targets`).

## Structure

```
hivsim_eswatini/
├── pyproject.toml          # uv project root
├── (existing model code: run_sims.py, plot_*.py, ...)
├── calibration_data/       # fitting targets + validation hold-out
├── external_data/          # raw PHIA/UNAIDS files
├── data/                   # processed inputs (deaths, age, ASFR, ...)
└── experiments/
    ├── 001-007             # legacy model-development experiments
    └── 008_calibration_targets/  # first calibration experiment (current)
```

Calibration experiments continue the existing `experiments/NNN_*`
numbering. The transition from model-dev (007) to calibration (008) is
deliberate and recorded here, not in the directory layout.

## Workflow notes

- Each calibration experiment gets `README.md` (question + plan) before
  any code runs, and `SUMMARY.md` (findings + decisions) before opening
  the next.
- Hand off to `project-workflow` for experiment scaffolding details.
- Hand off to `method-selection` after the coverage check (experiment 02).
- Hand off to `idm-azure` for VM lifecycle (already invoked).
