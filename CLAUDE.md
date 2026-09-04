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
- UNAIDS deaths (`data/eswatini_hiv_calib.csv`, column `hiv.new_deaths`)
  — **corrected 2026-08-19.** This entry, and the source table in
  `experiments/008_calibration_targets/README.md`, previously cited
  `data/eswatini_deaths.csv`. That is a different file: all-cause
  mortality *rates* by age/sex fed to `ss.Deaths` as background
  mortality, not the UNAIDS AIDS-death *counts* used as a target.
  008's README is closed and left as written; see exp 016.
- UNAIDS age distribution (`data/eswatini_age_1985.csv`)
- **SHIMS incidence, adult aggregates only — added 2026-09-04.**
  `calibration_data/incidence_by_age_sex.csv`, rows with `fit=True`:
  2011 (SHIMS1 cohort, 18-49) and 2016 (SHIMS2 Table 5.3.B, 15-49), by
  sex. Four rows. Built by `incidence_construction.py`.

  This entry was missing until 2026-09-04. Incidence was added to the
  target set during exp 023/024 and the list above was never updated,
  so the file said incidence was excluded while the code registered it.
  Recorded now.

  **The six age-banded incidence rows are retained but NOT fitted**
  (`fit=False`). The published female profile is essentially flat
  (1.67, 1.54, 2.09 across 15-24, 25-34, 35-49), which is the expected
  signature of false-recency bias: a recency assay credits a fraction
  FRR of long-standing infections as recent, and against a susceptible
  denominator the spurious incidence scales as FRR x P/(1-P)/MDRI —
  5x larger in women 35-49 than 15-24, 18x in men. See
  `incidence_construction.py`'s `FIT_POLICY`. **Consequence: the
  model's incidence age profile is unconstrained by data.** Exp 024's
  `figures/incidence_age_profile.png` shows what it does anyway.

  The 2011 rows carry an extra sigma term (f 0.126, m 0.194) because
  SHIMS1 publishes 18-49 and `PopByAgeSex` works in 5-year bands, so
  the model can bracket that range but not express it.

**Validation hold-out (NOT a fitting target):**
- 2021 incidence (`calibration_data/incidence_2021_VALIDATION_ONLY.csv`)

**Explicitly excluded from fitting now:** ART coverage, VMMC, the
age-banded incidence rows above, anything else not listed here.

**When the target set changes, update this list in the same commit.**
The incidence omission above went unnoticed for two experiments.

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

**Repo state.** This is its own GitHub repo (`starsimhub/hivsim_eswatini`).
No separate branch for calibration — calibration work continues the
existing `experiments/NNN_*` numbering (current:
`016_double_counted_mortality`).
See **Collaboration and branching** below for the branch model, which
changed on 2026-08-19.

## Collaboration and branching (decided 2026-08-19)

**Collaborator.** Daniel Citron. Both of us work on the same model; the
plan is one integrated model, then calibration, then downstream
workflows built on the calibrated model.

**`main` is the integration point** — not a personal branch. Until
2026-08-19 the real work lived on `adam/model-dev` while `main` sat 23
commits stale, which meant anyone cloning the repo got a model with no
VMMC fix and no network fix. PR #5 merged that work back. Rationale for
`main` over a shared personal branch: a personal branch used as an
integration point takes on trunk responsibilities without trunk
protections, and PRs into it have a base that moves under the
contributor.

**Per-experiment branches.** Each experiment gets a short-lived branch
(`exp/016-model-setup`) PR'd into `main`. Experiment folders are
append-only and self-contained, so they rarely conflict; the branch
exists to satisfy branch protection and to give Daniel a review point.
If this becomes friction, committing experiments directly to `main` is
acceptable — the discipline that matters is the README/SUMMARY pair and
the model tag, not the branch.

**Tag the model before calibrating.** A calibration is only valid for
the model version it ran against. Tag the integrated model and record
the tag in each experiment's `config.yaml`, alongside the starsim and
stisim versions already recorded there. Model changes during a
calibration go to `main` for the next cycle; they do not disturb a run
in flight.

**Shared code surface** — the files where conflicts actually happen:
`run_sims.py`, `interventions.py`, `vmmc.py`, `hiv_mortality.py`.
Everything under `experiments/NNN_*/` is append-only.

**Open decision.** PR #2 (`update-coverage`, Robyn Stuart) adds an ANC
HIV testing intervention. It changes diagnosis rates, hence ART uptake,
hence prevalence and deaths — so it must land *before* the model is
frozen for calibration, or explicitly wait for the next cycle. It will
conflict in `interventions.py` now that #5 has merged. Two of its three
changes (the ART coverage-API migration, the VMMC line) are already
superseded by work on `main`; only the ANC testing product is new.

**CodeQL false positive.** `py/clear-text-logging-sensitive-data` fires
on `plot_debut_check.py` because a loop variable is named `sex`. It is a
hardcoded label on aggregate simulated statistics, not personal data.
Expect this rule to keep firing on modelling code that uses `sex`, `age`
or `race` as model dimensions.

## Structure

```
hivsim_eswatini/
├── pyproject.toml          # uv project root
├── (existing model code: run_sims.py, plot_*.py, ...)
├── calibration_data/       # fitting targets + validation hold-out
├── external_data/          # raw PHIA/UNAIDS files
├── data/                   # processed inputs (deaths, age, ASFR, ...)
├── vmmc.py                 # VMMCPrevalenceTarget — in-repo model fix (exp 015)
├── hiv_mortality.py        # HIVMortalityMultiplier — calibration knob (exp 014)
└── experiments/
    ├── 001-007             # legacy model-development experiments
    ├── 008-015             # calibration targets, coverage checks, model fixes
    └── 016_double_counted_mortality/  # current
```

Model fixes live in this repo as subclasses, never as edits to the
editable starsim/stisim checkouts — an upgrade silently overwrites those
(this is how the exp-005 VMMC patch was lost, see exp 015).

Calibration experiments continue the existing `experiments/NNN_*`
numbering. The transition from model-dev (007) to calibration (008) is
deliberate and recorded here, not in the directory layout.

## Workflow notes

- **Every experiment produces the standard prevalence-fit figure.** Call
  `standard_figures.plot_prevalence_fit(df, label, outdir/'prevalence_fit_vs_phia.png')`
  from `run.py`. It is one shared implementation, so the figure cannot
  drift between experiments — which is exactly what happened before:
  `plot_dashboard.py` produced `dashboard_fit_*.png` for 003–015, then
  broke at 016 because it runs its own sims and could not consume an
  A/B harness's output. 016, 017, 018 and 020 ended up with no
  prevalence-fit figure at all, and 019/021 grew a second incompatible
  format. `plot_fit_progression.py` retrofits the standard figure onto
  past experiments where outputs survive, and calls the same function.
- **Save what the figure needs.** The `analyzers.PopByAgeSex` analyzer
  must be in the analyzer list, and `n_infected_*` per band must survive
  the `KEEP` column filter. 016 and 017 kept `n_alive` per band but only
  *aggregate* infected, so their age-stratified fit is permanently
  unrecoverable.
- Each calibration experiment gets `README.md` (question + plan) before
  any code runs, and `SUMMARY.md` (findings + decisions) before opening
  the next.
- Hand off to `project-workflow` for experiment scaffolding details.
- Hand off to `method-selection` after the coverage check (experiment 02).
- Hand off to `idm-azure` for VM lifecycle (already invoked).
