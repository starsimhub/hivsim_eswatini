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

- Each calibration experiment gets `README.md` (question + plan) before
  any code runs, and `SUMMARY.md` (findings + decisions) before opening
  the next.
- Hand off to `project-workflow` for experiment scaffolding details.
- Hand off to `method-selection` after the coverage check (experiment 02).
- Hand off to `idm-azure` for VM lifecycle (already invoked).
