# HIVsim eSwatini

Agent-based HIV model for eSwatini, built on STIsim/Starsim. Calibrated
to UNAIDS and PHIA data. Code is maintained by the starsimhub team.

## Intake

**Model:** Starsim/STIsim HIV model, ~10 000 agents, 1985–2031, monthly
timesteps. ~5–15 sec/sim. Local editable installs of starsim and stisim
at `~/GIT/starsim` and `~/GIT/stisim`.

**Question:** Explore PrEP targeting for FSW in eSwatini. Requires a
*posterior* over plausible epidemic trajectories — not a single best-fit
point — to support intervention comparisons under uncertainty.

**Data:** UNAIDS estimates (prevalence, incidence, PLHIV, deaths, ART
counts) and PHIA household surveys (age/sex-specific prevalence,
2007/2011/2016/2021). Calibration CSV: `data/eswatini_hiv_calib.csv`.

**Constraints:** 120-core VM available for heavy runs; local development
on Mac (10 CPUs). Target: solid calibration within 1–2 weeks. Repo is
PUBLIC and managed by a separate team — all calibration work lives on
named feature branches and must not be pushed to `main`.

**Method:** Switching from Optuna (point fit) to ABC to obtain a
parameter posterior. Coverage check first to confirm the prior can
envelop the data.

**Current fit status:** Poor — epidemic collapses to near-zero incidence
by 2020 while data shows sustained ~1.7%/year. `rel_init_prev` is fixed
at 0.1 and excluded from calibration; suspected root cause.

**Environment:** `uv` with `pyproject.toml`. Run scripts with
`uv run python <script>`. starsim and stisim are local editable installs
at `~/GIT/starsim` and `~/GIT/stisim`.

**Experiments:** `experiments/01_coverage_check/` — prior predictive
check.
