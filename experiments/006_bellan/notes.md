# Experiment 006 — Bellan 2015 acute HIV parameters

## Goal

Update the acute-phase HIV transmission parameters to Bellan 2015 central
estimates, replacing the prior values which produced an Excess Hazard Months
(EHM) of ~15 — substantially higher than Bellan's central estimate of 8.4.

## Changes

In `stisim/diseases/hiv.py`:
- `dur_acute`: `lognorm_ex(months(3), months(1))` → `lognorm_ex(months(1.7), months(1))`
- `rel_trans_acute`: `normal(loc=6, scale=0.5)` → `normal(loc=5.3, scale=0.5)`
- Resulting EHM ≈ (5.3-1)×1.7 ≈ 7.3 (vs Bellan's 8.4)

Cherry-picked from upstream `fix/396-bellan-acute-pars` (commit `a5e9ec1`).
Applied as `28adb68` on local main.

## Why

Per the upstream commit message: *"Previous overestimate was partly because
prior estimates failed to account for risk heterogeneity (Bellan 2015)."* When
high-risk individuals dominate transmission early in their infection, treating
all infected agents as equally infectious during acute (which is what the
classical RH-based estimates implied) overstates the population-level acute
contribution.

## Results (10-seed dashboard)

Reduces overall incidence by ~25-30% vs experiment 005. F:M ratio preserved.
Both incidence levels now further below PHIA targets — the Bellan correction
trades natural-history credibility for absolute incidence fit.

Re-calibration of `beta_m2f` (and ideally other free params, on a VM) is the
intended next step.

See `experiments/log.md` for the full results table and decision.
