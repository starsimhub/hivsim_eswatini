# SUMMARY — experiment 009 (prior predictive coverage check)

**Status:** complete. **Coverage failed** — diagnosis required before
any fitting.
**Date:** 2026-05-07.
**Compute:** raccoon (120 cores). 50 sims drawn from prior, ~38 s wall.

## Top-level result

**30/89 target rows inside the 5–95 % envelope (34 %).**

| Quantity | Inside | Total | % |
|---|---|---|---|
| AIDS deaths | 7 | 35 | 20 % |
| PHIA prevalence | 23 | 54 | 43 % |

Both miss "systematically," not at edges. Per `coverage-check`: **stop
and diagnose** — do not start calibration on this prior.

## Diagnosis

The two failures are mutually consistent and point to a single class of
issue: **the model lets people accumulate as PLHIV without enough
mortality flow.**

- **Deaths**: simulation envelope is below UNAIDS at the peak (median
  ~6 k/year, p95 ~10 k vs observed ~11 k) and dramatically below in the
  post-2010 plateau (sim → near zero, observed → ~3 k/year stable). See
  [coverage_deaths.png](outputs/coverage_deaths.png).
- **Prevalence**: 2007 fit is reasonable. By 2011 and 2016 the sim
  envelope overshoots — Female 2016 sim p95 reaches 60–70 % vs observed
  40–50 %. See [coverage_prevalence.png](outputs/coverage_prevalence.png).

If deaths are too few, PLHIV accumulate → prevalence drifts up over
time → 2016 overshoot. The two findings are not separate bugs.

## Three concrete suspects

### 1. `rel_init_prev` is fixed, not in the prior

`run_sims.py` hard-codes `rel_init_prev=0.2` (was 0.1, doubled in
experiment 007). `run_calibrations.py` has `rel_init_prev` commented
out of the prior. If the seed prevalence is below where it needs to be
for Eswatini's epidemic, no other parameter can compensate.

**Fix to test:** add `rel_init_prev` back to the prior with a wide
range (e.g. 0.05 – 0.5 or wider).

### 2. `rel_dur_on_art` (1.0 – 20.0) may be too narrow

Modern ART roughly normalises life expectancy — equivalent to
`rel_dur_on_art ≫ 20` for most ages. The post-2010 plateau in observed
deaths is real; it represents non-AIDS mortality among PLHIV plus the
small residual from late-presenters. The model dropping to zero deaths
suggests ART durations are pushing survival too far.

But this is an upper-bound symptom; widening alone won't fix the peak
miss in 1995–2005 (pre-ART era).

**Fix to test:** widen `rel_dur_on_art` upper bound (or switch to
log-scale).

### 3. HIV mortality / progression pars aren't in the prior

The Bellan 2015 acute pars (experiment 006) and `rel_sus_age`
(experiment 004) are fixed in `run_sims.py`. The CD4-progression /
time-to-death distribution is governed by stisim defaults. If
Eswatini's real-world HIV mortality runs faster than these defaults
(which is plausible — Eswatini has historically had high HIV mortality
even in pre-ART era), the model will systematically underproduce deaths.

**Fix to test:** identify the HIV-mortality parameters in stisim and
expose at least one (e.g. an overall mortality multiplier) to the prior.

## What did NOT explain the miss

- Beta range looks adequate — peak prevalence in 2007 is roughly the
  right magnitude, so transmission rate isn't grossly off.
- Network structure parameters (`prop_f0`, `prop_m0`, `m1_conc`)
  primarily move *who* gets infected, not *how many*; their range is
  unlikely to drive the systematic miss.
- Observation model is a plausible secondary suspect for deaths
  (UNAIDS Spectrum estimates are model-derived too), but not for PHIA
  prevalence (direct survey measurement).

## Next experiment

**010 — Diagnostic prior expansion.** Re-run the coverage check with:
- `rel_init_prev` added to the prior (suggested 0.05 – 0.5)
- `rel_dur_on_art` upper bound raised (suggested 1 – 50)
- One HIV-mortality multiplier added — TBD after a quick read of the
  stisim HIV module to identify the right parameter

If coverage passes, the failure was prior-narrowness (a fix) and we
proceed to method-selection. If coverage still fails, the model itself
needs structural work before any calibration is meaningful.

50 sims again, on raccoon. Single replicate. Same target set frozen
in experiment 008.

## Compute footnote

Workflow established with raccoon for the first time:

- Local stisim has 3 commits ahead of `starsimhub/stisim:main` plus
  uncommitted edits. `tar | ssh` from local → raccoon worked but is
  ad-hoc; long-term the local stisim commits should be pushed to a
  branch so raccoon can `git pull`.
- `experiments/*/outputs/` is gitignored (correct for ensembles), so
  the calibration_targets.csv from experiment 008 needed to be `scp`-ed
  manually. Worth deciding whether the canonical targets file should
  be tracked in git.

Sim wall time on 120 cores: 38 s for 50 sims. Suggests we can comfortably
run a 1 000-sim coverage check in ~10 min if we want a more thorough
test.
