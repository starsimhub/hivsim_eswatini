# Experiment 009 — Prior predictive coverage check

## Question

Can the model produce the observed data, anywhere in its prior? Before
spending compute on calibration, draw 50 parameter sets from the
implicit Optuna prior, simulate one trajectory per draw, and check
whether the frozen calibration targets (experiment 008) fall inside the
ensemble.

This is the cheapest sanity check in a calibration pipeline. A pass means
the model + prior + observation model are *compatible* with the data.
A fail tells us exactly which of three things to fix: prior too narrow,
model can't reach the data, or wrong observation model. Skipping it is
how people lose a week chasing convergence on a structurally impossible
fit.

## Method

**Quick coverage** (per `coverage-check` skill): 50 draws × 1 replicate.
The question is binary — *can* the model reach the data — and one rep
per draw is enough to answer it. (Replicate-count is a `model-setup`
question we'll address separately if needed.)

### Prior

Six uniform priors, lifted directly from `run_calibrations.py` — these
are the parameters the previous Optuna setup was searching over:

| Parameter | Module | Low | High |
|---|---|---|---|
| `beta_m2f` | hiv | 0.002 | 0.014 |
| `eff_condom` | hiv | 0.5 | 0.9 |
| `rel_dur_on_art` | hiv | 1.0 | 20.0 |
| `prop_f0` | structuredsexual | 0.55 | 0.9 |
| `prop_m0` | structuredsexual | 0.55 | 0.8 |
| `m1_conc` | structuredsexual | 0.05 | 0.3 |

Drawn with a fixed seed (LHS-style uniform, one seed) so the experiment
is reproducible.

### Simulation

`make_sim(seed, hiv_pars=..., network_pars=..., stop=2026)`. Default
`n_agents=10e3`. One replicate per draw (seed = draw index). Parallelized
with `sc.parallelize` across local cores; each worker pinned to a single
BLAS thread (`OMP_NUM_THREADS=1` etc., per the existing `run_calibrations.py`).

### Outputs

- `outputs/ensemble.parquet` — all 50 sims' yearly results, columns trimmed
  to what's needed to compare against targets (prevalence by age×sex×year,
  HIV deaths over time)
- `outputs/draws.csv` — the prior parameter draws indexed by `par_idx`
- `outputs/coverage_summary.csv` — per-target row: observation, simulated
  5–95 % range, in/out flag
- `outputs/coverage_prevalence.png` — ensemble vs PHIA prevalence (3 PHIA
  years × 2 sexes panels)
- `outputs/coverage_deaths.png` — ensemble vs UNAIDS deaths over time

### Reading the result

- **Fully covered** (all PHIA + deaths inside 5–95 % envelope): proceed
  to `method-selection` and re-identification.
- **Mostly covered, edge points outside**: usually fine; note in SUMMARY.
- **Systematic miss**: stop. Diagnose which of the three (prior too narrow,
  model structure, observation model) is the problem before any fitting.

## Expected runtime

50 sims × ~10 s/sim = ~500 s sequential. With 10 parallel workers
(your laptop's typical capacity per intake), wall time should be
~1–2 min plus import/startup overhead.

## Out of scope

- No likelihood. No fitting. No method choice. We are confirming the
  ensemble *covers* the data, not finding the best fit.
- No sensitivity to `n_agents`, replicate count, or run-time — that's
  `model-setup` if a follow-up is needed.
- No 2021 incidence (validation hold-out) and no age distribution
  (deferred per experiment 008 SUMMARY).

## Status

`README.md` written. `run.py` to follow.
