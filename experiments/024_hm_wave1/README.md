# Exp 024 — History matching, wave 1

> First calibration wave. Seven parameters from
> [023](../023_parameter_engineering/SUMMARY.md), targets settled there, run
> through the IDM `historymatching` package.

## Question

**Does the data rule out any of the parameter space, and does the emulator work
well enough to trust the cut?**

Wave 1 is not trying to produce a posterior. It removes regions of the box that
are inconsistent with a single macro-scale observation, checks the emulator
fits, and reports how much of the space survives. NROY ("not ruled out yet") is
a region, not a distribution — the Bayesian step comes after the waves converge.

## Parameter box

History matching takes **bounds only — no prior density.** Every point in the
box starts equally plausible. The package does not transform internally, so the
four parameters 023 sampled on a log scale are bounded in log space here and
un-transformed inside the simulator.

| bound | range | note |
|---|---|---|
| `log_beta_m2f` | ln(0.0096) – ln(0.025) | **floor raised from 0.008**, per 023: all three of its dead draws sat in 0.0082–0.0096, so this costs no coverage |
| `log_rel_beta_f2m` | ln(0.15) – ln(0.60) | |
| `log_s_f_young` | ln(0.8) – ln(3.0) | extends below 1.0 so the data can reject the mechanism |
| `age_gap_shift` | −2 to +3 years | linear |
| `log_age_gap_sd_mult` | ln(0.6) – ln(1.8) | |
| `prop_f0` | 0.45 – 0.85 | linear |
| `prop_m0` | 0.40 – 0.80 | linear |

`rel_init_prev` fixed at 0.2 and `conc_mult` at 1.0 — dropped in 023 on measured
evidence (max |ρ| of 0.10 and 0.19, and `conc_mult` confounded with
`rel_beta_f2m` at r = 0.81).

## Observation uncertainties — the part that decides whether this works

HM needs `(mean, std)` per observable. **There are no weights in history
matching, only standard deviations** — so "down-weight the deaths target"
translates into "give deaths a large σ", not into a weight coefficient. Getting
these wrong is the most likely way for wave 1 to fail, in either direction.

**σ from the surveys is small.** Derived from the published 95% CIs as
(ub − lb)/3.92, PHIA prevalence 15–49 by sex gives σ = 0.0043–0.0066.

**Replicate noise is 4–8× larger.** At one replicate and N = 10 000,
[020](../020_model_sizing/SUMMARY.md) measured a CV of 4.4% at high transmission
and 12.9% at low, so σ_rep ≈ 0.02–0.04 on a prevalence of 0.2–0.39.

**But replicate noise does not belong in σ_obs.** The emulator is fitted to
noisy simulator output, so its residual variance already absorbs simulation
noise; adding it to the observation σ would count it twice. Implausibility is

```
|z| = |emulator mean − observation| / sqrt(Var_emulator + σ_obs² + σ_disc²)
```

and the emulator term is handled by the package.

**So the extra allowance goes into model discrepancy, declared explicitly.**
Using σ_obs alone would put the model's known structural bias at 6–8σ and
collapse the NROY in wave 1 — the failure mode the skill warns about, and one
that would tell us nothing we did not already know.

We have unusually good evidence for what the discrepancy is. Seven experiments
(016–022) established a prevalence bias of about **−0.043 on the 15–49
aggregate** that no mechanism tested could remove. Setting

**σ_disc = 0.02 on prevalence**

places the best configuration measured to date (022 arm A) at roughly **2σ** —
not comfortable, but not automatically excluded. That is the anchor, and it is a
modelling judgement rather than a measurement. It is recorded here so that if
wave 1 collapses anyway, the first thing to revisit is this number and not the
model.

| target family | σ | source |
|---|---|---|
| prevalence 15–49 by sex | sqrt(σ_CI² + 0.02²) | PHIA CIs + declared discrepancy |
| prevalence ratios (F:M, young:old) | 0.15 relative | propagated, generous |
| incidence by sex and band | published CI, floored at 0.25 | SHIMS1 cohort / SHIMS2 Table 5.3.B |
| F:M incidence ratio | 0.35 | on a value near 2.0 |
| **peak AIDS deaths** | **2000** | **the down-weighting.** Model reaches 7 069 against UNAIDS 11 000, so σ = 2000 puts that gap at ~2σ. With a survey-like σ it would be 8σ and rule out everything |

**Implausibility threshold 4.0**, not the default 3.0 — the skill prescribes
raising it when the model has known structural misspecification, which is
precisely our situation.

## Feature selection

The emulator fits **one output at a time**, and the choice of which is the most
consequential decision in a wave. Wave 1 must use something macro-scale that a
coarse design can actually fit; fine temporal or age structure cannot be
emulated on 1000 points.

**Wave 1 emulates `prev_15_49_all_mean`** — overall adult prevalence 15–49,
averaged across the three PHIA years. Chosen manually rather than
automatically, per the skill's advice when the macro feature is obvious.

Everything else is *registered* but not emulated this wave. History matching is
a layered consistency check, and there is no double-counting, so later waves can
emulate different summaries of the same data:

| tier | wave | features |
|---|---|---|
| A | 1 | overall prevalence 15–49; peak AIDS deaths |
| B | 2–3 | F:M prevalence ratio, female young:old ratio, incidence by sex, F:M incidence ratio |
| C | 4+ | prevalence by age band × sex × year (24), if NROY is still large |

## Design

| | |
|---|---|
| points | 1000, Latin hypercube |
| replicates | **1 per point** — HM absorbs noise in the emulator; the 10 replicates 020 specified are for the coverage check and the later trajectory-selection step |
| N agents | 10 000 |
| emulator | `bayes_linear` (pure NumPy/SciPy) |
| waves this experiment | 1 |
| wall time | ~4.5 h on 8 workers |

**N = 10 000 rather than 20 000** deliberately. 020 recommended 20 000 because
two PHIA strata fall below 5 expected infected agents at 10 000 — but that only
matters for age-stratified targets, which are tier C. Wave 1 emulates an
aggregate, where 10 000 is ample. Later waves move to 20 000.

## Environment note

`historymatching` 2.0.1 requires TensorFlow, which has no wheels for Python
3.14. [`hm_shim.py`](../../hm_shim.py) makes the package importable by stubbing
`gpflow` (used only inside the GPR emulator's method bodies), so `bayes_linear`
works fully and `gpr` is unavailable. Two doc discrepancies worth recording: the
module is `historymatching`, not `history_matching`, and 2.0.1 replaced
`HistoryMatchingBuilder` with a single `HistoryMatching(...)` constructor.

## Success criteria

- **Clean:** emulator R² > 0.8 on the wave-1 feature, NROY between roughly 10%
  and 80% of the box, and the surviving region is interpretable. Wave 2 opens.
- **Uninformative:** NROY stays near 100%. The feature is not sensitive enough,
  or σ is too generous. Revisit the discrepancy allowance downward, or pick a
  sharper feature.
- **Collapse:** NROY near 0%. Either σ_disc = 0.02 is too small, or the model
  genuinely cannot reach the data — which is what 009 and 014 both found and
  what this experiment is partly re-testing with a better parameter set.
- **Emulator failure:** R² < 0.8. Try a different feature in the same wave
  before committing it (`drop_emulator_from_pending`).

## Accepted risk, recorded

**Re-identification (workflow step 7) has not been done.** Recovering known
parameters from synthetic data before fitting real data is the standard check
that the pipeline can recover anything at all, and 023 gave a specific reason to
want it: `beta_m2f` and `s_f_young` have effect signatures correlated at
r = 0.82, so a synthetic test would likely show them as mutually
unidentifiable. The researcher elected to proceed to calibration instead
(2026-09-03). Recorded so that if wave 1 returns an NROY that looks
suspiciously wide in those two directions, the cause is known in advance.

## Not in scope

- The Bayesian step. NROY is a region; converting it to a posterior via
  trajectory selection or MCMC comes after the waves converge, and needs a
  pseudo-likelihood that is deliberately not designed yet.
- The 2021 incidence hold-out.
- Tier B and C features — waves 2+.
