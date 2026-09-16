# Exp 028 — The gate passes and the marginals do not move: three features and a halved σ cut the joint NROY by half while leaving every parameter at its prior

**Date:** 2026-09-16. **Model:** model-v1.5 (starsim 3.5.2 / stisim 1.5.11).
**Compute:** raccoon, 120 workers, 43 min, `rc=0` on the first attempt.
**Outcome:** **the pipeline ran to completion for the first time.** 025 and 027
both died at wave 2 of the pre-flight; the non-negative CD4 invariant fixed that.

**Questions.** (1) Can the pipeline separate `beta_m2f` from `s_f_young`?
(2) Does incidence constrain what prevalence cannot? (3) Does the NROY tighten
when σ_disc is halved, or collapse?

**Result.** The re-identification gate **passes** — truth inside the NROY on
7/7 parameters, with both the product `beta_m2f × s_f_young` and the **ratio**
recovered. On the pre-registered reading that is the good branch: the marginals
in those two directions are interpretable, not merely jointly constrained.

**And it does not matter as much as it sounds, which is the actual finding.**

| parameter | 024 width | 028 width | change |
|---|---|---|---|
| `log_beta_m2f` | 0.8900 | 0.8790 | −1.24% |
| `log_rel_beta_f2m` | 1.3026 | 1.2880 | −1.12% |
| **`log_s_f_young`** | 1.2563 | 1.2560 | **−0.02%** |
| `age_gap_shift` | 4.7310 | 4.7390 | +0.17% |
| `log_age_gap_sd_mult` | 1.0353 | 1.0339 | −0.14% |
| `prop_f0` | 0.3800 | 0.3784 | −0.42% |
| `prop_m0` | 0.3778 | 0.3760 | −0.48% |
| **joint NROY volume** | **87.184%** | **45.208%** | **−48.1%** |

024 was **one wave, one feature, σ_disc = 0.02**. 028 is **three waves, three
features, σ_disc = 0.01**. Adding incidence and the female young:old ratio,
halving the discrepancy allowance, and running two further waves bought a **48%
cut in joint volume and no marginal narrower than 1.3%.** `s_f_young` — the
parameter the entire tier-B feature design existed to constrain — moved by
**0.02%**.

The cut is happening entirely in correlated directions. Marginally, every
parameter is still at its prior.

## Observations

**obs 1 — "Truth inside the NROY" is a weak test at these widths.**
`frac_of_prior_width` is **0.90–0.95 on all seven parameters**. An interval
spanning 95% of the prior contains the truth close to by construction. The gate
therefore establishes that the pipeline is not grossly broken — prior, features,
emulator and implausibility are wired up correctly — and very little beyond that.
It should not be cited as evidence the features are sharp.

**obs 2 — The ratio is recovered, but it is still the worse-constrained
direction.** Both degeneracy quantities contain truth, and the
`beta_m2f`/`s_f_young` correlation inside the NROY is **−0.34** (against −0.59 in
025's n=25 smoke). But the ratio's interval is wider than the product's
(sd **0.52** vs **0.38**), which is the 023/024 degeneracy still present — merely
not severe enough to exclude the truth. The pre-registered pass is real; the
underlying confounding is not resolved.

**obs 3 — Answering question 3: the NROY tightens in volume and stalls, and the
stall is early.** Fraction of the box by wave —

| wave | pre-flight (synthetic) | wave 2 (real data) |
|---|---|---|
| 1 | ~50% | 50.3% |
| 2 | 39.8% | 46.3% |
| 3 | 39.0% | 45.2% |

Wave 1 does nearly all the cutting; wave 2 adds ~4 points on real data, wave 3
adds ~1. **The three features saturate after two waves.** Neither the collapse
025 flagged as a live risk (024's scan predicted 8.4% surviving on prevalence
alone at σ_disc = 0.01) nor the uninformative branch (>70%) occurred.

**obs 4 — Answering question 2: incidence did not constrain what prevalence
could not.** That was the central hypothesis for adding `inc_all_mean` — that a
rate would bite where a stock confounds rate with mortality, ART and cohort
history. Against 024's prevalence-only wave, the marginals are unchanged. Either
the incidence aggregate carries little information beyond prevalence at this σ
(the 10% relative discrepancy was added on 2026-09-07 and is a candidate), or
adult-aggregate incidence is simply too coarse a summary to separate these
parameters.

**obs 5 — Real data cuts less than synthetic.** 45.2% against the pre-flight's
39.0%, with near-identical marginals. Expected in direction, since the real-data
σ carries discrepancy terms the synthetic target does not, but the closeness of
the two sets of marginals is itself a signal that the constraint is dominated by
the prior and σ rather than by the data.

**obs 6 — σ_disc 0.02 → 0.01 bought nothing marginally.** 025's README argued the
halving was worth "roughly 5× the cutting power" on the grounds that 0.02 had
been calibrated against a bias 024 showed did not exist. The joint volume did
fall; no marginal did. The argument was sound and the effect is confined to
directions that do not show up in marginals.

**obs 7 — The CD4 crash is fixed, and the cause is not what 027 hypothesised.**
The invariant held: 43 minutes, `rc=0`, no retries. Telemetry across all 3,538
points:

```
points that tripped : 3 / 3538
total clip events   : 12      total agents clipped: 12
worst cd4 seen      : -inf
state mix: latent 1.000 | on_art 0.000 | acute 0.000 | falling 0.000 | art_disc 0.000
```

**100% latent, 0% on-ART — so 027's `cd4_increase` / negative-`dur_art`
hypothesis is dead.** `step_state` assigns latent agents `cd4 = cd4_latent`
directly, so **`cd4_latent` itself is `-inf`** for one agent. What produces an
`-inf` draw from `ss.normal(loc=500, scale=50)` remains unidentified.

**The upstream gap is independent of that and is worth reporting on its own:**
`step_state` validates CD4 with `np.isnan(...)`, and **`np.isnan(-inf)` is
`False`** — an infinity passes the check untouched and is caught only later, by
the positional rate lookup, as an `IndexError`.

Three distinct worker processes tripped at identical timesteps (`ti=352–355`),
which points at a specific agent under starsim's per-agent common random numbers
rather than a rare stochastic event. Only 3 of 3,538 points saw it, so whether
that agent reaches latency off-ART evidently depends on the parameter draw. Not
verified in the RNG source.

**Impact: nil.** 12 clipped agent-timesteps across 3,538 runs of 20,000 agents,
each pinned to the highest-mortality bin.

## What this does not establish

**Emulator R² is unverified.** The pre-registered criterion was R² > 0.8 per
feature, with anything below dropped rather than carried. The engine logs
emulator *training* but no quality metric at INFO or DEBUG; the per-feature
metrics 025's README expected are presumably inside `outputs/hm/*/checkpoint.pkl`.
**Until that is read, the claim that three simultaneous features did not dilute
attribution is unsupported** — and that claim is what justified deviating from
the 1–2 feature guidance.

**No figures.** 028 produced none, including the standard prevalence-fit figure
the workflow requires of every experiment.

## Acceptance

Accepted. The pipeline works end to end and the gate passes, which is what 028
existed to establish. The scientific result — that the feature set cuts joint
volume without identifying parameters — is a real answer, not a failure.

## Next

**The wave-3 plan needs revisiting before it is run.** The plan was the F:M
prevalence and incidence ratios, to pin `rel_beta_f2m`. 028's evidence is that
adding features of this kind cuts joint volume and leaves marginals alone, so
running wave 3 as designed is likely to reproduce that outcome at the cost of
another 40 minutes and an experiment folder. The question to settle first is
**structural: can this target set identify these seven parameters at all, or is
the joint-but-not-marginal outcome the best available?** That is a
`parameter-engineering` question — 023 did the calibration-sensitivity half, not
the identifiability half.

**Read the emulator R² from the checkpoints first.** It is cheap, it closes the
one unverified criterion, and if any feature is under 0.8 it changes the reading
of everything above.

**For exp 026's decision analysis, the practical consequence is concrete.**
Sampling this NROY gives marginals ≈ the prior for every parameter, so
propagating it into the PrEP-versus-cascade contrast will produce wide intervals
that mostly reflect the prior rather than the data. That is still more honest
than 026's current seed-noise-only interval, but it should be described as such.

**Upstream, two separate items:** the acute-decline floor (patch drafted, PR not
opened) and the `np.isnan` check that lets `±inf` through — the second is a
one-line change (`np.isfinite`) and is independent of whatever produces the
`-inf`.

**Still outstanding:** `pyproject.toml` pins with `>=`, so the correct stack held
here only because it was set by hand on the VM.

## Artifacts

| | |
|---|---|
| `outputs/preflight_recovery.csv` | the gate: truth, NROY interval, `inside`, `frac_of_prior_width` per parameter |
| `outputs/preflight_degeneracy.csv` | product vs ratio recovery, and the `beta_m2f`/`s_f_young` correlation |
| `outputs/preflight_nroy_summary.txt` | synthetic-data NROY, 3 waves, 39.047% |
| `outputs/wave2_nroy_summary.txt` | real-data NROY, 3 waves, 45.208% |
| `outputs/cd4_guard_summary.txt` | 3/3538 points tripped, state mix of clipped agents |
| `outputs/ensemble.parquet` | 3,538 points consolidated, 148,596 rows × 140 cols |
| `outputs/hm/{preflight,wave2}/checkpoint.pkl` | engine checkpoints — **emulator R² is in here** |
| `outputs/e028.log` | full run log including CD4 guard lines |
