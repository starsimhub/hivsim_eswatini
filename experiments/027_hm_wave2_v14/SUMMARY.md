# Exp 027 — Same crash, fix verified active: the acute floor was necessary and not sufficient

**Date:** 2026-09-10. **Model:** model-v1.4 (starsim 3.5.2 / stisim 1.5.11).
**Compute:** raccoon, 120 workers, 22 min before the crash.
**Outcome:** **failed run, and the failure is informative.** None of the three
questions answered — they carry forward to 028 unchanged.

**Question.** Unchanged from [025](../025_hm_wave2/SUMMARY.md): can the pipeline
separate `beta_m2f` from `s_f_young`; does incidence constrain what prevalence
cannot; does the NROY tighten at halved σ_disc or collapse?

**Result. The identical crash, at the identical place:**

```
stisim/diseases/hiv.py:540  step_state
stisim/diseases/hiv.py:401  make_p_hiv_death
IndexError: index 6 is out of bounds for axis 0 with size 6
```

Wave 1 of the pre-flight completed; wave 2 died in the simulation phase, same as
025. 025's fix did not prevent it.

## The finding: the fix was active, so a second defect exists

This is the whole value of the run, and it needed establishing rather than
assuming. The model-v1.4 fix **was** loaded:

| check | result |
|---|---|
| `5f8e367` (acute floor) ancestor of launched tree `2e29007` | yes |
| `040f547` (name pin) ancestor of `2e29007` | yes |
| `run_sims.py` at `2e29007` | `hiv = (hiv_class or HIVCD4Floor)(**hiv_kwargs)` |
| `hiv_cd4_floor.py` at `2e29007` | `np.maximum(cd4_end, ...)` present in `acute_decline` |
| run log header | `model=model-v1.4` |

So CD4 still reached `make_p_hiv_death` negative with the acute decline floored.
**The unfloored acute decline was a real defect — demonstrated, one extra step
overshooting the set-point by ~350 CD4 and the second crossing zero — but it is
not the only route to a negative CD4.**

Recorded plainly because the 025 SUMMARY reads as though the cause had been
found. It had not. It had found *a* cause.

## What is ruled out, and what is left

**Ruled out by inspection:** `falling_decline` and `post_art_decline` both floor
at `cd4_end = 1`; the latent phase assigns `cd4_latent` directly; `acute_decline`
is now floored.

**Leading hypothesis — `cd4_increase`, the on-ART recovery logistic:**

```
cd4_now = 2*g/(1 + exp(-dur_art*growth)) - g + cd4_preart,   g = cd4_potential - cd4_preart
```

At `dur_art = 0` this returns `cd4_preart`; as `dur_art → +∞` it returns
`cd4_potential`. But as **`dur_art → -∞` it tends to `2*cd4_preart -
cd4_potential`, which is negative whenever `cd4_potential > 2*cd4_preart`.** A
negative `dur_art` means `self.ti < ti_art` — an agent flagged on-ART before its
scheduled start. Note `initial_cases_diagnosed` are given `ti_art = 0`.

**This is algebra, not evidence.** It has not been observed firing. 028 is built
to settle it from data rather than argument.

## Observations

**obs 1 — Guessing a third narrow fix is the wrong move.** Two attempts have now
cost ~45 min of raccoon and two experiment folders, each ending at the same
traceback. The response in 028 is to enforce the invariant at the lookup itself,
where *every* path must pass, and to log the offending agents' state flags so
the source is identified from data.

**obs 2 — The two mechanical corrections 027 carried both worked.** The log
header printed `emulated every wave (simultaneously): ...` instead of 025's stale
per-wave labels, and `point_key` hashing the model tag was verified to produce
different keys for identical parameters under v1.3 and v1.4. Neither is a
scientific result, but both are now known-good and carry into 028.

**obs 3 — Cost of the crash is bounded and the cache is not poisoned.** The run
died in wave 2's simulation phase with wave 1 committed. Because `point_key` now
includes the model version, nothing from 027 can be silently reused by a run on a
different model.

**obs 4 — Compute was lost to spot reclamation before results were retrieved.**
Raccoon was deallocated within six days, taking the run directory with it from
reach (disk persists, but the machine must be restarted to read it). 025's
partial ensemble had been pulled; 027's was not. No loss of value, since 027
produced no analysable output — but the lesson is to pull on crash, not only on
success.

## What survives

Nothing analysable. 027 produced no NROY, no recovery table, no figures. The run
log was read for the traceback and the diagnosis above is its only output.

## Acceptance

Accepted as a failed run whose failure narrowed the problem: the acute-phase
defect is fixed and confirmed insufficient, and the remaining search space is
small and named.

## Next

**028 runs the same experiment against model-v1.5**, which adds a non-negative
CD4 invariant enforced at both positional rate lookups (`make_p_hiv_death` and
`get_art_mortality_hazard`), clipping to `CD4_FLOOR = 1.0` — upstream's own
`cd4_end = 1` convention, not an invented number — and **records every clip**
with the agents' state flags (`acute` / `latent` / `falling` / `on_art` /
`art_discontinued`), plus `min_dur_art` and the fraction with negative
`dur_art`, which is the hypothesis above stated as a measurement.

So 028 both completes and diagnoses. If `guard_events` comes back empty, the
invariant was never needed and something else changed; if it comes back
dominated by on-ART agents with negative `dur_art`, the hypothesis is confirmed
and the upstream fix is a second one-liner.

**Unchanged from 025's Next:** the questions, the features, σ, the design, and
the gate reading order — including `age_gap_sd_mult`, per
[026's SUMMARY](../026_scenarios/SUMMARY.md).

**Still outstanding, unrelated to the crash:** upstream the acute-floor fix
(patch drafted, PR not opened), and pin `==` in `pyproject.toml` so a fresh
`uv sync` cannot resolve a different model stack.
