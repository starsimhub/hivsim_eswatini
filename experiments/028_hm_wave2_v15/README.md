# Exp 028 — History matching, wave 2 (+ re-identification pre-flight), on model-v1.5

> Third attempt at the same experiment. [025](../025_hm_wave2/SUMMARY.md) and
> [027](../027_hm_wave2_v14/SUMMARY.md) both died at the same traceback. Same
> questions, features, σ, design and success criteria as both; the model version
> is again the only scientific change.

## Question

Unchanged, and still unanswered after two attempts:

1. **Can the pipeline separate `beta_m2f` from `s_f_young` at all?**
2. **Does incidence — a rate — constrain what prevalence, a stock, cannot?**
3. **Does the NROY tighten when σ_disc is halved, or collapse?**

## Why a third attempt, and what is different

Both previous runs crashed in wave 2 of the pre-flight on

```
IndexError: index 6 is out of bounds for axis 0 with size 6
```

from `sti.HIV.make_p_hiv_death`, which indexes a 6-element `cd4_death_rates`
positionally via `np.digitize` against **descending** bins. Index 6 is reachable
only when `cd4 < 0`.

**025** identified and fixed a real defect: `acute_decline` was the only one of
three CD4-decline functions with no floor, so an agent staying acute past the
timestep rounding overshot its set-point by ~350 CD4 on the first extra step and
crossed zero on the second. **027 proved that fix was active and insufficient** —
the fix commits are ancestors of the launched tree, `make_sim` defaulted to the
patched class, and the crash was identical. So a second path to negative CD4
exists and has not been identified.

**028 stops guessing.** model-v1.5 enforces a non-negative CD4 invariant at both
positional lookups — `make_p_hiv_death` and `get_art_mortality_hazard` — which
every path must pass through, whatever produced the value. CD4 below the floor is
clipped to `CD4_FLOOR = 1.0`, which is upstream's own
`cd4_end = 1  # To avoid divide by zero problems`, not an invented constant.

## The diagnostic, which is the point of this experiment beyond the re-run

Clipping alone would hide the defect. So every clip is recorded with the
offending agents' state flags, and a per-point summary is written into each
point's parquet:

| column | meaning |
|---|---|
| `n_cd4_guard` | number of clip events in the run |
| `n_cd4_clipped` | total agent-clips |
| `cd4_guard_min` | most negative CD4 seen |
| `cd4_guard_f_on_art`, `_f_acute`, `_f_latent`, `_f_falling`, `_f_art_discontinued` | agent-weighted state mix of the clipped |
| `cd4_guard_f_dur_art_negative`, `cd4_guard_min_dur_art` | the hypothesis below, as a measurement |

**Hypothesis under test — `cd4_increase`, the on-ART recovery logistic:**

```
cd4_now = 2*g/(1 + exp(-dur_art*growth)) - g + cd4_preart,   g = cd4_potential - cd4_preart
```

At `dur_art = 0` this returns `cd4_preart` and as `dur_art → +∞` it returns
`cd4_potential` — both fine. But as **`dur_art → -∞` it tends to
`2*cd4_preart - cd4_potential`, negative whenever `cd4_potential > 2*cd4_preart`.**
A negative `dur_art` means `ti < ti_art`, an agent flagged on-ART before its
scheduled start. This is algebra, not evidence; `cd4_guard_f_dur_art_negative` is
how it gets confirmed or killed.

**Expected value of every guard column is zero.** A non-zero count is a model
defect to push upstream, not a result to accept. If the columns come back empty
and the run completes, that itself needs explaining — it would mean something
other than the invariant changed the outcome.

## Plan

Identical to 027 except the model version. Both parts in one job: the
re-identification pre-flight, which **gates** interpretation of wave 2, then the
waves.

| | |
|---|---|
| points | 1000 per wave, drawn from the previous NROY by the engine |
| N agents | 20,000 |
| replicates | 1 |
| emulator | `bayes_linear`, drop any feature with R² < 0.8 |
| threshold | 4.0 |
| waves | up to 3, all three features emulated **every** wave |
| truth | 024's design row 868, simulated at a different seed |
| compute | raccoon, 120 workers, ~40 min |

Features and σ exactly as 025 set them; the rationale is in
[025's README](../025_hm_wave2/README.md) and is not restated. `point_key` hashes
the model tag, so 027's cache cannot be reused here — 028 gets a cold cache.

## What the gate must be read for, in order

0. **The guard columns first.** If `n_cd4_guard` is non-zero, read the state mix
   before anything else — it names the second defect, which is the blocker that
   has cost two experiments.
1. **Truth inside the NROY on every parameter**, and the pre-registered question:
   is the `beta_m2f`/`s_f_young` **ratio** recovered as well as the **product**?
   Product-only means wave 2's marginals in those directions are reportable as
   jointly constrained only.
2. **`age_gap_sd_mult` recovery**, per [026's SUMMARY](../026_scenarios/SUMMARY.md)
   — the decision quantity is what PrEP buys AGYW 15–24 after the cascade closes
   while the coverage gap is in men 25–34, so the link is who AGYW acquire from
   across the age gap, and no feature in the set is an age-mixing statistic.
3. **Per-feature emulator R²** — drop anything under 0.8 rather than carrying it.
4. **NROY fraction** — under 5%, or `NROY space collapsed` at n = 1000, suspect
   σ_disc = 0.01 before the model.

## Success criteria

Carried over verbatim from 025 so this cannot be graded on a curve:

- **Clean:** pre-flight recovers truth; wave-2 emulator R² > 0.8 on both tier-B
  features; NROY roughly 5–50%; `s_f_young` marginal narrower than 024's.
- **Informative failure:** pre-flight recovers the product but not the ratio.
- **Collapse:** NROY < 5% — σ_disc is the first suspect, not the model.
- **Uninformative:** NROY > 70%.

## Not in scope

- The Bayesian step; NROY is a region, trajectory selection comes after.
- The 2021 incidence hold-out.
- The age-banded incidence rows (excluded by decision, 2026-09-04).
- Any model change beyond the invariant. 028 runs model-v1.5 unmodified.
- Re-analysing 025's surviving wave-1 ensemble (produced under v1.3).
