# Exp 027 — History matching, wave 2 (+ re-identification pre-flight), on model-v1.4

> A re-run of [025](../025_hm_wave2/SUMMARY.md) against a fixed model. Same
> questions, same features, same σ, same design, same pre-registered success
> criteria. The only scientific change is the model version; everything else is
> deliberately held constant so this is a clean second attempt rather than a new
> experiment wearing 025's clothes.

## Question

Unchanged from 025, which never answered any of them:

1. **Can the pipeline separate `beta_m2f` from `s_f_young` at all?** 024's PC1
   showed they enter the wave-1 feature as a *product*. If a synthetic test with
   wave-2's features can't recover them separately, a real wave-2 NROY that looks
   informative in those directions would be misleading.
2. **Does incidence — a rate — constrain what prevalence, a stock, cannot?**
3. **Does the NROY tighten when σ_disc is halved, or collapse?**

## Why now, and what changed

025 crashed 22 minutes in, at the simulation phase of wave 2 of the pre-flight,
on an `IndexError` from a positional death-rate lookup against a **negative CD4**.
Upstream `sti.HIV.acute_decline` was the one of three CD4-decline functions with
no floor, so an agent lingering in the acute state past the timestep rounding
declined through its latent set-point and across zero.

Fixed as `hiv_cd4_floor.HIVCD4Floor`, on by default in `run_sims.make_sim`, and
tagged **model-v1.4**. The fix is a no-op where the floor does not bind — an
identical seeded run gives `cum_infections` 82,943 with or without it.

**This was not a parameter-space problem.** None of the seven calibration
parameters touch CD4 dynamics; `s_f_young` is `rel_sus_age`, susceptibility, not
survival. 025's crash was exposure — wave 2 draws inside the NROY, which
concentrates on higher transmission, at N = 20,000 against 024's 10,000. So the
prior, the features and the σ budget all carry forward untouched. See 025's
SUMMARY for the full diagnosis.

## Plan

Identical to 025 except where noted. Both parts in one job: the re-identification
pre-flight (which **gates** interpretation of wave 2), then the waves themselves.

| | |
|---|---|
| points | 1000 per wave, drawn from the previous NROY by the engine |
| N agents | 20,000 |
| replicates | 1 |
| emulator | `bayes_linear`, drop any feature whose R² < 0.8 |
| threshold | 4.0 |
| waves | up to 3, all three features emulated **every** wave |
| truth | 024's design row 868 (the 48/48 draw), simulated at a different seed |
| compute | raccoon, 120 workers, ~40 min for both parts |

Features and σ exactly as 025 set them — adult prevalence
(`prev_15_49_all_mean`), the incidence aggregate (`inc_all_mean`), and the female
15–24 : 35–44 prevalence ratio (`prev_young_old_f_mean`); prevalence σ_disc 0.01,
incidence discrepancy 10% relative. The rationale for each is in
[025's README](../025_hm_wave2/README.md) and is not restated here.

### Three mechanical changes, none scientific

1. **Fresh cache, and the cache key now includes the model version.** 025's 2000
   cached points were produced under model-v1.3 and are invalid here. `OUT_DIR`
   resolves relative to `run.py`, so this folder already gets its own
   `outputs/sims/` — but 025 obs 6 noted the key hashes only parameters, N, stop
   year and seed, so a model change was invisible to it. The key now carries the
   model tag and the stisim version, which is what makes the folder separation
   belt-and-braces rather than the only defence.
2. **The stale per-wave print labels are gone.** 025 printed `wave 1 emulates:` /
   `wave 2 emulates:` while passing all three features as one list, so the log
   contradicted the code and the README. It now prints the actual feature list.
3. **`age_gap_sd_mult` is read explicitly in the gate output**, per exp 026's
   SUMMARY — see below.

## What the gate must be read for, in order

1. **Truth inside the NROY on every parameter**, and — the pre-registered
   question — is the `beta_m2f`/`s_f_young` **ratio** recovered as well as the
   **product**? Product-only means wave 2's marginals in those two directions are
   reportable as jointly constrained only, and wave 3's job becomes finding a
   feature that separates them. That is a legitimate outcome; misreading wave 2
   is not.
2. **`age_gap_sd_mult` recovery.** Added on the strength of
   [026's SUMMARY](../026_scenarios/SUMMARY.md): the decision quantity is what
   PrEP buys AGYW 15–24 after the cascade gap closes, while the gap being closed
   is in men 25–34, so what links them is who AGYW acquire from across the age
   gap. No feature in the set is an age-mixing statistic, which makes
   `age_gap_shift` and `age_gap_sd_mult` the parameters these features inform
   least while the decision depends on them most.
3. **Per-feature emulator R².** Anything under 0.8 is dropped rather than
   carried; three simultaneous features is already against the 1–2 guidance.
4. **NROY fraction.** Under 5%, or the engine's own `NROY space collapsed`
   warning firing at n = 1000: suspect σ_disc = 0.01 before suspecting the model.
   024's scan predicts only 8.4% of draws survive on prevalence alone at that σ.

## Success criteria

Carried over from 025 verbatim, so this cannot be graded on a curve after the
fact:

- **Clean:** pre-flight recovers truth; wave-2 emulator R² > 0.8 on both tier-B
  features; NROY between roughly 5% and 50%; `s_f_young` marginal narrower than
  024's. Wave 3 opens on the F:M ratios.
- **Informative failure:** pre-flight recovers the product but not the ratio.
  Wave 2 still runs, but its `beta_m2f` and `s_f_young` marginals are reported as
  jointly constrained only.
- **Collapse:** NROY < 5%. σ_disc = 0.01 is the first suspect, not the model.
- **Uninformative:** NROY > 70%. The features are not sensitive enough, or the
  incidence discrepancy is too generous.

## Not in scope

- The Bayesian step. NROY is a region; trajectory selection comes after the waves
  converge.
- The 2021 incidence hold-out.
- The age-banded incidence rows (excluded by decision, 2026-09-04).
- Any further model change. 027 runs model-v1.4 unmodified.
- Re-analysing 025's surviving wave-1 ensemble. It was produced under v1.3 and a
  one-wave NROY does not answer the gate's question anyway.
