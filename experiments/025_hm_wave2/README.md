# Exp 025 — History matching, wave 2 (+ re-identification pre-flight)

> Two parts. A synthetic-data recovery test that **gates** the wave, then wave 2
> itself: sex-specific incidence and the female young:old prevalence ratio, at
> half the model-discrepancy allowance [024](../024_hm_wave1/SUMMARY.md) used.
> Sex-specific incidence enters the target set here for the first time.

## Questions

1. **Can the pipeline separate `beta_m2f` from `s_f_young` at all?** 024's PC1
   showed they enter the wave-1 feature as a *product*. If a synthetic test with
   wave-2's features can't recover them separately, a real wave-2 NROY that looks
   informative in those directions would be misleading.
2. **Does incidence — a rate — constrain what prevalence, a stock, cannot?**
3. **Does the NROY tighten when σ_disc is halved, or collapse?**

## Part 1 — re-identification pre-flight (workflow step 7)

**Gating.** Wave 2 does not run until this reports. Recorded as a gate rather
than a companion because the failure mode is specifically a *misreading* of
wave 2, not a wave-2 failure.

**Why this is finally affordable.** 023 flagged it as the largest gap and
deferred it. 024 deferred it again by researcher decision. But 024 also made it
cheap: 1000 points in **456 s** on raccoon's 120 cores, against the 4.5 h the
024 README budgeted on 8 workers. There is no longer a cost argument for
skipping it.

**Design.** Take 024's design row 868 as known truth — the point that fit 48/48
targets, so it is both plausible and already characterised:

| parameter | truth |
|---|---|
| `beta_m2f` | 0.0142 |
| `rel_beta_f2m` | 0.3097 |
| `s_f_young` | 2.2240 |
| `age_gap_shift` | +0.07 |
| `age_gap_sd_mult` | 1.5001 |
| `prop_f0` | 0.6088 |
| `prop_m0` | 0.4742 |

Generate synthetic observations by running that point at a **different seed**
from the one that produced it, so the target carries one replicate's worth of
noise rather than being the model's own exact output. Then run the same
waves, the same features, the same σ, and ask what comes back.

**Success:** truth inside the NROY on every parameter, and the
`beta_m2f`–`s_f_young` ridge visibly narrower than 024's.

**The failure that matters:** truth recovered on the *product* `beta_m2f × s_f_young`
but not on the *ratio*. That would say tier B is insufficient to separate them,
and the response is to find a feature that can — or to accept the degeneracy,
report the correlation, and stop pretending the marginals are interpretable.
Either is a legitimate outcome; misreading wave 2 is not.

**Non-goal.** This is not a test that the model is right. Synthetic data comes
from the model, so structural error cancels by construction. It tests the
*pipeline* — prior, features, emulator, implausibility — and nothing else.

## Part 2 — wave 2

### Features — and a correction to how "wave" is used here

**Corrected 2026-09-07, after the smoke test.** This README originally described
wave 1 emulating the prevalence aggregate and wave 2 adding the two tier-B
features. **The package does not work that way.**
`ManualFeatureSelection(list)` emulates *every* listed feature in *every* wave —
it is not a per-wave schedule. The smoke-test log is unambiguous:

```
[Wave 1] Phase 3/5 FEATURE SELECTION: ['prev_15_49_all_mean', 'inc_all_mean', 'prev_young_old_f_mean']
[Wave 1] Phase 4/5 EMULATOR TRAINING: training 3 emulators...
```

So all three features are applied **simultaneously**, and successive waves
re-draw the design inside the surviving NROY rather than layering a new check.
A "wave" here means an NROY-refinement iteration.

**Kept rather than worked around, for a reason.** The question 024 left open is
whether the model can satisfy prevalence *and* incidence *and* the female
young:old ratio **jointly** — 024's headline was a joint 48/48 fit, and
simultaneous emulation tests exactly that. Sequential layering would answer
"does each successive check cut further", which is the weaker question here.

**The cost, stated.** This deviates from the history-matching guidance of 1–2
features per wave, on the grounds that more features dilute the signal. Two
mitigations: the package writes per-feature emulator metrics and z-scores, so
attribution of which feature did the cutting survives; and any emulator
returning R² < 0.8 gets dropped rather than carried.

The three features and their jobs:

| feature | job |
|---|---|
| **`prev_15_49_all_mean`** — adult prevalence, 024's feature | Carried so the cut remains comparable to 024's, which used it alone. |
| **`inc_all_mean`** — mean of the four fitted sex-specific incidence aggregates | Everything the transmission parameters govern is a **rate**. Prevalence is a **stock** that confounds rate with mortality, ART and cohort history — which is why 016–022 could not diagnose the age-shape defect from prevalence alone. This is the reason incidence was added to the target set. |
| **`prev_young_old_f_mean`** — female 15–24 : 35–44 prevalence ratio, over the three PHIA years | The one statistic that separates `s_f_young` from `beta_m2f`. 023 obs 4 measured ρ = +0.36 against women 15–24 and +0.04 against women 35–44 — specific, not a level parameter in disguise. Directly targets the degeneracy Part 1 tests. |

**Deferred to wave 3:** the F:M prevalence and incidence ratios, which pin
`rel_beta_f2m`. That is 024 obs 6's open question — 023 read 0.20–0.30 off the
F:M incidence ratio alone, while 024's joint fit across 48 targets preferred
0.30–0.49 (median 0.40). Worth its own wave rather than being crowded in here.

### Targets

**Incidence is four rows, sex-specific aggregates only** — 2011 (SHIMS1 cohort,
18–49) and 2016 (SHIMS2 Table 5.3.B, 15–49), by sex. The six age-banded rows
were removed from the target sheet on 2026-09-04 and live in
`data/shims2_incidence_by_age_REFERENCE_NOT_A_TARGET.csv`.

Reason, recorded because it is a judgement about the data and not about the
model: the published female age profile is essentially flat (1.67, 1.54, 2.09
across 15–24, 25–34, 35–49), which is the expected signature of **false-recency
bias**. A recency assay credits a fraction FRR of long-standing infections as
recent, and against a susceptible denominator the spurious incidence scales as
FRR × P/(1−P)/MDRI — 5× larger in women 35–49 than 15–24, and 18× larger in men.
At FRR = 0.3% that accounts for ~40% of the reported value in both older female
bands. See `incidence_construction.py`'s `FIT_POLICY`.

**Consequence to state plainly: the model's incidence age profile is now
unconstrained by data.** 024's `figures/incidence_age_profile.png` shows what it
does anyway (female peak at 22.5, male at 32.5).

### σ — the decision that will make or break this wave

| family | σ | change from 024 |
|---|---|---|
| prevalence | sqrt(σ_CI² + **0.01**²) | **halved** from 0.02 |
| incidence | sqrt(σ_CI² + (0.10 × target)²) | discrepancy term **added** |
| prevalence ratios | 0.15 relative | unchanged |
| F:M incidence ratio | 0.35 | unchanged |
| peak AIDS deaths | 2000 | unchanged (the deliberate down-weighting) |

**Prevalence σ_disc 0.02 → 0.01.** 024 set 0.02 to absorb a −4.3 pp prevalence
deficit that seven experiments had diagnosed as structural — and then showed it
was not structural at all. At σ_disc = 0 the best draw sat at 3.98σ on the six
PHIA 15–49 targets, so the deficit was a parameter-value deficit measured at
default transmission parameters. 0.02 was calibrated against a bias that does not
exist, and it dominated the observation variance (0.02² against PHIA's
0.0043–0.0066²). Halving it is worth roughly 5× the cutting power.

Not to zero, though: 024 obs 8 leaves male prevalence at 15–24 genuinely
unreachable, so some allowance is real.

**Incidence gets a 10% relative discrepancy** (researcher decision 2026-09-07).
A modelled incidence and a recency-assay or cohort estimate are not the same
quantity, and prevalence already carries an explicit allowance — without one,
incidence would bite harder than prevalence and the two families would not be on
comparable footing. On 2016 men (0.85, σ_CI 0.327) it adds 0.085 in quadrature,
so it is a modest widening, not a down-weighting.

### Design

| | |
|---|---|
| points | 1000 per wave, each wave drawn from the previous NROY by the engine |
| N agents | **20,000** (up from 10,000) |
| replicates | 1 per point |
| emulator | `bayes_linear`; fall back to `gpr` if R² < 0.8 |
| threshold | 4.0 |
| waves | up to 3 (`--max_waves`), all three features every wave |

**N = 20,000, and this is why it matters here and did not in 024.** 020
recommended 20,000 because two PHIA strata fall below 5 expected infected agents
at 10,000. Wave 1 emulated an aggregate, where 10,000 was ample. Wave 2's
features are age-stratified (a 15–24 : 35–44 ratio) and incidence numerators are
small, so the noise floor now bites.

## Why this cannot resume from 024, and what that means for the record

**Halving σ_disc means wave 1's cut must be redone.** 024's NROY was computed at
σ_disc = 0.02; layering a wave 2 at 0.01 on top of it would mix two uncertainty
budgets in one NROY. So 025 runs as a fresh two-wave job.

**024's 87.2% NROY is therefore superseded, not built upon.** Stated here so the
025 SUMMARY does not read as contradicting 024. Everything else 024 established
— the emulator quality, the coverage result, the 48/48 joint fit, the
constrained-directions diagnostic — stands, because none of it depended on the
NROY volume.

Cost is small: wave 1's 1000 points re-simulate in ~8 min at N = 20,000.

## Corrections landing with this experiment

Three, all found while planning rather than while debugging:

1. **`run.py` must read `sigma` and `model_age_basis` from the incidence CSV.**
   024's version recomputes σ from `(ub−lb)/3.92` and defaults the model side to
   ages 15–50. Left alone it silently undoes the 2026-09-04 decision and compares
   SHIMS1's 18–49 estimate against a model 15–49 one — which flipped the sign of
   the 2011 male residual, from −0.67 to +0.26.
2. **The point cache key must include `N_AGENTS`.** `point_key()` hashes only the
   seven parameter values. Wave 1 ran at N = 10,000; at 20,000 the cache would
   return **stale 10,000-agent results** for any repeated point. Invisible, and
   it would have corrupted the wave. 025 also uses its own `outputs/sims/`.
3. **Tier-B features added to `summarise_point`**: incidence aggregates computed
   on each target's own age basis (pro-rating the 15–19 band to 2/5 for the 2011
   rows), and the F:M ratios for wave 3.

## The engine's own collapse tripwire

The smoke test surfaced this and it is worth recording, because it changes how a
short run should be read: **the engine stops early if fewer than 20 NROY samples
can be found**, logging

```
NROY space collapsed: only N samples found (need 20+). ... The model may be
over-constrained -- consider relaxing the implausibility threshold or increasing
observation uncertainty (model discrepancy).
```

At n_samples = 12 that fired trivially and is meaningless. At n_samples = 1000
it would be a genuine collapse signal, and per the success criteria below the
first suspect is σ_disc = 0.01, not the model.

## Success criteria

- **Clean:** Part 1 recovers truth; wave-2 emulator R² > 0.8 on both features;
  NROY between roughly 5% and 50%; `s_f_young` marginal narrower than 024's.
  Wave 3 opens on the F:M ratios.
- **Informative failure:** Part 1 recovers the product but not the ratio. Wave 2
  still runs, but its `beta_m2f` and `s_f_young` marginals are reported as
  jointly constrained only, and wave 3's job becomes finding a feature that
  separates them.
- **Collapse:** NROY < 5%. σ_disc = 0.01 is the first suspect, not the model —
  and 024's scan already predicts 8.4% under threshold on prevalence alone, so
  this is a live possibility rather than a remote one.
- **Uninformative:** NROY > 70%. The two features are not sensitive enough, or
  the added incidence discrepancy is too generous.

## Not in scope

- The Bayesian step. NROY is a region; trajectory selection comes after the
  waves converge.
- The 2021 incidence hold-out.
- The age-banded incidence rows (excluded by decision, see above).
- Any change to the model. 025 is a calibration wave against model-v1.3.
