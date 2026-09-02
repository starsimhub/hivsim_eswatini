> **Intended as the last experiment before calibration.** After this, 023 is
> parameter engineering and 024 is coverage check v3.

# Exp 022 — A pivoting survival gradient, and suppression that can improve for people already on ART

## Question

Two model changes, both prompted by checking this project's inputs against the
**Akullian et al. 2020** EMOD model of eSwatini (Lancet HIV; appendix at
`external_data/mmc1.pdf`) — the researcher's own prior published work on the
same epidemic.

**1. Does a *pivoting* age gradient on untreated survival fix the age-shape
defect that six previous experiments could not?**

[019](../019_age_dependent_survival/SUMMARY.md) concluded an age gradient was
not worth a calibration parameter, because its half-strength and full-strength
arms were indistinguishable (78.2% vs 78.1% of the UNAIDS death peak). **That
conclusion was drawn from too narrow a family.** Every arm in 019 had
multipliers ≤ 1.0 — all four *shortened* survival at every age. The EMOD model
does something different:

```
λ = 21.182 − 0.2717 × age_at_infection,   Weibull shape 2
mean survival = λ·Γ(1.5)  →  13.96 y at 20, 11.55 at 30, 9.14 at 40, 6.73 at 50
```

Against stisim's flat 13.1 y nominal, that **pivots**: young survival is
*longer*, old survival much shorter. Converted to a latency multiplier at the
019 band midpoints, holding the non-latent portion (`dur_acute` 1.7 months +
`dur_falling` 3 y = 3.14 y) fixed:

| age at infection | 019 arm D | **EMOD-derived** |
|---|---|---|
| 15–24 | 0.94 | **1.086** |
| 25–34 | 0.84 | 0.845 |
| 35–44 | 0.64 | 0.604 |
| 45+ | 0.44 | 0.363 |

This matters because it pushes **both** of the model's prevalence defects the
right way at once. The retrofit decomposition
([`plot_fit_progression.py`](../../plot_fit_progression.py)) shows young men at
−65% and young women at −49% relative to PHIA, while women 35–44 sit **+11%
above** it. Longer young survival adds prevalent infections where the model is
low; shorter old survival removes the excess where it is high. Every previous
attempt moved the whole curve in one direction and so could only trade one
defect against the other — which is exactly what 019 measured as "the trade-off
bites".

**2. Can viral suppression improve for people already on ART?**

`sti.ART(vls_coverage=...)` accepts a time-varying series, but applies it **only
at initiation** — `on_effective_art` is written in `start_art()` and cleared
only on death or ART discontinuation. So an agent who started in 2010 keeps
their 2010 suppression status for life.

That makes the actual mechanism unrepresentable. Eswatini's third-90 gain came
substantially from **better regimens and adherence support** — the dolutegravir
(TLD) transition from ~2019 — which improve suppression among *existing*
patients, and by 2021 existing patients are most of the treated population.
021 measured the symptom: realized suppression lagged its own input, 0.956
against 0.967 in 2021.

Structurally this is the VMMC defect from [015](../015_vmmc_prevalence_target/README.md)
again — coverage as a per-step hazard on new entrants rather than a stock
target — and [`vls_stock_target.py`](../../vls_stock_target.py) takes the same
shape of fix, ranking the treated stock by a persistent per-agent adherence
propensity exactly as `sti.VMMC` ranks by willingness.

Measured in a 3 000-agent smoke test, the difference is real:

| | 2016 | 2021 | 2025 |
|---|---|---|---|
| input | 0.905 m / 0.918 f | 0.967 / 0.959 | held flat |
| flow only | 0.935 / 0.935 | 0.956 / 0.955 | 0.981 / 0.963 |
| stock target | 0.912 / 0.921 | 0.966 / 0.959 | 0.966 / 0.959 |

Flow-only drifts *above* the input early (survivorship: unsuppressed patients
have higher on-ART mortality, so the surviving stock is enriched for suppressed)
and *below* it later (the target rises faster than the stock turns over).

## Plan

**2 × 2 factorial**, fixed parameters and seeds — the lineage of 015–021. Ten
seeds (018 obs 5: CV 4.4% at this point), N = 10 000, 1985–2026,
high-transmission parameters for continuity and low variance.

| arm | survival | VLS | isolates |
|---|---|---|---|
| **A** `base` | flat (upstream) | flow only | model-v1.2, the control |
| **B** `pivot` | EMOD gradient | flow only | **the survival pivot** |
| **C** `stock` | flat | stock target | **the suppression capability** |
| **D** `both` | EMOD gradient | stock target | interaction, and the candidate v1.3 |

A factorial rather than a ladder, because the two changes could interact:
shortening old-age survival removes treated patients, which changes the treated
stock the suppression target acts on.

## Metrics

The standard fit figure via
[`standard_figures.plot_prevalence_fit`](../../standard_figures.py) for every
arm, plus:

1. **Prevalence by age and sex vs PHIA** — the target. Specifically whether the
   young-age deficit narrows *without* inflating women 35–44. Reported as the
   bias decomposition by band, not just MAE, since MAE hides the shape.
2. **AIDS deaths vs UNAIDS.** 019's arms bought deaths at the cost of
   prevalence; the pivot should buy deaths at the *older* ages while giving
   prevalence back at the younger ones. If that holds, the trade-off 019
   identified is an artefact of only testing compressions.
3. **Realized survival by age at infection**, to confirm the pivot does what the
   EMOD parameterisation says — ~14 y at 20 falling to ~7 y at 50.
4. **Cascade vs PHIA** (`analyzers.Cascade`): population VLS among PLHIV, and
   suppression among the treated against the input series. Arm C and D should
   track the input; A and B should drift.
5. **Rare-event floor** — 020 found the thin strata are thin partly because the
   model under-infects young men. If the pivot raises young prevalence, the floor
   improves too, which would be a second reason to prefer it.

## Success criteria

- **Clean:** the pivot narrows the young-age deficit and reduces or holds the
  35–44 excess, moving age-stratified MAE below the 0.0586 that has been flat
  since 018. Adopted as model-v1.3, and calibration opens against a model whose
  age structure is right in shape.
- **Partial:** the pivot helps deaths and young prevalence but overshoots
  somewhere else. Still valuable — it converts a structural mismatch into a
  parameter with known direction, which is what `dur_latent_mult` becomes in
  023.
- **Null:** the pivot behaves like 019's compressions and moves nothing beyond
  noise. Then untreated survival is exhausted as an explanation, the age-shape
  defect belongs to the network (age mixing, `rel_sus_age`), and 023 opens those
  instead.
- **Either way the suppression capability is adopted**, on the same argument
  018 used for PrEP and 021 for the VLS default: it is a mechanism the model
  should have regardless of whether it moves the fit, because the decision
  question is about raising suppression.

## Not in scope

- The VLS definitional basis. Checked and dropped: SHIMS2's self-report and
  biomarker figures differ by under a point at matched 15+ (91.3/92.2 vs
  90.5/91.8). `vls_construction.py` now uses the biomarker 15+ values for
  consistency with SHIMS3, but it is a tidy, not a correction, and not an arm.
- `rel_sus_age` and `rel_beta_f2m`. The EMOD appendix gives female
  susceptibility 4.894 (<25) and 2.844 (≥25) — a ratio of 1.72, which is where
  this model's 1.7 came from, but at levels ~1.4× lower than this model's
  implied 6.8 and 4.0. Reconciling that is **023**, parameter engineering, not
  here: opening it in the same experiment as the survival pivot would confound
  two age-shape mechanisms.
- The 2016 young-male ART coverage discrepancy (0.360 in the input vs 0.556
  derived from SHIMS2 Table 10.3.B). Both SHIMS2 inputs are parenthesised
  small-denominator estimates and the Akullian appendix supplies no age/sex ART
  coverage data — that model *fitted* coverage via `ART_Link_Max/Mid/Rate`. Left
  open with the provenance recorded in
  [`art_coverage_construction.py`](../../art_coverage_construction.py).
- Any calibration run.
