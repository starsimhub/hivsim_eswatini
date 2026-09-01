# Exp 021 — Viral suppression as an input: removing the `vls_coverage = 1.0` default

## Question

**What does the model do when ART patients are allowed to be unsuppressed?**

`sti.ART` gained a `vls_coverage` argument in stisim 1.5.11 (PR #577,
`art_vls_ingest`) — the fraction of ART initiators who achieve viral
suppression rather than non-suppressive ART.
[017](../017_version_bump/SUMMARY.md) flagged it in its porting table and
explicitly scoped it out ("any use of `vls_coverage` … beyond making the model
run"). It has never been set. The default is **1.0: every ART initiator is
virally suppressed**, and the consequence is not cosmetic:

| state | transmission efficacy | source |
|---|---|---|
| `on_effective_art` | 0.99 | PARTNER/PARTNER2 (U=U) |
| `on_nonsuppressive_art` | 0.35 | Quinn, Rakai 2000 |

So an unsuppressed patient retains **65× the residual transmission** of a
suppressed one. Running at 1.0 means the model has been treating every treated
agent as effectively non-infectious for the whole calibration window.

This is structurally the same defect [018](../018_adopt_and_size/SUMMARY.md)
removed for PrEP: an inherited default that is not neutral and that nobody
chose. Unlike PrEP, this one is *observed three times over* by PHIA, so it
belongs as an input rather than a calibration parameter.

## What the data says, and why the expected effect is modest

Verified against primary sources and transcribed in
[`vls_construction.py`](../../vls_construction.py) →
[`data/eswatini_vls.csv`](../../data/eswatini_vls.csv):

| suppression **among those on ART**, 15+ | men | women | source |
|---|---|---|---|
| 2016–17 | 0.913 | 0.922 | SHIMS2 Table 9.3.A |
| 2021 | 0.967 | 0.959 | SHIMS3 Table 8.1 |

The rise in *population* suppression over the same period is far larger —
62–75% to 82–89% among all PLHIV — but that is almost entirely the **ART
coverage ramp, which the model already takes as an input** from
`data/art_coverage.csv`. The error introduced by the 1.0 default is therefore
the missing **4–9% of treated patients**, not the 60-point population gap.

Expected effect: mean residual transmission from ART patients rises from 0.01
to roughly 0.92 × 0.01 + 0.08 × 0.65 ≈ **0.06**, a sixfold increase on a small
base. Whether that is visible in prevalence depends on how large the treated
population is, which by 2021 is most of the PLHIV — so the effect should grow
across the window and be largest in the years that matter most for the
decision question.

**No age stratification.** The derivation in `vls_construction.py`
(`derive_given_art_by_age`) divides age-specific VLS-among-PLHIV by
age-specific ART coverage and finds suppression-given-ART **flat in age** to
within 0.3–0.8 percentage points across 15–44 in 2016, with the derived mean
matching PHIA's directly measured sex figure to 0.15 pp. The steep age gradient
in population suppression is a *coverage* gradient, already represented. So
`vls_coverage` is stratified by sex and year only, and that is an evidence-based
choice rather than a data limitation.

## Plan

**Two arms, fixed parameters and seeds** — the lineage of 015–019.

| arm | `vls_coverage` | isolates |
|---|---|---|
| **A** `vls_1.0` | `None` (upstream default) | current model behaviour, the control |
| **B** `vls_phia` | PHIA sex × year series | the effect of the correction |

10 seeds (018 obs 5: CV 4.4% at this parameter point), N = 10 000, 1985–2026,
high-transmission parameters — continuity with 016–019 and the low-variance
point, so effect sizes are readable.

**Pre-2016 assumption, stated rather than buried.** The surveys begin in
2016; ART begins in the model around 2004. `vls_coverage` is held flat at the
2016 value back to the start. There is no measurement before SHIMS2, and
early-ART-era suppression was plausibly *worse* than 0.92, so holding flat is
conservative with respect to the change being tested — it understates the
effect rather than overstating it. If the result turns on this, it becomes an
arm in a follow-up rather than a silent choice here.

**`vls_coverage` defaults any stratum it is not given to 100%**, which would
silently reintroduce the very default being removed. The input table must cover
both sexes across the full adult age range, and arm B asserts this before
running.

## Metrics

1. **Population viral suppression vs PHIA** — the validation that matters, and
   it is not circular in the fitting sense: coverage and suppression-given-ART
   are both inputs, so VLS *among all PLHIV* is an emergent output. If arm B
   reproduces SHIMS2's 70.8% and SHIMS3's 86.6% (ages 15–49), the cascade is
   internally consistent. Needs the new `Cascade` analyzer, since stisim
   carries no result for the effective/non-suppressive split.
2. **HIV prevalence 15–49, and by age and sex vs PHIA.** 019 left the model
   biased low by 4 points; more transmission from unsuppressed patients should
   push against that. Direction and size both matter.
3. **Incidence trajectory.** The mechanism — this is where extra transmission
   shows up first. **The 2021 incidence hold-out is not used**; it is the
   validation set and this experiment does not touch it.
4. **AIDS deaths vs UNAIDS.** 019 closed 39% of the deficit at a prevalence
   cost. Non-suppressive ART plausibly worsens survival as well as raising
   transmission, so this could move without the trade-off — or not.
5. **Cascade composition** — counts on effective vs non-suppressive ART over
   time, as a check that the input is doing what it claims.

## Success criteria

- **Clean:** arm B reproduces PHIA population VLS within its confidence
  intervals, prevalence moves toward the PHIA points, and the correction is
  adopted into model-v1.2.
- **Null and useful:** the effect is below seed noise. Then the 1.0 default was
  harmless, we know the size of what we are ignoring, and it is adopted anyway
  on provenance grounds — the same argument 018 made for PrEP.
- **Awkward:** the model's population VLS misses PHIA badly even with both
  inputs set, which would mean the ART coverage input, the diagnosis rates, or
  the ART initiation logic is off — a cascade problem the decision question
  cannot afford to leave unresolved.

## Not in scope

- Non-suppressive ART's effect on *mortality* as a separate mechanism. Metric 4
  observes whatever the existing `get_art_mortality_hazard` does with it; this
  experiment does not modify it.
- Re-opening the age gradient. Settled above on evidence.
- ART coverage itself, which remains an input. But note the provenance problem
  recorded in `vls_construction.py`: the 2021 rows of the source ART file are
  *back-derived* from SHIMS3 VLS divided by the national suppression rate
  (`Count` is NaN for 2021, values carry full float precision, and
  0.793535 × 0.959 reproduces SHIMS3's women 15–24 VLS exactly). That deserves
  its own experiment.
- Any calibration run, prior bounds, or parameter list — that is 022.
