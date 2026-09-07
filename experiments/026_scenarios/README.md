# Exp 026 — Decision scenarios: long-acting PrEP vs the treatment cascade

> For the CROI abstract, due 1 October 2026. Scenario arms layered on the
> calibrated model via `scenarios.py` and
> `make_sim(extra_interventions=...)`.

## Question

From CLAUDE.md: the relative value of scaling up long-acting PrEP versus
improving the treatment cascade to raise population viral suppression — **and
where prevention still matters when suppression is high.**

That last clause is the novel part and it is a subtraction:
(PrEP + high suppression) − (high suppression alone).

## The finding that reframes the whole design

**Eswatini has already achieved the third 95. The cascade headroom is entirely
in coverage.** Measured from the project's own target tables at 2021:

| | value | headroom to 0.95 |
|---|---|---|
| suppression given ART, women | **0.959** | −0.009 |
| suppression given ART, men | **0.967** | −0.017 |
| ART coverage, **men 25–34** | **0.650** | **+0.299** |
| ART coverage, women 15–24 | 0.793 | +0.156 |
| ART coverage, men 15–24 | 0.833 | +0.117 |
| ART coverage, women 25–34 | 0.894 | +0.056 |
| ART coverage, all strata | 0.650–0.970 (mean 0.870) | |

Implied population viral suppression ≈ 0.870 × 0.963 = **0.838** of PLHIV,
against 0.95² = 0.902 on this conditional chain. **The +0.065 gap is essentially
all coverage.**

Three consequences, all of which change what we run:

1. **"Improve the third 95 to 95%" is not a scenario, it is a regression.**
   Suppression among the treated is already 96–97%, so ramping it *to* 0.95
   would lower it. `scenarios._ramp_long` now raises on this rather than running
   it silently — a mistake that would have produced a "cascade improvement" arm
   with worse outcomes than baseline and no error.
2. **The cascade arm has to be ART coverage**, not regimen or adherence. In
   Eswatini, "raise population viral suppression" means *find and treat people*,
   not *suppress the treated better*.
3. **The largest single gap is men 25–34 at 65%** — the same band carrying this
   model's most persistent prevalence deficit through experiments 016–022. That
   is a coherent story rather than a coincidence, and it makes a men-targeted
   cascade arm worth running on its own.

Sex coding verified: `art_coverage.csv` uses stisim's convention, **0 = female,
1 = male**, cross-checked against the 2016 men 15–24 value of 0.556 corrected in
`art_coverage_construction.py`. The opposite convention holds in
`prevalence_by_age_sex.csv`; getting this backwards would have swapped the
headline finding onto women.

## Arms

Scenario window opens in **2026** — after the last calibration target (2016) and
after the 2021 validation hold-out, so no scenario can contaminate either.
All ramp linearly to their final value by 2030 and run to 2040.

| arm | PrEP | cascade | what it answers |
|---|---|---|---|
| `baseline` | none | status quo, held flat | the counterfactual |
| `prep_agyw` | Len, 30% of AGYW 15–24 | status quo | prevention alone |
| `prep_agyw_risk` | Len, 30% of *higher-risk* AGYW | status quo | does targeting within AGYW matter? |
| `prep_fsw` | Len, 60% of FSW | status quo | targeting comparison |
| `art_95` | none | ART coverage → 95%, all strata | treatment alone |
| `art_95_men` | none | ART coverage → 95%, men only | is the men's gap the whole story? |
| `both` | Len 30% AGYW | ART → 95% | additive or redundant? |
| `prep_at_high_art` | Len 30% AGYW | ART **already** at 95% from 2026 | **the actual question** |

`prep_at_high_art` minus `art_95` is the headline quantity: what prevention buys
*after* the cascade gap is closed. The first five arms are table stakes.

## PrEP product

Lenacapavir, 6-monthly, from Akullian et al. (LenOptim, *Lancet HIV* R1),
Methods and Table 1:

- efficacy **0.95** base case, DSA range **0.90–1.00** (refs 29, 30)
- `prep_dur = ss.months(6)`
- `prep_adh = 1.0` — adherence is not the binding constraint for an injectable;
  discontinuation is, and that is expressed by coverage, not adherence

## The uptake-heterogeneity problem, and why two AGYW arms

**Within an eligible group, `sti.Prep` fills its coverage target by ranking on
`willingness` — a *random* per-agent score** (`hiv_interventions.py:830`). So by
default PrEP reaches a random subset of the eligible: relative incidence in
users ≈ 1× the group average.

The LenOptim manuscript's own uptake-heterogeneity analysis found the **>2×
quartile** "most representative of self-selecting PrEP adopters consistent with
the ECHO PrEP implementation study". Random allocation across a broad group
therefore **understates** impact relative to a real programme, and the size of
that understatement is a delivery assumption rather than biology.

With the existing API the way to express self-selection is to narrow
*eligibility* onto higher-risk strata, not to raise coverage — hence
`prep_agyw` (random within AGYW) and `prep_agyw_risk` (AGYW in risk groups ≥ 1).
Reporting both brackets the assumption instead of hiding it inside one number.

## Parameter basis — and the limitation the abstract must state

**Point fit now, ensemble if 025 lands** (researcher decision 2026-09-07).

The calibration will not be finished before 1 October: exp 025 is blocked on
compute, there is no posterior, and wave 3 is unscoped. So arms run at **024's
design row 868** — the draw that fit 48 of 48 targets:

| `beta_m2f` | `rel_beta_f2m` | `s_f_young` | `age_gap_shift` | `age_gap_sd_mult` | `prop_f0` | `prop_m0` |
|---|---|---|---|---|---|---|
| 0.0142 | 0.3097 | 2.2240 | +0.07 | 1.5001 | 0.6088 | 0.4742 |

`run.py` takes a **list** of parameter sets and loops over it, so swapping in the
wave-2 NROY sample is a one-line change if 025 finishes in time. That is the
whole reason for the list-shaped interface.

**Limitation to state in the abstract:** with a single parameter set, differences
between arms carry seed noise only — not parameter uncertainty. The comparison
is internally valid (same parameters, same seeds, one intervention changed) but
has no credible interval. Do not report one.

## Two traps this design avoids, both previously live

1. **Never construct `sti.Prep()` without explicit `coverage`.** With
   `coverage=None` it applies a built-in ramp reaching 80% of FSW **starting in
   2004**. Every experiment from 001 to 017 ran with that phantom programme.
   In a counterfactual comparison it would contaminate the baseline, which is
   worse than contaminating a fit.
2. **The baseline suppression is only correct because of the exp 022 fix.**
   stisim's `vls_coverage=None` means 100% of ART initiators suppressed, which
   exp 021 measured as overstating population suppression by **8.8 pp**.
   Measuring a suppression-raising scenario against a baseline that already
   overstates suppression understates what cascade improvement buys, and so
   **systematically favours PrEP** — the exact bias that would flatter this
   paper's PrEP arm. The PHIA-derived `vls_coverage` plus the stock target is
   what makes the comparison valid, and it belongs in the methods.

## Outputs

| | |
|---|---|
| primary | cumulative HIV infections 2026–2040, per arm; and infections averted vs `baseline` |
| secondary | population viral suppression, incidence 15–49 by sex, AIDS deaths |
| stratified | infections averted by sex and age band, since "where prevention still matters" is an age/sex question |
| figures | infections-averted bar chart with seed spread; incidence trajectories by arm; the `prep_at_high_art` − `art_95` difference |

Analyzers: `PopByAgeSex` (per-band `n_infected`, `n_alive`, `new_infections`)
and `Cascade` (95-95-95 by sex, including the effective/non-suppressive split).

## Not in scope

- Costs, DALYs, cost-effectiveness. The LenOptim manuscript covers that ground
  with a static model; this is the dynamic-transmission complement.
- Parameter uncertainty (see above).
- Any model change. 026 runs model-v1.3 unmodified.
