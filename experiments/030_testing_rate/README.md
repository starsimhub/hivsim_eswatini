# Exp 030 — Slow the model's testing down to the observed awareness, and find out what breaks

**Question.** [029](../029_cascade_audit/SUMMARY.md) found the model diagnoses
people too fast: it over-states awareness by up to 14.6 pp (men 25-34: 0.894
against SHIMS3's 0.748) while reproducing ART coverage among PLHIV to a mean of
1.4 pp, because ART coverage is an input and the linkage step silently absorbs
the error. This experiment asks two things:

1. **Can sex-specific routine testing rates reproduce the observed awareness
   profile?** Two parameters, eight targets — so the fit is over-identified and
   the answer is informative either way.
2. **Does ART coverage survive the fix?** This is the question that makes the
   experiment worth running, and 029 could not ask it.

**Why (2) is the real question.** Awareness is a hard ceiling on ART coverage in
stisim — the stratified target fills only from `diagnosed & ~on_art`
(`hiv_interventions.py:573`). The two SHIMS3 quantities leave almost no slack:

| stratum | aware \| PLHIV | on ART \| PLHIV | headroom |
|---|---|---|---|
| women 50+ | 0.968 | 0.966 | **0.002** |
| women 35-49 | 0.983 | 0.968 | **0.015** |
| men 50+ | 0.970 | 0.953 | 0.017 |
| women 15-24 | 0.836 | 0.807 | 0.029 |
| men 35-49 | 0.934 | 0.903 | 0.031 |
| women 25-34 | 0.931 | 0.900 | 0.031 |
| men 15-24 | 0.911 | 0.877 | 0.034 |
| men 25-34 | 0.748 | 0.650 | 0.098 |

In women 50+ the survey says 96.8% know and 96.6% are treated — **two tenths of
a percentage point of slack.** The model currently hits its ART coverage target
comfortably *because* its awareness is inflated. Bring awareness down to the
observed level and the model must thread a gap it has never had to thread. If it
cannot, that is a finding about the joint feasibility of the two data sources
under this model's structure, not a tuning failure — and it would matter for
every scenario built on the cascade.

**The prior analysis this rests on, and its limits.** A Little's-Law estimate
(mean time unaware = unaware fraction / new-infections-per-PLHIV) implies roughly
**3.8 years** for men 25-34 and **2.8-3.0** for women, against the model's ~1.3.
That calculation is a sanity check, not an input:

- it assumes steady state, which a ten-fold decline in incidence violates, and
  the bias inflates the estimate;
- correcting for infected agents *ageing into* a band moved the women's figure
  from 3.9 to 2.9 and the men's barely at all, so the estimate is sensitive to an
  assumption that is easy to get wrong;
- it ignores the competing risks the model actually has — `low_cd4_testing` at
  0.85-0.95/yr, and mortality.

So the rates are **fitted against the awareness data**, not derived. The
Little's-Law numbers are recorded here to be checked against the fitted answer
afterwards. If the fit lands far from ~2-4 years mean time-to-diagnosis,
something is wrong with one of them and that is worth knowing.

**What 029 ruled out.** An earlier hypothesis — that a never-testing subgroup was
needed — is not supported. The sex difference in awareness (25.2% vs 6.9% unaware
at 25-34) is mostly a stock-age effect: men acquire HIV later, so their PLHIV
stock at 25-34 is epidemiologically young (new infections are 6.45% of the stock
per year against women's 1.76%). Once that is accounted for, the residual
difference in unaware *duration* is about 30%, not 3.6x. A uniform-within-sex
hazard is therefore the right thing to try first.

## Plan

**Model change.** Split `other_testing` in `get_testing_products()`
(`interventions.py`) into male and female interventions, each with its own
eligibility and its own rate. No subclass is needed — `HIVTest` already takes an
`eligibility` callable, so this is a change to the existing factory, not new
machinery. `fsw_testing`, `low_cd4_testing` and `anc_testing` are unchanged;
ANC testing is a genuine sex-specific route and should stay.

**Parameters.** Two multipliers, `k_m` and `k_f`, scaling the existing
`gp_prob` ramp rather than replacing it — so the historical scale-up *shape*
(0 to 0.5 across 1990-2020, then to 0.6) is preserved and only its level moves.
Scan range 0.2-1.4, which brackets 1.0 (no change) and the ~0.5 the
Little's-Law estimate suggests.

Because the ramp is linear from zero, `linspace(0, 0.5) * k` is
`linspace(0, 0.5k)` — **fitting the multiplier is fitting the 2020 plateau**,
with the ramp still anchored at 0 in 1990. Early years scale too but barely in
absolute terms (1995: 0.08 -> 0.04 at k = 0.5), so the change is concentrated in
the recent era.

### Does the testing rate's change over time need fitting too?

Raised by Adam, and it is the weakest assumption here. Three parts to the answer.

**2021 awareness has little leverage on the historical shape.** Awareness
saturates: an agent infected in 2000 has had twenty years of testing hazard and
is ~100% aware by 2021 under almost any plausible ramp. The observed 97% in men
50+ is consistent with a wide range of histories. So the eight strata constrain
mainly the **recent** level, and shape misspecification contaminates that only
mildly. The corollary is the uncomfortable one: **2021 data alone cannot
identify the trajectory.** SHIMS1 (2011) and SHIMS2 (2016) awareness would, and
this is a second independent reason to want those reports.

**Where the shape does bite is historical ART reachability.** A slower ramp
means less awareness in earlier years, and awareness is a hard ceiling on ART
coverage. If 2011 awareness falls below that year's ART target (women 0.18, men
0.08) the target stops filling. Those targets are low so it likely holds — but
029 only checked 2011, 2016 and 2021 because those are the survey years, and
there is no reason not to check every year. **Added to the plan below.**

**Decision: fit the level, diagnose the shape, do not fit the shape.** Adding a
shape parameter on 2021 data alone would be fitting something the data cannot
see. Instead the age-band residual pattern is recorded as a flag (below).

**Design.** Men's and women's rates are nearly separable — each sex's testing
affects the other's awareness only weakly, through transmission — so a 1-D scan
per sex first (~8 values x 3 seeds), then a joint confirmation run at the
selected pair. Stop year **2022**, not 2031: only 2021 is needed for the fit,
which roughly halves the runtime. A final run to 2031 at the chosen values
regenerates 029's ceiling and coverage figures for comparison.

**Targets.** The 8 awareness strata (`aware_of_all_plhiv`, sex x {15-24, 25-34,
35-49, 50+}) from `data/eswatini_cascade_95s.csv`. **Not a calibration target
in the CLAUDE.md sense** — this is a fit of two model-configuration parameters,
not of the epidemiological parameters 023-028 work on, and prevalence/incidence
/deaths are unaffected by construction.

**Instrumentation.** `CascadeByAge` (029) plus `PopByAgeSex`, and the
target-vs-achieved ART coverage check from 029's `analyse.py` — which is now the
primary diagnostic, not a side check. Two extensions to it:

- **Every year, 2004-2030, not just survey years.** A slower testing ramp could
  make an early ART target unreachable, and 029 would not have seen it.
- **Age-band residual pattern**, recorded as the shape diagnostic. If the eight
  residuals are systematically ordered by age band, the assumed ramp shape is
  wrong rather than just its level.

  **Read that flag cautiously.** An age-ordered residual is ambiguous between
  "the testing history is wrong" and "the model's age-at-infection distribution
  is wrong" — and the latter is unconstrained in this project, because the
  age-banded SHIMS incidence rows are `fit=False` for false-recency bias
  (CLAUDE.md). This is the same confounding 028 hit. Report the flag; do not
  diagnose from it.

## Success criteria

**Good.** A `(k_m, k_f)` pair reproduces all eight awareness strata to within a
few points, ART coverage is still met in every stratum, and the implied mean
time-to-diagnosis is in the 2-4 year range the Little's-Law estimate suggests.

**Informative failure — the more likely and more interesting outcome.** Either:

- *ART coverage breaks in the tight strata* (women 50+, women 35-49). That would
  say the model cannot simultaneously reproduce SHIMS3's awareness and SHIMS3's
  ART coverage — most plausibly because the model's PLHIV age/duration structure
  differs from the survey's, so the people it has diagnosed are not the ones the
  coverage target needs. This is the outcome the headroom table predicts.
- *Two rates cannot hit four age bands per sex.* A single hazard per sex implies
  a specific relationship between awareness at 25-34, 35-49 and 50+. If no `k`
  reproduces all three, that is evidence for the heterogeneity 029 ruled
  out on the level data alone — and this is the over-identification paying off.

Either failure is a result. Neither is a reason to re-run with more parameters
inside this folder.

## Notes

- `interventions.py` is on CLAUDE.md's shared-code list (conflict risk with
  Daniel). The change is additive within one function.
- If adopted, this is a model change and needs a new tag (**model-v1.6**) and a
  merge to `main`, per CLAUDE.md's branching note. It invalidates 028's
  calibration — though 028's marginals were the prior, so little is lost.
- Awareness data exists for **2021 only**; SHIMS2's and SHIMS1's reports are not
  in the repo. The fit is therefore anchored at one point in time and the ramp
  *shape* is assumed, not fitted — see the section above for why 2021 data
  cannot test it. Adam is looking for the earlier PDFs;
  `cascade_construction.py` is already set up to extract the same tables if they
  turn up, which would turn this into a trend fit and make the shape
  identifiable. **If they arrive before this experiment closes, stop and
  re-scope** — a two-point or three-point trajectory changes the design, not
  just the precision.
