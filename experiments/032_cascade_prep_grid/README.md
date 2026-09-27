# Exp 032 — A cascade x PrEP surface, with the cascade axis measured rather than labelled

**Question.** [031](../031_scenarios_v16/SUMMARY.md) found that the `art_95` arm
could not reach its target in four of eight strata, and that the model's
Eswatini **already meets 95-95-95 in aggregate** (0.895 of PLHIV suppressed at
2030 against the 0.857 product target) while concealing a 31-point spread
between strata. Discrete arms labelled by what they *asked for* are therefore
the wrong instrument: `art_95` is a scenario whose name is not true of it.

This experiment replaces the discrete arms with a **surface**. Every cell is
plotted at the viral suppression it actually *achieved*, not the one it was
asked for — so an intervention that falls short simply lands further left, and
031's labelling problem disappears by construction.

- **x — achieved VLS among all PLHIV** (15+, at 2030), the outcome of a
  progressively more complete cascade.
- **y — PrEP coverage**, expanding outward from the highest-risk group.
- **z — two panels:** cumulative infections averted 2026-2040, and % reduction
  in HIV incidence from baseline.

The decision question this answers, which no arrangement of discrete arms can:
**for a level of suppression that is actually achievable, how much is adding
PrEP worth — and what PrEP coverage substitutes for what suppression gain?**

## Design

**Cascade ladder (x), nested and cumulative.** Each rung adds one lever, in
increasing order of difficulty:

| rung | lever added | how |
|---|---|---|
| C0 | — | baseline |
| C1 | linkage | ART coverage target → 0.95; fills only from the diagnosed, so this maximises on-ART-given-aware |
| C2 | suppression | + VLS-given-ART → 0.99 **by sex** |
| C3 | testing x2 | + routine testing rate doubled |
| C4 | testing x3 | |
| C5 | testing x4 | testing saturates here — the plateau probability hits 1.0/yr |

**Suppression is raised by sex only, not by age.** `VLSStockTarget` is sex-keyed
(`groupby('Gender')`, and it ranks agents with `for sex, is_sex in (('f',
ppl.female), ('m', ppl.male))`); it has no age dimension. Adding `AgeBin` rows
would give each sex several values per year and the lookup would silently
resolve to one of them. Age-stratifying it is a code change to the class exp 022
introduced and model-v1.3 depends on, and it perturbs a calibrated model — so it
is **deliberately out of scope here** and recorded in Next. 031 obs 10's caveat
stands: the model's suppression carries no age gradient and so understates the
youth cascade gap.

The VLS target must exceed the observed 0.959 (f) / 0.967 (m), or
`scenarios._ramp_long` raises — a deliberate guard against arms that silently do
nothing. 0.99 is used, roughly the best any stratum achieves in SHIMS3
(50+ reads 0.989 / 0.992), so it is a defensible ceiling rather than a round
number.

**PrEP ladder (y), nested and disjoint.**

| rung | groups covered |
|---|---|
| P0 | none |
| P1 | FSW @ 60% |
| P2 | + higher-risk AGYW @ 30% |
| P3 | + remaining AGYW @ 30% |
| P4 | + remaining women 25-34 @ 30% |

**Eligibility is made disjoint by construction.** `scenarios.agyw()` is
`female & 15-24` and does **not** exclude FSW, so the shipped helpers overlap.
Each rung here subtracts the groups already covered by lower rungs, and each
`sti.Prep` instance gets its own name — otherwise an agent could be enrolled
twice and the person-year denominator would be wrong, which would corrupt the
efficiency numbers that were 031's most robust output.

**Grid.** 6 cascade rungs x 5 PrEP rungs = **30 cells x 10 seeds = 300 sims**,
1985-2041. Common random numbers across cells, so every comparison is paired —
031's seed SD was 4,440 on a mean of 66,178, and paired contrasts were far
tighter than that implies.

**Compute.** Raccoon (120 cores), ~3 waves.

## What the surface is expected to show, and the trap in it

**The x axis is an outcome, not a control.** Two cells can land at the same
achieved VLS|PLHIV and avert different numbers of infections, because *who* is
suppressed matters — men 25-34 sit at 0.623 suppressed and carry far more
onward transmission per person than women 50+ at 0.939.

**Vertical scatter at equal x is therefore the finding, not noise.** It is 031
obs 10 reappearing on the x axis: the aggregate cascade metric does not
determine impact. **Plot the individual cells; do not fit a smooth surface over
them and report only the fit.**

Expect the axis to be strongly non-linear. Per 031, linkage is already
0.95-0.98 everywhere except men 25-34 (0.810), and the suppression deficit sits
in 15-24 who hold 4.4% of PLHIV — so C1 and C2 should move x very little, and
**testing should do nearly all the work.** If that is what happens, it is a
result: it says the first 95 is the only treatment-side lever left in Eswatini.

## Success criteria

**Good.** The surface is populated across a usable range of achieved VLS|PLHIV,
the PrEP ladder separates cleanly at each x, and iso-impact contours can be read
off to give a substitution rate between suppression and PrEP coverage.

**Informative failure.** The cascade axis barely moves even at C5 — i.e. the
model cannot be pushed much past its baseline suppression at any plausible
testing intensity. That would say the residual burden is structurally beyond
treatment, which is the strongest possible form of the "prevention still
matters" claim and is directly abstract-worthy.

**Stop condition.** C0/P0 must reproduce 031's baseline (66,178, seed range
61,863-76,693) within seed spread. Same model, same parameters, same seeds — a
mismatch means the grid harness differs from 031's in a way that has not been
accounted for.

## Notes

- Deadline-driven scope. The abstract goes to co-authors before the 1 Oct
  submission, so age-stratified suppression, a perfect-cascade bound, and a
  retention/long-acting-ART arm are all deferred to a later experiment.
- This subsumes the FSW-PrEP-on-cascade combination that 031 identified as
  missing: it is simply cell (C1..C5, P1).
- Parameters are 024 design row 868 throughout, as 026/029/031. One point, so
  every interval is Monte Carlo error across seeds and **not** a credible
  interval.
