# Exp 034 — A cascade ladder Eswatini could actually climb, against the PrEP gradient

**Question.** [033](../033_cascade_ceiling/SUMMARY.md) obs 9 established that its
own top rungs are not policy: they demand 0.999 awareness, 0.998 linkage and
0.999 suppression from *every* age-sex group, against SHIMS3's measured
**0.937 / 0.973 / 0.962** and a best-ever single stratum of 0.983 / 0.998 /
0.992. The model only reaches them because it has no never-testing subgroup, no
linkage delay in the stratified ART path, and no age gradient in suppression.

So: **what does prevention add across a cascade ladder built from targets
Eswatini could plausibly hit?**

**The framing that makes this experiment different.** Eswatini's problem is not
the aggregate — the national figures already clear the second and third 95s
(0.973 linkage, 0.962 suppression). **It is equity.** Men 25-34 sit at 0.792
aware and 0.810 linked while women 50+ are at 0.999 and 0.998. Every target
below is therefore applied **per stratum**, so `_ramp_long`'s per-stratum `max()`
raises the laggards and leaves the leaders untouched. "Achieving 95-95-95" here
means achieving it *in every group*, which is a real and unmet goal, rather than
on average, which is already done.

## The ladder

| rung | testing | ART target | VLS target | what it represents |
|---|---|---|---|---|
| **S0** status quo | x1 | — | — | current trajectory |
| **S1** testing only | **x3** | — | — | find everyone, expand nothing else |
| **S2** UNAIDS 95-95-95 | x3 | **0.9025** | — | 95 aware x 95 linked, in every stratum |
| **S3** 99-96-98 | x4 | 0.96 | 0.98 | near-universal testing; 96% of PLHIV on ART; 98% suppressed |

Each rung is reported with the cascade it actually achieves (aware / on-ART
given aware / suppressed given ART, 2030), because the target and the outcome
are not the same thing — S2 targets 95-95-95 per stratum and lands at 99-95-96
in aggregate, since raising the laggards pulls the average past the target.

**The theoretical bound has been removed.** It sat here as S4 (ART 0.999 / VLS
0.999) so the plausible rungs could be read against a limiting case. Exp 033
already establishes it, nobody could pursue it, and carrying it in a table of
scenarios invited it being read as one.

**S3 is named for what it achieves, not for a judgement about it.** It was
previously labelled "best-in-class", which asserted a comparison to other
programmes that this experiment does not make.

**Why testing x3 is the "first 95" rung.** Measured from 033's runs, per
stratum at 2030: at x1 five of eight strata sit below 0.95 awareness (worst
0.781); at x2, three (worst 0.887); at **x3, one** (women 15-24 at 0.933, the
rest ≥0.96). Tripling the routine rate — 0.75/yr testing probability among
undiagnosed adults against 0.25 now — is roughly the price of the first 95
everywhere. Ambitious but not counterfactual: it is universal annual testing
with good uptake, plus index and self-testing.

**Why S1 exists.** 033 obs 3 found that with the ART target fixed, scaling
testing cuts the undiagnosed by 1.18 pp while raising diagnosed-but-untreated by
0.94 pp — ART coverage stays pinned and suppression barely moves. S1 isolates
that: **diagnosis with no expansion of treatment.** If it delivers little, that
is the headline programme lesson, cleanly demonstrated rather than inferred.

**Why the ART target is 0.9025 in S2.** The model's ART input is coverage among
*all PLHIV*, not among the diagnosed. UNAIDS 95-95 is 0.95 x 0.95 = **0.9025 of
PLHIV on ART** — the honest translation. Setting it per stratum raises men 25-34
(currently ~0.64 of PLHIV on ART) and leaves women 50+ alone.

**Why the third 95 is left at baseline in S1 and S2.** Eswatini's measured
suppression-given-ART is 0.962, already above 95%. Holding it at baseline is the
honest statement that the third 95 is met; raising it would credit the scenario
with a gain that is not available. It rises only at S3, and even there to 0.98 —
below the 0.992 the best stratum achieves. Note the model cannot represent the
youth suppression deficit that would be the real target (029 obs 6), so S3's
third-95 gain is likely optimistic.

**S4 is not a scenario and must never be reported as one.** It is carried so the
plausible rungs can be read against the limiting case.

**PrEP axis** — unchanged from 032/033, nested and disjoint: none → FSW 60% →
+ higher-risk AGYW 30% → + remaining AGYW 30% → + remaining women 25-34 30%.

**Grid.** 5 x 5 = 25 cells x 10 seeds = **250 sims** on raccoon.

## Success criteria

**Good.** The plausible rungs (S2, S3) land in the 0.91-0.96 suppression range,
PrEP's increment is resolvable at each, and the surface supports a substitution
statement in the zone where policy actually operates.

**The informative result.** S1 delivering little — confirming that testing
without treatment expansion relocates people within the cascade rather than
moving them through it.

**Stop condition.** S0/P0 must reproduce 66,178 within seed spread, as 032 and
033 both did exactly.

## Notes

- Deadline-driven: the abstract goes to co-authors before 1 Oct. Age-targeted
  testing, a never-testing subgroup, and age-stratified suppression are all
  deferred — each is a model change, and all three would *lower* the achievable
  ceiling rather than raise it, so their absence makes these rungs optimistic.
- One parameter point (024 row 868). Intervals are seed noise, never credible
  intervals.
