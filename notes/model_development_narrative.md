# From the first experiment to wave 1

*Experiments 001–025, April to September 2026. Written 2026-09-07 from
`experiments/log.md` and the individual SUMMARYs.*

Part one covers the model-development arc, 001–009, and the point at which it
turned into a calibration. Part two covers 011–025, the four months between the
first failed coverage check and the first history-matching wave.

**Two things about the numbering.** There is no experiment 010 — it was planned
in 009's SUMMARY as a diagnostic prior expansion, deferred when the network bug
surfaced, and eventually absorbed into 012. And 015 ran a month *before* 014
(10 July against 11 August): 013 identified a broken VMMC that had to be fixed
before the coverage check could mean anything, and the folder was numbered when
it was queued rather than when it ran.

---

## The starting point

In April the Eswatini model ran end to end and produced a plausible-looking
epidemic, but almost nothing about its inputs had been checked against Eswatini
specifically. ART coverage came in as a national total. Sexual debut used the
stisim package defaults. There was no circumcision. The acute-phase
transmission parameters were the package's, which were known to be high. The
seven experiments that followed were each a local fix to one of those, and each
one is defensible on its own terms. What they added up to is the more
interesting story.

## 001–002 — ART coverage, and a model that could not reproduce its own input

The first change swapped the national ART total for PHIA's coverage
stratified by age band, sex and year. The motivation was straightforward: the
unstratified input applied one coverage to everybody, while PHIA shows very
large differences — 90% among women 35–45 against 36% among men 15–25 in 2016.

It did not work. The model absorbed the stratified input and produced male and
female ART curves that nearly overlapped. **The lesson was that a model can
fail to reproduce its own inputs, and that this is worth checking explicitly.**
Reading the stisim source found the cause: the allocation summed every
stratum's target into a single national number and then filled it by CD4
priority across the whole population, so older and sicker agents were treated
at the expense of younger ones and the sex differential washed out entirely.

Fixing it to allocate within each age × sex stratum produced the largest single
improvement of the whole period. The bug was filed upstream.

## 003 — Sexual debut, and the first sign that one lever would not be enough

Package defaults put debut at 20 (women) and 21 (men). DHS Eswatini puts the
median nearer 17–18 and 18–19; the EMOD model of the same setting used 16.3 and
17.5. Late debut delays network entry by three or four years and was suppressing
prevalence in young women, so debut moved to 17.5 and 18.5.

A follow-up tightened the distribution (SD 2.5 → 1, because the wide version
overlapped the old setting so heavily that the change was invisible) and swept
debut age across 13–20 at ten seeds each. **The sweep gave the first real
dose-response curve of the project: roughly 0.3–0.5 additional infections per
100 person-years at the epidemic peak for each year earlier.** Clear, monotone,
and comparable for both sexes.

It was also the first clear negative result. Female incidence still came in at
0.7 against a target of 1.7 for 2016, and the roughly 2× female-to-male ratio in
the survey data was not reproduced. **Debut age alone could not close the gap,
which pointed at the transmission parameters.**

## 004 — Three levers at once, and a caveat that came back later

Experiment 004 attacked the female-to-male gap three ways. A sweep of the
female-to-male transmission multiplier showed 0.25 reproduces the observed ~2×
ratio, but that lowering it also drops absolute incidence, because fewer
infected men means fewer onward infections to women — **the ratio and the level
are not independently controllable by that parameter.** Scaling condom use to
half the stisim default brought 2016 incidence closest to the surveys, on the
argument that "used at last sex" overstates consistent use. And a 1.7×
susceptibility multiplier for women 15–24 was cherry-picked from an unmerged
stisim PR.

All three were adopted, and the experiment recorded its own caveat plainly:
halving the condom input without re-fitting left the transmission probability
and condom efficacy calibrated against inputs that no longer existed. **That
caveat is the thread that runs through everything after it.**

## 005 — Circumcision, and stock versus flow

VMMC was added with coverage read as a *prevalence target* — the share of men in
an age band who are circumcised — rather than as a per-step hazard. This matters
because circumcision is irreversible and the surveys measure a cross-sectional
stock. The data were triangulated across SDHS 2007, SHIMS2 2016 and SHIMS3 2021
(Table 12.5: 47.2% medical plus 1.1% traditional for 15–49).

Two lessons. First, **an intervention's data semantics have to match its
mechanism**; this same distinction reappeared later when an upstream rewrite
reverted it. Second, interpolation outside the data range is not free: with only
2007/2016/2021 anchors the smoother propagated the 2007 values back to 1985,
implying 20% circumcision in the 1980s, and a 1990 anchor had to be added.

The cost was that male incidence, already low, dropped further — to about 0.3
against a target of 0.85.

## 006 — Acute-phase transmission, and paying for credibility

Bellan 2015's central estimates replaced the package defaults: acute duration
3 months → 1.7, relative transmissibility 6 → 5.3, which takes the excess hazard
attributable to acute infection from roughly 15 months down to 7.3. Incidence
fell 25–30% across the board.

**This is the clearest case in the period of a change that made the model more
credible and the fit worse**, and it was adopted anyway, on the grounds that the
natural history should be right and the transmission parameters were the things
that ought to absorb it.

## 007 — Seeding, and the difference between noise and level

Snapshotting the model at initialisation found only about 11 infected agents out
of ~4,440 active adults in 1985, with most strata expecting fewer than one. The
early epidemic was therefore a Bernoulli lottery. Doubling the seed prevalence
took it to about 20.

**The finding worth keeping is what it did not do.** Early-epidemic uncertainty
bands narrowed, as expected, but the 2016 and 2021 levels were untouched — after
thirty years the initial seed count is irrelevant. Seeding controls early
certainty, not late level. (A later experiment complicated this: when the prior
was allowed to sample *below* the establishment threshold, seeding became the
switch determining whether an epidemic happened at all, and briefly looked like
one of the most important parameters in the model. It is not.)

## The transition — 008 and 009

By 007 the model was better in seven specific, well-documented ways and its
transmission parameters had been fit against a version of the model that no
longer existed. Three of the seven experiments recorded "re-calibration needed"
as their decision, and each time it was deferred to the next one.

**008 stopped adding changes and wrote down the target set instead:** 89 rows,
54 of PHIA prevalence by age × sex × year and 35 of UNAIDS AIDS deaths. Doing
that surfaced things that informal tuning had not — that the deaths series had
no published uncertainty and was carrying a ±15% placeholder, that PHIA 2011
only covers 15–49 while the other rounds reach 65, and that a decision was
needed on whether 2021 was a target or a hold-out.

**009 then ran the first prior predictive check** — 50 draws from the prior,
asking not "does the model fit?" but the prior question, "can the model produce
this data at all?" It came back at 30 of 89 rows inside the 5–95% envelope, 20%
for deaths and 43% for prevalence, and both missed systematically rather than at
the edges. That is a stop-and-diagnose result, and it ended the practice of
tuning one parameter at a time.

## What the first arc taught

Every one of the seven changes was individually justified and most were
improvements. But each one moved transmission, none was followed by a re-fit,
and the deferral compounded: by 007 the model's most important parameters were
inherited from a configuration that had been superseded five times. **A sequence
of locally correct changes does not compose into a calibrated model.**

The second lesson is about the value of writing the targets down. 008 produced
no new simulation output at all, and it was still one of the most useful
experiments of the period, because a target set with explicit uncertainties is
what makes it possible to say a model has failed.

One honest coda. 009's diagnosis — that the model produced too little mortality
flow — turned out to be wrong; experiment 014 opened those exact parameters and
measured them to be nearly inert, and 024 eventually showed the residual bias
was a parameter-value problem rather than a structural one. **The check was
right to stop us. The explanation attached to it was not, and it took five more
experiments to find that out.**

---

# Part two — 011 to 025, the calibration proper

Most of the next fourteen experiments produced no improvement in fit. That is
not a complaint about them; it is what the period was for. Several returned
something more useful than an improvement, which was a correct diagnosis.

## 011–013 — The network detour, and an attribution that was wrong twice

009's failure had several suspects, and one was the sexual network. **011 asked
a narrower, blocking question: does the age-gap parameter actually control the
realised age gap?** It did not. Partners were matched by rank within the
surviving pool, which flattens the configured distribution into a function of
the marginal age supplies. Sweeping the configured mean 0 → 7 → 14 years moved
the realised mean only 0.7 → 1.4 → 5.3, and women 35–49 came out with a
biologically impossible **−2.7-year** gap. **A prior over a parameter that does
not transmit is meaningless**, so calibrating the network was off the table
until this was fixed. Upstream PR #477 fixed it; our own local patch turned out
to flatten identically to the broken version and was discarded.

**012 re-baselined on the upgraded stack** and found peak prevalence had moved
*down*, which flipped the direction of the problem: 009 had diagnosed prevalence
drifting too high, and now it undershot. The experiment attributed the drop to
the network fix.

**013 tested that attribution and refuted it.** Holding everything else fixed
and flipping only the matcher, the corrected version *raises* 2021 prevalence
(4.09% → 6.46%); it does not lower it. So the drop came from somewhere else in
the upgrade. Inspection found the real culprit: the stisim 1.5.6 rewrite of the
VMMC class had silently overwritten our exp-005 patch, and circumcision was
ratcheting to ~100%. **The lesson was about attribution discipline — 012 had
assigned a cause without an experiment that could separate it from the
alternatives, and 013 cost three weeks establishing that the cause was wrong.**

## 015 — Fixing VMMC, and learning where in-repo fixes belong

Prevalence-target circumcision was re-implemented as an **in-repo subclass
rather than a patch to the editable stisim checkout** — because that is exactly
how the 005 version had been lost. Coverage went from a 99.3% overshoot to 45.4%
against SHIMS3's 48.3%, and the epidemiological effect was large: prevalence
15–49 rose 6.69% → 11.94%, with male prevalence roughly doubling.
Over-circumcision had been suppressing male acquisition hard.

## 014 — The second coverage check, and a diagnosis inverted

With the network and VMMC both corrected, the check was re-run with the three
suspects 009 had named opened into the prior. **It came back at 4 of 89 rows —
worse than the 34% it was trying to beat.**

The raw number was misleading and the experiment said so. 85 of 89 rows had the
*entire* ensemble below the observation and not one above it, but the best
single draw reached 0.502 against an observed 0.507. **The model could reach the
data; the prior was in the wrong place.** Two things had gone wrong: the prior
now sampled below the epidemic establishment threshold, so 22% of draws were
dead epidemics dragging the envelope to the floor, and it spread 50 draws across
nine dimensions instead of six.

The scientific return was the inversion. The mortality multiplier that 009 had
named as its leading suspect measured **ρ = −0.01 over a six-fold range** — it
does essentially nothing, because AIDS deaths fire from a separate pathway it
does not touch. **009's explanation was dead, and the best draws were piled
against the transmission parameter's upper bound with the data on the far side
of it.** The problem was transmission and seeding, not mortality.

**Also learned, and worth keeping: widening a prior can lower coverage.** Added
volume in directions that do not matter dilutes the draws that land where the
data lives.

## 016, 017, 019, 022 — Four experiments chasing the AIDS-death deficit

**016** found the background mortality rates fed to the demographic module
already contained AIDS deaths — their adult rates rise 5.5× to a 2005 peak,
which is the shape of the epidemic. Reconstructing a non-AIDS counterfactual
from all-cause data alone, using no HIV information as input, reproduced the
UNAIDS AIDS-death curve to within about 10%, and correcting it took the 2015
population from 13.8% below target to 0.1% above. Real, and adopted. But it
bought only about 1.5 percentage points of prevalence where roughly 6 were
needed, **and in the process exposed the actual problem: the HIV module supplies
only ~65% of peak AIDS deaths.**

**017** tested the most attractive explanation — that nobody on ART could die of
HIV — by taking the stisim release that made on-ART mortality nonzero. It closed
none of the deficit. The reason is a design invariant: the on-ART hazard is
anchored to the off-ART hazard *at the same CD4*, and agents on ART have
restored CD4, where that hazard is tiny. **The property that makes the design
safe is what makes it inert.** The same experiment refuted a second candidate,
a new female mortality multiplier, by setting it back to 1.0 and measuring
nothing. It also confirmed the whole version bump was behaviourally neutral,
which is the good outcome for an adoption decision.

**019** tested untreated survival, which was flat at 13.1 years for a 17-year-old
and a 55-year-old alike. Shortening it closes 39% of the deficit — and makes the
prevalence fit monotonically worse, from 0.0590 to 0.0782. **The two targets are
in structural tension.** It also measured, rather than estimated, that 70–74% of
deaths flow through a route no mortality parameter touches, which explains why
every knob tried so far had been inert.

**022** came at the same thing from the opposite direction, with EMOD's
*pivoting* gradient — longer survival at young ages, shorter at old. It landed on
the same frontier and destroyed the 45–64 fit doing it. **Three independent
routes to the same trade-off is not a coincidence; it is a structural
property**, and it is why deaths were subsequently down-weighted rather than
fitted hard. About 22% of the deficit survives every mechanism tested.

## 018, 020, 021 — Housekeeping that turned out to matter

**018** adopted the corrected stack as `model-v1.1` and removed an inherited PrEP
default nobody had chosen: a fabricated ramp to 80% of sex workers starting in
2004, a decade before efficacy evidence, which every experiment from 001 onward
had been carrying. It also mapped where the epidemic establishes at all.

**020** ran the population-size sweep that had been deferred three times, and
**reframed the question it was asked**. Going from 10,000 to 50,000 agents still
leaves 2 of 54 target strata below 5 expected cases. The stratum that will not
clear is young men in 2007, where the model puts prevalence at 0.003 against a
measured 0.019. **No amount of population fixes a fit error** — past 20,000
agents, more agents mostly buy precision on a bias.

**021** supplied viral suppression from PHIA instead of the default of 1.0, which
had every ART initiator suppressed. Epidemiologically null. But the model had
been overstating its own cascade by 8.8 percentage points, and **a counterfactual
that asks "what if we raised suppression?" measured against a baseline that
already overstates it understates the value of the treatment arm** — a bias in
the decision, invisible in the fit, that no calibration would have caught.

## 023 — Choosing what to open

400 draws, measuring which parameters actually move which targets. Nine
candidates pruned to seven, both cuts on evidence: seed prevalence measured
**ρ = 0.10**, the weakest in the set, and a concurrency parameter was both weak
and redundant with a stronger one.

**The seed-prevalence result reversed 014's own ranking**, which had put it
second of nine. The explanation is clean and worth remembering: 014's prior
sampled below the establishment threshold, so the parameter was acting as an
on/off switch for whether an epidemic happened at all, and that is what its
correlation was measuring. **A parameter can look important because it controls
whether the model works, not because it controls the outcome.**

The experiment also caught a flaw in its own design. The orthogonality check as
written correlated the prior draws, which are an independent uniform sample, and
so confirmed only that the sampler works. The question that matters — do two
parameters do the same *thing* to the model — needed correlating their effect
signatures across targets, and that found two pairs genuinely confounded.

## 024 — Wave 1, and the reversal

1000 design points in 456 seconds on 120 cores. The emulator came in at
R² = 0.923. **The coverage question that had blocked 009 and 014 was answered
affirmatively: 45 of 48 targets inside the 5–95% envelope, and one draw of 1000
within 3σ on all 48 at once** — joint reachability, which covering each target
separately does not imply.

Then the headline. Experiments 016 through 022 had all measured a −4.3 percentage
point prevalence deficit and converged on calling it structural misspecification;
wave 1 budgeted a discrepancy allowance to absorb it. **With that allowance set
to zero — survey sampling error alone — the best draw still sat at 3.98σ. The
deficit was not structural. It was a parameter-value deficit, and opening
transmission largely removed it.** Seven experiments had diagnosed a
model-structure problem while all holding the same parameters fixed.

Two smaller lessons. The 87% NROY looked uninformative and was not: one
observation constrains roughly one direction, and the cut was diagonal in the
box, so **reporting per-parameter intervals would have made the wave look like it
did nothing.** And the one defect the wave reported turned out to be in the
target construction rather than the model — a top age band built as 45–65 and
compared against survey data that stops at 49.

## 025 — Where it stands

Wave 2 is written, committed and smoke-tested: half the discrepancy allowance,
incidence and the female young:old ratio added as emulated features, 20,000
agents. **It is gated on a synthetic-data recovery test** — recovering known
parameters before touching real data — which 023 flagged as the largest
remaining gap in the sequence and which two experiments deferred. Blocked on
compute.

## What the second arc taught

**The most valuable experiments were the ones that returned nothing.** Six
explanations for the AIDS-death deficit were closed, and closing them is why the
remaining candidate can now be stated precisely rather than guessed at.

**Diagnosis is harder than measurement, and it failed repeatedly here.** 009
blamed mortality; 014 showed mortality was inert. 012 blamed the network; 013
showed the network moved prevalence the other way. 016–022 blamed model
structure; 024 showed it was a parameter value. In each case the *measurement*
was sound and the *explanation attached to it* was wrong. The pattern is that an
explanation formed while several things are held fixed tends to indict whatever
is varying.

**Corrections were caught by re-reading committed records, not by remembering.**
The sex mapping inverted in 018 and caught in 019; the age-band mismatch in 024;
the binning bug in 016 that briefly reversed that experiment's recommendation.
Every one surfaced because a later experiment recomputed something an earlier one
had reported — which only works if the earlier one wrote down what it did.

---

## Beyond this write-up

This account stops at 025 because that is where the calibration itself stops.
Work has since moved on to **026 — decision scenarios**, which is what the
calibration was for: long-acting PrEP against improvements to the treatment
cascade, aimed at a CROI abstract.

Its opening finding is worth knowing even if you read nothing else of it.
Eswatini has effectively already achieved the third 95 — suppression among
people on treatment sits at 0.959 for women and 0.967 for men. The remaining
headroom in the cascade is almost entirely in **coverage**, not suppression,
and it is concentrated in men 25–34, at 0.650 against a 0.95 target. That
reframes the comparison the model was built to make: "improve the cascade"
turns out to mean "find and treat men", not "suppress people better".

025 itself is still open at the time of writing — wave 2 is written and
smoke-tested but waiting on compute — so the scenario work in 026 runs on the
wave-1 configuration rather than on a converged calibration. Read its results
with that in mind.
