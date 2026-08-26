# Review of stisim v1.5.11 (PR #580, `rc1.5.11 -> main`)

Working copy of the comment posted to
<https://github.com/starsimhub/stisim/pull/580>. Kept here because it motivates
our adoption decision (a future experiment) and records what we asked for.

**Context at time of writing (2026-08-24):** 1.5.11 not yet released (latest tag
v1.5.10). Release PR #580 open. Contents: #556 (ART adherence), #561 (ART
mortality simplification), #576 (traditional circumcision). Our model runs
stisim 1.5.8 / starsim 3.5.0, but `pyproject.toml` on `main` now pins
`stisim>=1.5.10, starsim>=3.5.2` (from Robyn's PR #2), so declared and actual
already disagree.

---

Reviewing as a downstream user — we run an Eswatini HIV model and are about to
calibrate against this release, so these are from the perspective of someone
who has to live with the defaults.

The ART adherence work is a real improvement for us: on-ART mortality and
non-suppressive ART were the two things most limiting what we could say about
the treatment cascade. Anchoring on-ART mortality to the off-ART CD4 hazard, so
"ART is never worse than no ART" holds by construction rather than by parameter
luck, is a nice piece of design.

Comments below, roughly in order of how much I'd want them addressed *before*
the release is cut.

## 1. Sequencing: #577 isn't in this release

1.5.11 ships the adherence machinery — effective vs non-suppressive states,
adherence-dependent mortality and transmission — but `p_effective_art` is a
scalar (`ss.bernoulli(p=1.0)`). PR #577 (`art_vls_ingest`), which lets
`vls_coverage` take the same time-varying and age/sex-stratified formats as
`coverage`, is still open against `rc1.5.11`.

If #580 merges first, 1.5.11 has adherence-dependent *mortality* without the
means to drive suppression from data — which for us is the main reason to adopt
it. Viral suppression among people on ART in Eswatini rose substantially between
the 2011 and 2021 PHIA rounds (~95% by SHIMS3); a fixed scalar can't represent
programme improvement, which is exactly what a treatment-cascade scenario is
about.

Any chance of landing #577 before cutting 1.5.11? If not, worth saying in the
release notes that suppression is scalar-only until #577 lands.

## 2. The 1.5.11 CHANGELOG is empty

`## Version 1.5.11 (TBC)` currently has no entries, but the release contains at
least three things that will surprise people upgrading:

- **On-ART mortality is nonzero by default** (`rel_art_mortality_effective =
  0.25`). Previously agents on ART had zero HIV-specific mortality. Not opt-in —
  which is presumably why `baseline.yaml` was regenerated.
- **Circumcision state and efficacy moved from `VMMC` onto `HIV`** (#576).
  `self.circumcised` / `self.pars.eff_circ` on the intervention are gone,
  replaced by `HIV.circumcise()` and `HIVPars.eff_circ`. Any `VMMC` subclass
  touching those breaks. (Ours did — easy to fix once spotted, silent if not.)
- **`never_art` -> `art_naive`.**

Models calibrated against 1.5.8 will need recalibration. The 1.5.9 notes said
exactly that for the VMMC change; it applies here at least as strongly.

## 3. ART mortality — what did the simplification preserve?

#556 says the mortality follows `ARTMortalityTable` from EMOD-HIV, and #561
simplified it to an age/CD4 function. Could the notes or a docstring record
where the simplification diverges from EMOD? Specifically, which of
`rel_art_mortality_effective = 0.25`, `unsupp_m = 0.7`, `unsupp_f = 0.35`,
`rel_death_f = 0.74`, `art_death_age = [1.0, 1.10, 1.21, 1.32]` are
EMOD-derived and which are new?

The commit introducing the 2x male:female non-suppressive ratio says the number
"is not based on anything empirical, and should be revisited later" — which is
fine as a starting point, but downstream that distinction matters a lot. These
are precisely the parameters we'd want either to fix on evidence or open into a
calibration prior, and we can't tell which is which from the code.

## 4. The on-ART <= off-ART invariant is enforced by comment, not code

The guarantee currently rests on a docstring computing `0.7 x 1.32 = 0.924 <= 1`
and warning users not to break it. Anyone calibrating
`rel_art_mortality_unsupp_m` above `1/1.32 ~= 0.758` silently inverts the
treatment effect. Since making that impossible is the whole point of the
anchoring design, an assertion at `init_pre` seems worth it:

```python
max_mult = max(m for *_, m in self.pars.art_death_age)
worst = max(self.pars.rel_art_mortality_effective,
            self.pars.rel_art_mortality_unsupp_m,
            self.pars.rel_art_mortality_unsupp_f) * max_mult
assert worst <= 1, 'on-ART mortality can exceed off-ART at the same CD4'
```

Related: there's no test asserting the invariant, or the sex/adherence ratios.
The 1.5.9 notes observe that `test_vmmc_specs` only checking "some
circumcisions occurred" is why the VMMC overshoot went undetected — the same
gap exists here, on this design's central claim.

## 5. `rel_death_f` applies only on ART

`rate[~male] *= self.pars.rel_death_f` sits inside `get_art_mortality_hazard`,
so women get a 26% mortality advantage on treatment and none off it —
discontinuation removes it discontinuously. Women have better HIV survival in
both states. Is the on-ART-only application deliberate?

## 6. Transmission: is `nonsupp_art_efficacy = 0.35` derivable?

No citation, and the right value depends entirely on what "non-suppressive"
means — VL > 1000, > 200 or > 50 imply very different transmissibility.

One way to anchor it: take a VL distribution among the non-suppressed and apply
Quinn (Rakai 2000, ~2.45x per log10) relative to untreated chronic infection at
~4.5 log10. With log10 VL ~ N(3.75, 0.8^2):

```
2.45^(3.75 - 4.5) x exp(ln(2.45)^2 x 0.8^2 / 2) = 0.51 x 1.29 ~= 0.66
-> efficacy ~= 0.34
```

which lands almost exactly on the shipped 0.35. (The second term is the Jensen
correction — the hazard is exponential in log VL, so integrating over the
distribution differs materially from evaluating at the mean.)

So the value looks defensible. It would help to state the assumed VL
distribution, because the answer is sensitive: mean 4.0 gives ~0.16, mean 3.5
gives ~0.47.

## 7. `effective_art_efficacy = 0.96` vs the docs' own U=U claim

The HIV user guide says "undetectable = untransmissible", then the parameter
retains 4% of transmission. PARTNER, PARTNER2 and Opposites Attract found zero
linked transmissions from virally suppressed partners across thousands of
couple-years. Suggest ~0.99-1.0.

This is pre-existing rather than introduced here, but #556 makes it more
consequential: suppression status now drives both transmission and mortality, so
the model will increasingly be used to evaluate "improve viral suppression"
against alternatives. A 4% residual systematically understates the benefit of
suppression, which biases that comparison in a specific direction.

## 8. The efficacy ramp

`rel_trans` ramps linearly to full efficacy over `time_to_art_efficacy`
(6 months). Real viral load falls roughly 1 log10/month after initiation, with
most people suppressed by ~3 months — so a linear ramp puts an agent at ~48% of
full efficacy at month 3, where reality is near-fully suppressed. An exponential
decay would track VL kinetics better. This compounds #7 in the same direction.

Also: should non-suppressive agents ramp at all? They currently ramp linearly to
0.35 over the same 6 months.

Minor/latent: the ramp mixes units. `new_on_art` uses `time_to_full_eff/self.dt`
(timesteps) but the interpolation divides by `time_to_full_eff.value` (6, in
months). These agree only when `dt` is exactly 1 month — which is the HIV
module's default, so it's currently fine. At weekly resolution
`efficacy_to_date` exceeds 1 and `rel_trans` goes negative.

---

Happy to open issues for any of these rather than hold up the release — #1 and
#2 are the ones I'd most want resolved before cutting, and the rest are fine as
follow-ups.
