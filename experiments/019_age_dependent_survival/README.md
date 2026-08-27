# Exp 019 — Age-dependent untreated survival: is that the missing AIDS mortality?

> **Pre-registered ahead of its turn.** [018](../018_adopt_and_size/README.md)
> runs first and is not yet closed. No code runs in this folder until it is.
> Written now because the design question was live; the seed count below is
> deliberately left to 018's replicate measurement.

## Question

[016](../016_double_counted_mortality/SUMMARY.md) found the HIV module supplies
~65% of peak UNAIDS AIDS deaths and ~50% post-2010.
[017](../017_version_bump/SUMMARY.md) eliminated two candidate explanations —
on-ART immortality (the bump moved peak deaths less than seed noise) and
`rel_death_f` (arm D, |z| ≤ 0.9). After both the version bump and 016's
mortality correction, arm C reached **64% of the UNAIDS peak**: ~4 000
deaths/year still missing.

017 observation 7 located where the remaining candidates must live. The model
kills HIV-positive agents two ways, and they are wildly unequal:

| route | mechanism | tunable by | share of untreated deaths |
|---|---|---|---|
| `p_hiv_death` | per-timestep Bernoulli, CD4-stratified (0.003/y at CD4 ≥ 500 rising to 0.30/y below 50) | `rel_death`, `rel_death_f`, on-ART multipliers | **~15–20%** (estimated) |
| `ti_zero` | death date drawn at infection: `dur_acute + dur_latent + dur_falling` | **nothing** | **~80–85%** (estimated) |

**Every mortality knob in stisim acts on the minority route.** That is why 017's
two candidates were both inert — they were tuning ~20% of deaths. The route
that decides the rest has no age dependence, no sex dependence and no
multiplier of any kind: untreated survival is `lognorm_ex(10 y, 3 y)` latency
plus `lognorm_ex(3 y, 1 y)` late-stage, mean **13.1 y**, applied identically to
a 17-year-old and a 55-year-old.

Real untreated survival falls steeply with age at seroconversion — ALPHA-network
pooled estimates run ~12.5 y at 15–24 down to ~7.5 y at 45+. **Two things are
confounded in that comparison and this experiment separates them:**

1. **Level.** Is 13.1 y simply too long? Cohort medians sit nearer 11–12 y.
2. **Gradient.** Is the *absence* of age dependence the problem, independent of
   the mean?

## Plan

**Design: four arms at fixed parameters and seeds**, the lineage of 015/016/017.
The only thing varying is a multiplier on `dur_latent` as a function of age at
infection.

| arm | multiplier by age at infection | mean survival | isolates |
|---|---|---|---|
| **A** `flat_13` | 1.0 everywhere (stisim default) | 13.1 y | baseline |
| **B** `flat_11.5` | 0.84 everywhere | 11.5 y | **level** (A→B) |
| **C** `grad_mild` | 0.89 / 0.84 / 0.73 / 0.61 | ~11.5 y | gradient, half strength |
| **D** `grad_alpha` | 0.94 / 0.84 / 0.64 / 0.44 | ~11.5 y | **gradient** (B→D), full ALPHA |

Age bands: 15–24, 25–34, 35–44, 45+.

**The design point.** A flat 0.84 gives ~11.5 y for everyone, which is
approximately the *population mean* of the ALPHA gradient once weighted by the
model's age-at-infection distribution. So B, C and D sit at roughly the same
level and differ only in shape — **B→D is the gradient effect at constant
level**, which is the contrast that matters. A→B is the level effect at zero
gradient. Realized mean survival is measured per arm (metric 3) rather than
assumed, since the age distribution of incident infections is endogenous and
will itself shift between arms.

**Implementation.** An in-repo `HIV` subclass overriding `set_prognoses()`:
call `super()`, then rescale the latent interval by age while preserving the
drawn `dur_falling`, so only the latency stretches or compresses. Injected via
`make_sim(hiv_class=...)`, which [run_sims.py:24](../../run_sims.py#L24)
already supports — the same pattern as `vmmc_class`. **In-repo, not a patch to
the editable stisim checkout**, so a `git pull` cannot silently wipe it.

**Parameters and seeds.** High-transmission set only. 017 established the
default-parameter arm cannot carry quantitative claims (CV ≈ 47% vs ≈ 4%), and
this experiment is entirely about effect sizes. **Seed count comes from 018's
replicate measurement** — placeholder 10, to be set before running.

## Metrics

1. **Deaths by route — `ti_zero` vs `p_hiv_death`.** The first measurement, and
   it gates the interpretation of everything else. The ~80/20 split above is
   integrated from the hazard table, not measured; if it is wrong, the premise
   of the experiment changes. Needs an analyzer counting each route separately,
   since `hiv.new_deaths` pools them.
2. **AIDS deaths vs UNAIDS** — total trajectory and peak-to-peak. The target.
   Baseline is 017 arm C at 64% of the UNAIDS peak.
3. **Realized survival distribution per arm** — mean and median time from
   `ti_infected` to death, overall and by age at infection. Validates the
   subclass does what it claims and quantifies the level/gradient separation.
4. **Age distribution of AIDS deaths**, against
   `../016_double_counted_mortality/outputs/excess_deaths_by_age_sex.csv` —
   016's implied AIDS deaths by age band, reconstructed from all-cause mortality
   with no HIV information as input. An independent age-stratified target the
   model has never been compared against.
5. **HIV prevalence 15–49 and by age/sex vs PHIA.** Shortening survival removes
   PLHIV, so prevalence will *fall* — against the direction 014's coverage gap
   needs. The size of that trade-off is the experiment's main cost, and may be
   what decides whether this is adoptable.
6. **Population**, as a demographic sanity check.

## Success criteria

- **Level explains it:** arm B closes a substantial share of the deficit and the
  gradient adds little. Simple fix — `dur_latent` becomes a calibration
  parameter, no subclass needed, no upstream conversation.
- **Gradient explains it:** B→D closes the gap where B alone does not, and
  improves the age distribution of deaths (metric 4). This is the interesting
  outcome, and it is **a structural gap in stisim, not a parameter we can
  tune** — it becomes an upstream feature request with evidence attached.
- **Neither explains it:** survival duration is not the answer and the deficit
  lies elsewhere — most likely in the acute/falling transmission multipliers,
  ART cascade timing, or the targets themselves. A clean null redirects 020.
- **The trade-off bites:** an arm closes the death gap but drops prevalence
  below the PHIA targets. That is a real finding about model structure — the
  model may not be able to match deaths and prevalence simultaneously with the
  current natural history, which is exactly the kind of thing to know *before*
  a calibration tries to split the difference.

## Why this matters beyond the deficit

The decision question in [CLAUDE.md](../../CLAUDE.md) turns on the value of
raising viral suppression against scaling prevention. **Survival on ART versus
survival off it is the quantity that comparison is made of.** If untreated
survival is 1.5 years too long and flat in age where it should be steep, the
model systematically understates what treatment averts — and it does so
differentially by age, which is precisely the axis a targeting decision runs on.
This is not a tidy-up; it is upstream of the answer.

## Not in scope

- The `p_hiv_death` hazard table itself. 014 measured ρ = 0.11 for `mort_mult`,
  and metric 1 will confirm it is the minor route. If metric 1 surprises us,
  that changes.
- Sex differences in survival. `rel_death_f` was refuted by 017 arm D, and it
  acts on the minor route anyway.
- On-ART survival. 017 observation 8 established on-ART mortality is inert by
  construction, and untreated survival is the prior question.
- Any calibration run, or re-opening prior bounds — `parameter-engineering`
  after this.

## Open question for the design

Whether to scale total survival (acute + latent + falling together) rather than
latency alone. Scaling everything is arguably more faithful — late-stage
progression also accelerates with age — but `dur_falling` sets the CD4 decline
rate, which feeds `rel_trans_falling` (8× transmissibility in late infection),
so compressing it changes transmission as well as mortality. Latency-only keeps
the arms interpretable. Revisit if the gradient turns out to matter.
