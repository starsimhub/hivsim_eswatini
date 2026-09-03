# Exp 023 — Parameter engineering: nine parameters, and a target set that includes incidence

> First experiment of the calibration proper. 022 closed the model-development
> sequence; this settles what is free, what is fixed, and what the calibration
> is fitted to. Coverage check v3 is 024.

## Question

**Which parameters should wave 1 open, on what priors, against which targets —
and does the prior actually move those targets?**

Seven experiments (016–022) failed to move the model's age-shape defect in
prevalence. The residual is now well characterised:

| band | women | men |
|---|---|---|
| 15–24 | −8.4 pp (−49%) | −3.7 pp (−65%) |
| 25–34 | −1.0 pp | **−12.7 pp (−47%)** |
| 35–44 | **+4.5 pp (+11%)** | −7.2 pp |
| 45–64 | −0.5 pp | −4.4 pp |

## The parameter set: four free, and why each earns its place

Per-act acquisition risk in this model is a product:

```
women:  beta_m2f × s_f(band)
men:    beta_m2f × rel_beta_f2m × s_m(band)
```

**Identifiability first.** If `beta_m2f` and the `s` multipliers are all free,
the model is unidentified: halving `beta_m2f` and doubling every `s_f` gives an
identical model, and `rel_beta_f2m` is only identified as a product with the
male multipliers. Two redundant dimensions, and history matching would wander
along that ridge reporting an NROY space wider than the data warrants. So the
25–34 band is **anchored at 1.0 for both sexes** — a choice of reference point,
not a biological claim, exactly as fixing sea level at zero says nothing about
mountains.

### Transmission and susceptibility

| parameter | prior | scale | mechanism |
|---|---|---|---|
| `beta_m2f` | 0.008–0.025 | log | per-act M→F risk for the reference band. Floor from [018](../018_adopt_and_size/SUMMARY.md)'s establishment map (below 0.008 epidemics die and pollute the ensemble — 014's 22% dead draws dragged its envelope to the floor); ceiling well past 0.014, where [014](../014_prior_expansion/SUMMARY.md) obs 4 found the best draws pinned and still rising. |
| `rel_beta_f2m` | 0.15–0.60 | log | the male:female per-act direction ratio. Currently 0.25. EMOD's implied value for its ≥25 band is ~0.35, and Boily's meta-analysis puts the ratio nearer 0.5. Per partnership-year, HIVsim currently gives men **26% less** acquisition risk than EMOD (18.1% vs 24.5%), so 0.25 looks low. |
| `s_f(15–24)` | 0.8–3.0 | log | young women's elevated susceptibility. Biologically supported — cervical ectopy, immature genital mucosa, HSV-2 acquisition around debut. Currently 1.7; EMOD fitted 1.72 (4.894/2.844). Prior extends **below 1.0** so the data can reject the mechanism rather than being assumed into it. |
| `rel_init_prev` | 0.1–0.5 | natural | 1985 seed prevalence multiplier. 014 obs 1 ranked it second only to `beta_m2f` (ρ = 0.59–0.62). Floor at 0.1 per 018's establishment map. |

### Mixing and risk structure — added after the four-parameter set was rejected as too narrow

The four-parameter version was **wrong, and wrong in a specific way**: it argued
at length that the male 25–34 residual is an *exposure* problem — age mixing,
risk-group composition, migration — and then opened only susceptibility
parameters while fixing every exposure one.

It also missed the mechanism entirely. `StructuredSexual.age_diff_pars` is the
direct parameterisation of age-disparate partnerships — mean and SD of partner
age gap, by the woman's age band and risk group — and `make_sim` never touches
it, so it has been at stisim defaults throughout. It is what experiment 011 was
about and what `references/Ott_Age_Gaps_AHRI.pdf` covers.

| parameter | prior | scale | mechanism |
|---|---|---|---|
| `age_gap_shift` | −2 to +3 yr | linear | additive shift on all nine `age_diff_pars` **means**. Sets *which* men women seek |
| `age_gap_sd_mult` | 0.6–1.8 | log | multiplier on all nine **SDs** — the assortativity. A tight SD locks cohorts together, a wide one mixes them |
| `prop_f0` | 0.45–0.85 | linear | female low-risk share. ρ = **0.45** against age-at-peak female prevalence in the shape re-analysis of 014 — the strongest signal after `beta_m2f`, and dropped from the four-parameter set on 014's *level*-only ranking, the exact error `plot_fit_progression.py` exists to prevent repeating |
| `prop_m0` | 0.40–0.80 | linear | male low-risk share. Weaker (ρ ≈ 0.19–0.24) but it is the male exposure side, where the worst residual is |
| `conc_mult` | 0.5–2.0 | log | one multiplier on `f1/f2/m1/m2_conc`. Included on request; expectations low — `m1_conc` measured ρ = 0.01 on level and 0.14 on shape, the weakest thing in 014's set |

Two scalars rather than eighteen numbers for `age_diff_pars`, and one rather
than six for concurrency: the "one calibration parameter modifying many model
parameters through a scaling relationship" pattern, which preserves the relative
structure while giving each hypothesis a single degree of freedom.

### How much of the male age residual these can actually reach

Verified in the matcher: male partnering eligibility is exactly
`over_debut & male & (partners < concurrency)`. No male age taper, no male upper
age cut, no male age term anywhere — where **women** get an explicit
seeking probability that declines to zero at 55. Male risk is structurally
age-flat given debut and concurrency.

But matching is **demand-driven**: women draw a desired partner age
(`own_age + gap`) and take the closest available man at or above it. So the
model does produce an emergent male age-risk profile — a convolution of the
female age distribution with the age-gap distribution. `age_gap_shift` and
`age_gap_sd_mult` are therefore genuine levers on it.

What they **cannot** produce is a *non-smooth* male profile. A convolution is
smooth by construction, so a trough localised to 25–34 while 35–44 is left
alone is unreachable at any parameter value. **That is the pre-registered
failure mode**: if wave 1 closes most of the male residual, the defect was
mixing; if a 25–34-specific trough survives, male age-dependent risk behaviour
is a genuine structural gap and that defines the next model experiment.

### Two warnings from the matcher's own docstring

> *"in a calibration setting, age_diff_pars is not orthogonal to stable_dur_pars
> (concurrency), as concurrency factors influence the available pool of females
> and males looking for partners"*

So expect `age_gap_*` to correlate with `conc_mult` in the orthogonality check.
Forewarned, and an argument for keeping concurrency to one scalar.

> *"a small downward bias of female-male relationship age gaps of < 1 year ...
> driven by the matching of young men"*

The **realized** mean gap runs slightly below the parameter. A known offset, not
a bug, but it means `age_gap_shift` should be read against realized gaps rather
than assumed to equal them.

## What stays fixed, and why — the harder half of the decision

**No male age bands.** This was proposed and then withdrawn, and the reasoning
matters more than the conclusion. **Male per-act susceptibility is essentially
age-flat biologically.** The candidate modifiers all fail: circumcision is
already modelled explicitly (`eff_circ`, the VMMC intervention) so putting it
here would double-count it; ulcerative STI co-infection is real and
age-varying, but it is a co-factor for an unmodelled mechanism rather than
intrinsic susceptibility; foreskin surface area and hygiene have no adult age
gradient.

So `s_m` by age would not be a biological parameter — it would be a **residual
sink**, absorbing misspecification from age mixing, risk-group composition by
age, and male labour migration (substantial in Eswatini, and entirely absent
from the model). Male incidence is driven by *whose* prevalence men are exposed
to, which is exposure, not susceptibility.

That matters here more than in a generic calibration: **age-specific
susceptibility is the quantity that determines who benefits from age-targeted
PrEP.** Fitting a migration pattern with a susceptibility multiplier would make
the calibration look better and the decision worse, and no calibration
diagnostic would reveal it.

**And there is no data to identify them from anyway.** Of the eight incidence
target rows, four have confidence intervals reaching zero — *all three* male
2016 age bands among them. The age-stratified male incidence carries
essentially no information.

**`s_f(35+)` also stays fixed.** Post-menopausal changes plausibly *increase*
susceptibility if anything, so a fitted value below 1.0 would have no story. The
+4.5 pp excess at women 35–44 is far more likely cohort and mortality dynamics,
which [019](../019_age_dependent_survival/SUMMARY.md) and
[022](../022_survival_and_vls/SUMMARY.md) showed the model gets wrong.

| fixed | value | authority |
|---|---|---|
| `mort_mult` | 1.0 | 014 (ρ = −0.01 over a 6× range); 019 obs 1 measured the route it acts on at ~27% of deaths |
| `dur_latent_mult` | 1.0 | **022 reverses 019's recommendation.** Its death gain is generic to any late-life shortening and always costs prevalence, so it would trade targets rather than fit them |
| age gradient on survival | none | 019 (shortening) and 022 (pivoting) both rejected |
| `eff_condom` | 0.85 | literature-supported (~80–90%); condom *use* is a data input and matters more than efficacy |
| `rel_dur_on_art`, `rel_death_f` | as is | 014 ρ ≤ 0.14 on both level and shape; 017 arm D refuted `rel_death_f` |
| `s_m` (all ages), `s_f(35+)` | 1.0 | above |

Leaving the male age pattern and the female 35–44 excess **unfittable is
deliberate.** If wave 1 cannot reproduce them with four well-motivated
parameters, that localises the defect to mixing or behaviour and defines 025.
If instead the calibration is handed six free multipliers it will fit them,
and we will learn nothing about their cause.

## Targets

| target | rows | status | source |
|---|---|---|---|
| PHIA prevalence by age/sex/year | 52 | **hard** | 2007/2011/2016; 2 of 54 dropped |
| **Incidence by age/sex** | **8** | **soft, CI-weighted** | SHIMS1 2011 cohort; SHIMS2 Table 5.3.B |
| F:M incidence ratio 15–49 | 2 | soft, derived | robust to recency-assay MDRI assumptions |
| UNAIDS AIDS deaths | trajectory | **down-weighted** | Spectrum |
| UNAIDS 1985 age distribution | — | hard | — |
| 2021 incidence | — | **HOLD-OUT, untouched** | SHIMS3 Table 5.1 |

**Why incidence is added.** Everything these four parameters govern is a *rate*;
prevalence is a *stock* that confounds incidence with mortality, ART and cohort
history. That is precisely why seven experiments could not separate the
candidate mechanisms — they all produce similar prevalence deviations but
different incidence profiles. It also makes the 2021 hold-out worth more: after
fitting 2011 and 2016 incidence, 2021 becomes a genuine out-of-period
prediction. See [`incidence_construction.py`](../../incidence_construction.py)
for provenance and the three cautions the likelihood must respect.

**Why the two dropped prevalence strata.** 2007 M 15–19 and 2007 F 60–64 hold
fewer than 5 expected infected agents even at N = 50 000
([020](../020_model_sizing/SUMMARY.md) obs 1–2). A coverage fraction computed
over strata the model cannot resolve is not a measure of the model. Stated as a
limitation, not done quietly.

**Why deaths are down-weighted.** The model reaches at most 62% of the UNAIDS
peak anywhere in the transmission prior (018's grid), and 019 and 022 each
bought ~9 points of it at a prevalence cost. As a hard target it would force the
calibration into a trade the data cannot resolve, and likely return an empty
NROY space.

## Plan

**A 400-draw prior sensitivity sweep**, 1 replicate per draw, N = 10 000, at the
four parameters above. Deliberately not the coverage check — 024 runs that at
N = 20 000 with 10 replicates per 020's sizing. This is the smaller,
cheaper question: *does the prior move the targets, and are the parameters
separable?*

Measured:

1. **Spearman ρ of each parameter against each target family** — prevalence by
   age band and sex, incidence by band and sex, the F:M ratio, peak deaths.
   Against *shape* statistics as well as levels, because 014's ranking used a
   single stratum and would have scored a shape parameter at ~0 (the error
   `plot_fit_progression.py` was built to avoid repeating).
2. **Pairwise |ρ| between parameters among draws that fit well** — flagging any
   pair above 0.8 as practically non-identifiable.
3. **Does `s_f(15–24)` move the young-women band specifically**, or does it just
   rescale everything? If the latter it is redundant with `beta_m2f` and should
   be dropped.
4. **Does `rel_beta_f2m` move the F:M incidence ratio** toward the observed
   1.90–2.04, and at what value?
5. **Failed-epidemic fraction.** Must be ~0 given the floors; any failures mean
   the establishment map does not transfer to this parameter set.

## Success criteria

- **Clean:** every parameter has |ρ| > 0.3 against at least one target it is
  meant to control, no pair exceeds |ρ| = 0.8, and no draws produce dead
  epidemics. The prior goes to 024 unchanged.
- **A parameter is redundant:** `s_f(15–24)` turns out to act like `beta_m2f`, or
  `rel_beta_f2m` cannot move the sex ratio. Drop it and record why — a
  three-parameter wave 1 is a better outcome than four with a passenger.
- **The prior cannot reach the data:** visible here as targets lying outside the
  ensemble. That is 024's job to quantify, but if it is obvious at 400 draws the
  bounds get revisited before spending 20 000-agent compute on it.

## Not in scope

- The coverage fraction itself — 024, at N = 20 000 with 10 replicates.
- Likelihood functional form and component weights — `likelihood-design`, once
  these targets are settled.
- Method choice. CLAUDE.md names history matching + trajectory selection as the
  leading method; `method-selection` confirms it after the coverage check.
- Age mixing, risk-group composition and migration. These are the live
  candidates for the male 25–34 residual, and the reason it is left unfittable
  here. **025**, informed by whether wave 1 can fit without them.
