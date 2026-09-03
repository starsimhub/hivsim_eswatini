# Proposals for stisim, from the Eswatini model development

**Against:** stisim 1.5.11 (`97b8bc1`), starsim 3.5.2 (`ceddb877`).
**From:** `starsimhub/hivsim_eswatini`.

Four proposals and two bug reports surfaced by experiments on the Eswatini model
dev.

**These should land as separate PRs**, not one — 1 and 2 are small and
independent, 3 is a design conversation, 4 needs nothing. They are collected here
because they came out of one piece of work.

Every workaround referenced lives as a subclass in our repo, never as a patch to
the stisim checkout — we lost a VMMC fix that way once, and the checkout is
currently a pinned dependency of a running calibration.

Ordered by importance:

---

## 1. `vls_coverage` should be a stock target

`on_effective_art` is written only in `start_art()`, so suppression is fixed at
ART initiation for life and the TLD transition is unrepresentable. stisim has
already applied stock-target semantics to VMMC (1.5.9) and to PrEP — ART
suppression is the last coverage-like quantity still applied as a flow.

Your own comments, in `interventions/hiv_interventions.py`:

```python
:690   # Coverage is a circumcision *prevalence* (stock) target, matching ...
:878   # Coverage is a PrEP *prevalence* (stock) target, matching how PrEP ...
```

The VMMC change is why we were able to delete our `VMMCPrevalenceTarget`
subclass. This is the same fix, one intervention over.

**Why it matters beyond tidiness.** Eswatini's rise in the third 95 came
substantially from better regimens and adherence support among people *already
on treatment*. By 2021 the treated stock is most of the treated population, so
flow-only semantics misses nearly all of it. Measured against the input series
it was given:

| | 2016 m | 2016 f | 2021 m | 2021 f |
|---|---|---|---|---|
| input | 0.905 | 0.918 | 0.967 | 0.959 |
| flow only (current) | 0.901 | 0.917 | **0.949** | **0.949** |
| stock target | 0.911 | 0.922 | 0.967 | 0.959 |

Flow-only lags its own input by 1.8 points for men in 2021 — exactly when
suppression is improving fastest.

**Implementation.** Ours is `vls_stock_target.VLSStockTarget`: each step, rank
agents on ART by a persistent per-agent `suppression_propensity` and mark the top
`target × n` suppressed. It deliberately mirrors `sti.VMMC`'s `willingness`.
Ranking by a persistent score rather than redrawing matters twice over — agents
don't churn between states, and a rising target *adds* to the suppressed pool
rather than reshuffling it. Because `update_transmission()` recomputes
`rel_trans` from `on_effective_art` each step and `get_art_mortality_hazard()`
reads the same state, flipping those booleans is sufficient. Natural home is
inside `sti.ART` itself, so `vls_coverage` simply *means* a stock target,
consistent with `coverage`.

---

## 2. No result reports the effective / non-suppressive ART split

Population viral suppression — the third 95 — is not observable from standard
outputs, despite a huge difference in residual transmission between the two
states (`effective_art_efficacy = 0.99` against `nonsupp_art_efficacy = 0.35`).
`hiv.p_on_art` pools them. Anyone comparing a model to a PHIA cascade has to
write their own analyzer.

Small addition to the existing `analyzers.art_coverage`, which already does
age × sex banded results and is the natural home:

```python
results.append(ss.Result(f'n_vls_{sex}_{lo}_{hi}', dtype=int))
results.append(ss.Result(f'p_vls_{sex}_{lo}_{hi}', scale=False))       # of PLHIV
results.append(ss.Result(f'p_vls_given_art_{sex}_{lo}_{hi}', scale=False))
```

The distinction between the last two is worth documenting there: `p_vls` (of all
PLHIV) is the GAM 1.3 indicator and a cascade *outcome*; `p_vls_given_art` is
what `vls_coverage` sets as an *input*. Conflating them double-counts the
coverage ramp — we did exactly that for half a day.

Same for per-band `new_infections`, without which age-stratified incidence
cannot be computed — only the 15–49 and 18–49 aggregates. Same file, same fix
shape. We needed it to fit incidence targets published in 5- and 10-year bands.

---

## 3. Untreated survival has no age dependence

`dur_latent = ss.lognorm_ex(ss.years(10), ss.years(3))` applies identically to a
17-year-old and a 55-year-old, against ALPHA-network estimates (~12.5 y for
infection at 15–24 falling to ~7.5 y at 45+) and EMOD's explicit
parameterization — Akullian et al. 2020 (Lancet HIV), appendix Table A2:

```
lambda = 21.182 - 0.2717 * age_at_infection,  Weibull shape 2
  -> mean survival 13.96 y at 20, 11.55 at 30, 9.14 at 40, 6.73 at 50
```

To be straight about the evidence: we implemented this twice and it did not
improve our fit — shortening gradients bought death-peak height at a monotone
cost in prevalence, and the EMOD profile above made our age-stratified fit
worse. So this is offered as *stisim cannot express a documented feature of HIV
natural history*, not as *this improves fit*. If the flat assumption is
defensible for stisim's intended uses, that is a reasonable answer and this can
be closed.

Ours is `hiv_survival.AgeDependentSurvival`, which rescales the latent interval
by age at infection while preserving the drawn `dur_falling`, and consumes no
random numbers — so a multiplier of 1.0 is bit-identical to upstream at the same
seed.

**It comes with useful instrumentation, which may be worth having regardless of
the survival question.** `hiv.new_deaths` pools two death routes, and
`step_die()` clears `ti_zero` and `ti_infected`, so the route is unrecoverable
afterwards. Reading it before calling `super()` splits the counts: **70–74% of
HIV deaths flow through the `ti_zero` route**, which no multiplier touches, and
only 26–30% through `p_hiv_death`, which is what `rel_death`, `rel_death_f` and
the on-ART multipliers act on. Which is why `rel_death` and `rel_death_f` appear
inert.

---

## 4. WITHDRAWN — CD4 death rates

The CD4 death rates were hard-coded literals in 1.5.8 but are
`pars.cd4_death_bins` / `pars.cd4_death_rates` in 1.5.11
(`diseases/hiv.py:33-34`, read at line 401). Already fixed upstream, so our
`hiv_mortality.py` subclass is obsolete and should be retired.

Recorded only because it appears in our experiment 014 and someone will
otherwise wonder why.

---

## 5. Design comment: several defaults are not neutral

| default | what it silently does |
|---|---|
| `sti.Prep()` with `coverage=None` | reaches 80% of FSW starting in 2004 (`{'year': [2004, 2005, 2015, 2025], 'value': [0, 0.01, 0.5, 0.8]}`) |
| `sti.ART(vls_coverage=None)` | 100% of ART initiators suppressed |
| a stratified `vls_coverage` omitting a stratum | that stratum defaults to 100% suppression |

**Suggestion: `None` should mean off, not invent a scale-up scenario.** Or
raise. Where a non-trivial default is genuinely wanted, a warning naming the
assumed values would be enough — the failure mode is silence, not the values.
Defaulting unlisted strata to 100% likewise turns a partial specification into a
confidently wrong one, where raising would catch it immediately.

Both of the first two cost us real work. Every experiment in our project from
001 to 017 ran with a PrEP programme nobody chose, beginning a decade before
PrEP had efficacy evidence; we found it only while porting to 1.5.11. And the
suppression default is worse in kind, because 100% is not a conservative choice
— it biases a *counterfactual comparison* rather than a fit. A scenario that
raises suppression, measured against a baseline that already overstates it by
8.8 percentage points, understates what treatment scale-up buys. No calibration
diagnostic would surface that.

---

## Bug reports (no code attached)

**`on_prep` appears dead in 1.5.8 and 1.5.10.** We did not chase it further
because we removed PrEP from our model entirely, so we have no reproducer beyond
the original observation — happy to build one if useful.

**`rel_death_f` has no documented backing** for its ~26% female survival
advantage at equal CD4. New in 1.5.11. The seroconverter-cohort literature we
are aware of does not clearly support an advantage of that size conditional on
CD4, and we found the model insensitive to it anyway (|z| ≤ 0.9) — because it
acts on the minority death route, see 3. Not arguing it is wrong; asking where
the number came from, and suggesting the source belongs in a comment beside it.

---

## Appendix — not a stisim issue: `historymatching` imports gpflow eagerly

Different repo (`InstituteforDiseaseModeling/history_matching`), kept here so it
is not lost.

`emulators/__init__.py` imports every emulator eagerly, including
`from .gpr import GPR`, and `emulators/gpr.py:7` is a bare `import gpflow`. But
gpflow is referenced only inside GPR's *method bodies* (lines 81–99 and 194),
never at class-definition time. So the package cannot be imported at all without
TensorFlow, even to use `bayes_linear`, which is pure NumPy/SciPy and the
documented first choice. That blocks any Python version TensorFlow has not caught
up to — we are on 3.14 locally, where TensorFlow publishes no wheels, and had to
shim `gpflow` into `sys.modules`. Making the import lazy would fix it in one
line.

Two doc notes for the same package: the module is `historymatching`, not
`history_matching`, and 2.0.1 replaced `HistoryMatchingBuilder` with a single
`HistoryMatching(...)` constructor.
