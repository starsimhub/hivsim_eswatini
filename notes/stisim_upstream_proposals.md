# Proposals for stisim, from the Eswatini HIV calibration

**From:** Adam Akullian (Gates Foundation), `starsimhub/hivsim_eswatini`
**Date:** 2026-09-03
**Against:** stisim 1.5.11 (`97b8bc1`), starsim 3.5.2 (`ceddb877`)

Four proposals and two bug reports, all surfaced by calibrating the Eswatini HIV
model over experiments 016–024. Each names the gap, the evidence, and what we
did locally to work around it, so you can judge whether the workaround belongs
upstream.

**These should be separate PRs**, not one. They are collected here because they
came out of one piece of work; #1 and #2 are small and independent, #3 is a
design conversation, #4 is a one-line observability addition. Ordered by how
strong we think the case is.

Every workaround referenced lives as a subclass in our repo, never as a patch to
the stisim checkout — we lost a VMMC fix that way once and don't intend to
repeat it.

---

## 1. `vls_coverage` should be a stock target, like VMMC and PrEP already are

**Strongest of these, and it closes an inconsistency in your own pattern.**

`sti.ART(vls_coverage=...)` accepts a time-varying series, but it is applied
**only at ART initiation**. `HIV.on_effective_art` is written in exactly three
places: `start_art()`, and cleared on death (`step_die`) and on ART
discontinuation. So an agent who starts treatment in 2010 keeps their 2010
suppression status for the rest of their life.

That makes the dominant real-world mechanism unrepresentable. Eswatini's rise in
the third 95 came substantially from **better regimens and adherence support** —
the dolutegravir (TLD) transition from around 2019 — which improve suppression
among people *already on treatment*. By 2021 the existing treated stock is most
of the treated population, so flow-only semantics misses nearly all of it.

**You have already accepted this argument twice.** In 1.5.11:

```python
# interventions/hiv_interventions.py:690
# Coverage is a circumcision *prevalence* (stock) target, matching ...

# interventions/hiv_interventions.py:878
# Coverage is a PrEP *prevalence* (stock) target, matching how PrEP ...
```

VMMC moved to stock-target semantics in 1.5.9 — which is why we were able to
delete our `VMMCPrevalenceTarget` subclass in our experiment 018. PrEP followed.
**ART suppression is the one remaining coverage-like quantity still applied as a
flow.**

### Evidence of the size of the effect

Measured in our experiment 022 (4 arms × 10 seeds, N = 10 000), suppression
among the treated, against the input series it was given:

| | 2016 m | 2016 f | 2021 m | 2021 f |
|---|---|---|---|---|
| input | 0.905 | 0.918 | 0.967 | 0.959 |
| flow only (current) | 0.901 | 0.917 | **0.949** | **0.949** |
| stock target | 0.911 | 0.922 | **0.967** | **0.959** |

Flow-only lags its own input by 1.8 points for men in 2021 — precisely when
suppression is improving fastest and the treated stock is largest. It also
*overshoots* early (0.935 against an input of 0.905 in a 3000-agent run),
because unsuppressed patients have higher on-ART mortality and the surviving
stock is enriched for the suppressed.

### Suggested implementation

Ours is `vls_stock_target.VLSStockTarget`, an `ss.Intervention` that runs after
ART and, each step, ranks agents currently on ART by a **persistent per-agent
`suppression_propensity`** and marks the top `target × n` as suppressed. It
deliberately mirrors `sti.VMMC`'s `willingness` pattern.

Ranking by a persistent score rather than redrawing matters twice over: agents
don't churn between suppressed and unsuppressed, and a rising target *adds* to
the suppressed pool rather than reshuffling it. It also reads as a plausible
individual trait — adherence propensity is persistent, not redrawn monthly.

Because `update_transmission()` recomputes `rel_trans` from `on_effective_art`
each step and `get_art_mortality_hazard()` reads the same state, flipping those
booleans is sufficient — transmission and on-ART mortality both follow.

Natural home upstream is inside `sti.ART` itself rather than a separate
intervention, so `vls_coverage` simply *means* a stock target, consistent with
`coverage`.

---

## 2. No result reports the effective / non-suppressive ART split

`HIV` maintains `on_effective_art` and `on_nonsuppressive_art`, and the
difference between them is large — `effective_art_efficacy = 0.99` against
`nonsupp_art_efficacy = 0.35`, a 65× difference in residual transmission. But
**neither state is reported in any result.** `hiv.p_on_art` pools them.

The consequence is that **population viral suppression — the headline PHIA
indicator and the third 95 — is not observable from standard outputs.** Anyone
comparing a model to a PHIA cascade has to write their own analyzer, as we did.

This looks like a small addition to the existing `analyzers.art_coverage`, which
already does age × sex banded results with configurable bins and is the natural
place for it:

```python
results.append(ss.Result(f'n_vls_{sex}_{lo}_{hi}', dtype=int))
results.append(ss.Result(f'p_vls_{sex}_{lo}_{hi}', scale=False))       # of PLHIV
results.append(ss.Result(f'p_vls_given_art_{sex}_{lo}_{hi}', scale=False))
```

The distinction between the last two matters and is worth documenting there:
`p_vls` (of all PLHIV) is the GAM 1.3 indicator and a cascade *outcome*, while
`p_vls_given_art` is the quantity `vls_coverage` sets as an *input*. Conflating
them double-counts the coverage ramp — we did exactly that for half a day.

**Related, smaller:** `art_coverage` gives banded counts but there is no
per-band `new_infections`, so **age-stratified incidence cannot be computed**
from standard outputs either — only the 15–49 and 18–49 aggregates. We added it
to our own analyzer for experiment 023 to fit the SHIMS2 incidence targets,
which are published in 5- and 10-year bands. Same fix shape, same file.

---

## 3. Untreated survival has no age dependence

`dur_latent = ss.lognorm_ex(ss.years(10), ss.years(3))` is applied identically
to a 17-year-old and a 55-year-old, so total untreated survival has no
dependence on age at infection.

The external evidence is reasonably strong. ALPHA-network pooled estimates put
untreated survival at roughly 12.5 y for infection at 15–24 falling to ~7.5 y at
45+. EMOD parameterises it directly — Akullian et al. 2020 (Lancet HIV),
appendix Table A2:

```
HIV_Adult_Survival_Scale_Parameter_Intercept   21.182
HIV_Adult_Survival_Scale_Parameter_Slope       -0.2717
HIV_Adult_Survival_Shape_Parameter             2
    lambda = 21.182 - 0.2717 * age_at_infection,  Weibull shape 2
    mean survival = lambda * Gamma(1.5)
      -> 13.96 y at 20, 11.55 at 30, 9.14 at 40, 6.73 at 50
```

**We want to be straight about our evidence here: we tested this twice and it
did not help our fit.** Experiment 019 tried shortening gradients (multipliers
≤ 1.0) and found they bought ~9 points of the UNAIDS death peak at a monotone
cost in prevalence. Experiment 022 tried the EMOD profile above, which *pivots*
— longer at young ages, much shorter at old — and it made our age-stratified fit
substantially worse (MAE 0.0584 → 0.0817; strata inside the PHIA CI 23 → 9),
because the realised effect is asymmetric: +0.67 y at 15–24 against −4.03 y at
45+.

So this is offered as **"stisim cannot express a documented feature of HIV
natural history"**, not as "this improves fit". If you think the flat assumption
is defensible for stisim's intended uses, that is a reasonable answer and this
proposal should be closed.

If it is wanted, ours is `hiv_survival.AgeDependentSurvival`, which overrides
`set_prognoses()` to rescale the latent interval by age at infection while
preserving the drawn `dur_falling` (so the CD4 decline rate, and hence
`rel_trans_falling`, is untouched). It consumes no random numbers of its own, so
a multiplier of 1.0 is bit-identical to upstream at the same seed — verified on
deaths, prevalence and infections.

**A useful piece of instrumentation came with it, which may be worth having
regardless of the survival question.** `hiv.new_deaths` pools two death routes,
and `step_die()` clears `ti_zero` and `ti_infected`, so the route is
unrecoverable afterwards. Our subclass reads it before calling `super()` and
splits the counts. Measured over 1985–2026 at N = 10 000, **70–74% of HIV deaths
flow through `ti_zero`** — the route with no multiplier of any kind — and only
26–30% through `p_hiv_death`, which is what `rel_death`, `rel_death_f` and the
on-ART multipliers act on. That explains why those knobs appear inert in
calibration, and it is not currently visible to users.

---

## 4. `hiv_mortality.py` — withdrawn, already fixed upstream

Recorded for completeness because it appears in our experiment 014. In 1.5.8,
`make_p_hiv_death()` embedded the CD4-stratified annual death rates as literals
in the method body, so they could not be calibrated. In 1.5.11 they are
`pars.cd4_death_bins` and `pars.cd4_death_rates` (`diseases/hiv.py:33–34`), read
from `pars` at line 401.

**Fixed upstream; no action needed.** Our `hiv_mortality.HIVMortalityMultiplier`
subclass is now obsolete and we are retiring it.

---

## Comment: several defaults are not neutral, and silently invent behaviour

This is a design observation rather than a bug, and it cost us real work, so it
seems worth raising as a class rather than as individual instances.

Three defaults in stisim do something substantive when the user supplies
nothing:

| default | what it silently does |
|---|---|
| `sti.Prep()` with `coverage=None` | falls back to `{'year': [2004, 2005, 2015, 2025], 'value': [0, 0.01, 0.5, 0.8]}` on FSW — a PrEP programme reaching 80% of FSW, **starting in 2004** |
| `sti.ART(vls_coverage=None)` | **100%** of ART initiators are virally suppressed |
| a stratified `vls_coverage` that omits a stratum | that stratum defaults to **100%** suppression |

The PrEP one is the costly example. **Every experiment in our project from 001
to 017 ran with a fabricated PrEP programme nobody chose**, beginning a decade
before PrEP had efficacy evidence. We only found it while porting to 1.5.11. Our
experiment 018 measured it at about 4% of 2021 prevalence and ~0.67 realised
protection of uninfected FSW by 2021 — not enormous, but it was silently
absorbing transmission parameters we were trying to calibrate, which is exactly
the kind of bias that does not announce itself.

The suppression default is worse in kind, because 100% suppression is not a
conservative choice: it makes every treated agent effectively non-infectious at
`effective_art_efficacy = 0.99`. Our experiment 021 measured it as overstating
population viral suppression by **8.8 percentage points** in 2016. And critically
that is a bias in a *counterfactual comparison* rather than in the fit — a
scenario that raises suppression, measured against a baseline that already
overstates it, understates what treatment scale-up buys. No calibration
diagnostic would surface that.

**Suggestion.** For any parameter representing programme coverage, `None` should
mean **off**, or raise. Not "invent a plausible-looking programme". Where a
non-trivial default is genuinely wanted, emitting a warning naming the values
being assumed would be enough — the failure mode is silence, not the values
themselves. The same applies to unlisted strata in a stratified table:
defaulting them to 100% turns a partial specification into a confidently wrong
one, where raising would catch it immediately.

---

## Bug reports (no code attached)

**`on_prep` appears dead in 1.5.8 and 1.5.10.** Recorded in our experiment 017,
observation 4. We did not chase it further because 018 removed PrEP from our
model entirely, so we have no reproducer beyond the original observation — happy
to build one if useful.

**`rel_death_f` has no documented provenance.** New in 1.5.11, it encodes "at
the same CD4, women die roughly 26% slower". The seroconverter-cohort literature
we are aware of does not clearly support a female survival advantage of that
size conditional on CD4, and our experiment 017 arm D found the model insensitive
to it anyway (|z| ≤ 0.9) — because it acts on the minority death route (see #3).
Not arguing it is wrong; asking where the number came from, and suggesting the
source belongs in a comment beside it.

---

## Not a stisim issue: `historymatching` imports gpflow eagerly

Different repo (`InstituteforDiseaseModeling/history_matching`), noted here so
it is not lost.

`historymatching/emulators/__init__.py` imports every emulator eagerly,
including `from .gpr import GPR`, and `emulators/gpr.py:7` is a bare
`import gpflow`. But gpflow is referenced **only inside GPR's method bodies**
(`gpflow.models.GPR`, `gpflow.kernels.SquaredExponential`, `gpflow.Parameter`,
`gpflow.optimizers.Scipy`, `gpflow.utilities` — lines 81–99 and 194), never at
class-definition time.

Consequently the package cannot be imported at all without TensorFlow, even to
use `bayes_linear`, which is pure NumPy/SciPy and the documented first choice.
That matters on any Python version TensorFlow has not caught up to — we are on
3.14 locally, where TensorFlow publishes no wheels, and had to write a shim that
stubs `gpflow` in `sys.modules` before importing.

Making the import lazy (inside the methods, or a `try/except ImportError` around
the `from .gpr import GPR` line) would fix it in one line.

Two smaller documentation notes for the same package: the module is
`historymatching`, not `history_matching`, and 2.0.1 appears to have replaced
`HistoryMatchingBuilder` with a single `HistoryMatching(...)` constructor — some
docs still describe the builder chain.
