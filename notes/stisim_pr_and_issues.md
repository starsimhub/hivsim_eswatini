# stisim upstream: the PR, and two issues to file

**Branch pushed:** `feat/vls-stock-target-and-cascade-results` on
`starsimhub/stisim`, commit `45a8b21`, based on `main` at `5d35252` (v1.6.1).

**Open the PR here:**
https://github.com/starsimhub/stisim/pull/new/feat/vls-stock-target-and-cascade-results

Paste the block below as the PR description. The two issues after it are not
code and should be filed separately.

Design context lives in [`stisim_upstream_proposals.md`](stisim_upstream_proposals.md).
Note that document was written against 1.5.11; both proposals in this PR were
re-verified against 1.6.1 before any code was written, and proposal 4 (CD4 death
rates) was already fixed upstream so it is not here.

---

## PR description — paste from here

### Summary

Two related gaps found while calibrating an HIV model for Eswatini against PHIA
and SHIMS data.

1. **`vls_coverage` is applied only at ART initiation**, so viral suppression is
   fixed for an agent's lifetime and a time-varying target has no effect on the
   existing treated stock. Made a stock (prevalence) target, matching what
   `sti.VMMC` and `sti.Prep` already do.
2. **Nothing reports the `on_effective_art` / `on_nonsuppressive_art` split**, so
   population viral suppression — the third 95 — is not observable from standard
   outputs. Added to the `art_coverage` analyzer, along with per-band
   `new_infections` so age-stratified incidence is computable.

### 1. `vls_coverage` as a stock target

`HIV.on_effective_art` is written only in `start_art()` (and cleared on death or
ART discontinuation). An agent who starts treatment in 2010 keeps their 2010
suppression status for life, so a time-varying `vls_coverage` only ever affects
new initiators.

That makes the dominant real-world mechanism unrepresentable. Suppression rises
mostly through **better regimens and adherence support among people already on
treatment** — the dolutegravir/TLD transition from around 2019 being the obvious
case — and by then the existing treated stock is most of the treated population.

**This is the same fix you have already applied twice.** From this file:

```python
:690   # Coverage is a circumcision *prevalence* (stock) target, matching ...
:878   # Coverage is a PrEP *prevalence* (stock) target, matching how PrEP ...
```

VMMC moved to stock semantics in 1.5.9 and PrEP followed. ART suppression was
the last coverage-like quantity still applied as a flow. (The VMMC change is
also why we were able to delete our own `VMMCPrevalenceTarget` subclass — this
PR is the equivalent one intervention over.)

**How.** `ART.vls_stock_correction()` runs each step after coverage correction,
ranks agents on ART by a persistent per-agent `suppression_propensity`, and
corrects the suppressed set to the target. Ranking by a persistent score rather
than redrawing matters twice over: agents don't churn between states every step,
and a rising target *adds* to the suppressed pool rather than reshuffling it. It
also reads as a plausible individual trait. Deliberately mirrors VMMC/PrEP
`willingness`.

Because `update_transmission()` recomputes `rel_trans` from `on_effective_art`
each step and `get_art_mortality_hazard()` reads the same state, flipping those
booleans is sufficient — transmission and on-ART mortality both follow.

**One intentional divergence from VMMC:** this moves agents in *both* directions.
Viral suppression is genuinely reversible (treatment failure, interrupted
adherence) where circumcision is not, so `_suppress_to_target` is not clamped
the way `_circumcise_to_target` is.

Stratified `vls_coverage` is corrected **per stratum** — a single global ranking
would wash out the age/sex differentials in the input data, the same reason
`art_coverage_correction` and the VMMC correction already work per stratum.

**Behaviour change, stated plainly.** `vls_coverage` now means a stock target, so
anyone relying on flow semantics will get different results. That is the point,
and it makes `vls_coverage` consistent with `coverage` beside it.
`vls_coverage=None` still means 100% suppressed, so the default is untouched. If
you'd prefer a deprecation path, a `vls_coverage_mode='stock'|'flow'` flag is a
small addition — say the word and I'll add it.

**Size of the effect.** Measured in our project at N=10,000 over 4 arms × 10
seeds, suppression among the treated against the input series it was given:

| | 2016 m | 2016 f | 2021 m | 2021 f |
|---|---|---|---|---|
| input | 0.905 | 0.918 | 0.967 | 0.959 |
| flow only (current) | 0.901 | 0.917 | **0.949** | **0.949** |
| stock target (this PR) | 0.911 | 0.922 | 0.967 | 0.959 |

Flow-only lags its own input by 1.8 points for men in 2021 — precisely when
suppression is improving fastest and the treated stock is largest. It also
*overshoots* early, because unsuppressed patients have higher on-ART mortality
and the surviving stock is enriched for the suppressed.

### 2. Cascade and incidence results

`effective_art_efficacy = 0.99` against `nonsupp_art_efficacy = 0.35` is roughly
a 65× difference in residual transmission, but neither state is reported and
`hiv.p_on_art` pools them. So the third 95 — the headline PHIA/GAM indicator —
can't be read off standard outputs.

Added to `art_coverage` (aggregate, per sex, per age×sex band):

- `n_vls`, `p_vls` — of all PLHIV. This is the GAM 1.3 indicator and a cascade
  **outcome**.
- `p_vls_given_art` — of those on ART. This is what `vls_coverage` sets as an
  **input**.

The distinction is documented inline because it's easy to get wrong: conflating
them double-counts the coverage ramp. We did exactly that for half a day.

Also added per band: `new_infections`, `n_susceptible`, `incidence`. Without a
per-band numerator there's no way to compute age-stratified incidence from
standard outputs — only the 15–49/18–49 aggregates — yet survey incidence is
published in 5- and 10-year bands. `n_susceptible` is included because it's the
correct denominator; using the whole band understates incidence by a factor of
(1 − prevalence), a ~30% error at Eswatini's prevalence.

All additive. No existing result changes.

### Testing

Verified a rising `vls_coverage` of 0.50 → 0.70 → 0.95 is tracked (0.699 at
2010, 0.950 at 2029), a flat 0.60 holds (0.599, 0.601), and `vls_coverage=None`
still gives 100%.

**I could not get the existing suite green locally, and it is not this PR's
fault.** `tests/test_hiv.py` gives 22 failed / 4 passed and
`tests/test_hiv_interventions.py` 17 failed / 1 passed, all
`TypeError: Module.define_states() got an unexpected keyword argument`. The
counts are **identical on clean `main`**, so these changes introduce no new
failures. The cause is a starsim version mismatch — stisim 1.6.1 requires
`starsim>=3.6.0` and I'm pinned to 3.5.2 for an in-flight calibration. **Please
let CI confirm on a supported starsim.**

## PR description — paste to here

---

# Issue 1 — Several defaults are not neutral and silently invent behaviour

Design observation rather than a bug, raised as a class rather than as three
separate reports because the pattern is what matters.

| default | what it silently does |
|---|---|
| `sti.Prep()` with `coverage=None` | applies `{'year': [2004, 2005, 2015, 2025], 'value': [0, 0.01, 0.5, 0.8]}` to FSW — a programme reaching 80% of FSW, **starting in 2004** |
| `sti.ART(vls_coverage=None)` | **100%** of those on ART are virally suppressed |
| a stratified `vls_coverage`/`coverage` omitting a stratum | that stratum defaults to **100%** |

**The PrEP one cost us real work.** Every experiment in our project from 001 to
017 ran with a fabricated PrEP programme nobody chose, beginning a decade before
PrEP had efficacy evidence. We only found it while porting to 1.5.11. We measured
it at about 4% of 2021 prevalence and ~0.67 realised protection of uninfected FSW
by 2021 — not enormous, but silently absorbing transmission parameters we were
trying to calibrate, which is the kind of bias that doesn't announce itself.

**The suppression default is worse in kind**, because 100% suppression is not a
conservative choice — it makes every treated agent effectively non-infectious at
`effective_art_efficacy = 0.99`. We measured it as overstating population viral
suppression by **8.8 percentage points** in 2016. Critically that biases a
*counterfactual comparison* rather than a fit: a scenario that raises suppression,
measured against a baseline already overstating it, understates what treatment
scale-up buys. No calibration diagnostic would surface that.

**Suggestion.** For any parameter representing programme coverage, `None` should
mean **off**, or raise. Not "invent a plausible-looking programme". Where a
non-trivial default is genuinely wanted, a warning naming the assumed values
would be enough — the failure mode is silence, not the values themselves. The
same applies to unlisted strata: defaulting them to 100% turns a partial
specification into a confidently wrong one, where raising would catch it
immediately.

---

# Issue 2 — Untreated survival has no age dependence, and `rel_death_f` has no cited source

Two things, related through the death routes.

## Untreated survival is age-independent

`dur_latent = ss.lognorm_ex(ss.years(10), ss.years(3))` applies identically to a
17-year-old and a 55-year-old, so total untreated survival doesn't depend on age
at infection. ALPHA-network pooled estimates put untreated survival at roughly
12.5 y for infection at 15–24 falling to ~7.5 y at 45+. EMOD parameterises it
directly — Akullian et al. 2020 (Lancet HIV), appendix Table A2:

```
lambda = 21.182 - 0.2717 * age_at_infection,  Weibull shape 2
  -> mean survival 13.96 y at 20, 11.55 at 30, 9.14 at 40, 6.73 at 50
```

**Being straight about our evidence: we implemented this twice and it did not
improve our fit.** Shortening gradients bought death-peak height at a monotone
cost in prevalence; the EMOD profile above made our age-stratified fit
substantially worse (MAE 0.0584 → 0.0817), because its realised effect is
asymmetric (+0.67 y at 15–24 against −4.03 y at 45+).

So this is filed as *stisim cannot express a documented feature of HIV natural
history*, **not** as *this improves fit*. **No PR attached deliberately** — we
don't use the code, and shipping code we rejected would be bad faith. If the flat
assumption is defensible for stisim's intended uses, that's a reasonable answer
and this should be closed.

## Most HIV deaths bypass every mortality multiplier

Found while instrumenting the above, and useful independently. `hiv.new_deaths`
pools two death routes, and `step_die()` clears `ti_zero` and `ti_infected`, so
the route is unrecoverable afterwards. Reading it before calling `super()` splits
the counts: over 1985–2026 at N=10,000, **70–74% of HIV deaths flow through the
`ti_zero` route**, which no multiplier touches, and only 26–30% through
`p_hiv_death`, which is what `rel_death`, `rel_death_f` and the on-ART
multipliers act on.

That explains why those knobs look inert in calibration, and it isn't currently
visible to users. A result splitting deaths by route would be a small, useful
addition.

## `rel_death_f` provenance

New in 1.5.11, encoding "at the same CD4, women die roughly 26% slower". The
seroconverter-cohort literature we're aware of doesn't clearly support a female
survival advantage of that size *conditional on CD4*, and we found the model
insensitive to it anyway (|z| ≤ 0.9) — because it acts on the minority death
route above. Not arguing it's wrong; asking where the number came from, and
suggesting the source belongs in a comment beside it.

## Also, not a stisim issue — `historymatching` imports gpflow eagerly

For `InstituteforDiseaseModeling/history_matching`, recorded here so it isn't
lost. `emulators/__init__.py` imports every emulator eagerly including
`from .gpr import GPR`, and `emulators/gpr.py:7` is a bare `import gpflow` — but
gpflow is referenced only inside GPR's *method bodies* (lines 81–99 and 194),
never at class-definition time. So the package can't be imported at all without
TensorFlow, even to use `bayes_linear`, which is pure NumPy/SciPy and the
documented first choice. That blocks any Python version TensorFlow hasn't caught
up to — we're on 3.14, where TensorFlow publishes no wheels, and had to shim
`gpflow` into `sys.modules`. Making the import lazy would fix it in one line.
