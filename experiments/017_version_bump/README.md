# Exp 017 — Version bump: stisim 1.5.8 → 1.5.11, starsim 3.5.0 → 3.5.2

*The population-size and replicate sweep that [016](../016_double_counted_mortality/SUMMARY.md)
queued for this slot moves to 018. stisim 1.5.11 released 2026-08-24 and changes
HIV mortality by default — measuring model setup on 1.5.8 would mean measuring
it twice.*

## Question

We run stisim 1.5.8 / starsim 3.5.0. `pyproject.toml` declares
`stisim>=1.5.10, starsim>=3.5.2`, so declared and actual have disagreed since
PR #2. Three releases have landed in between, and 1.5.11's headline breaking
change is this:

> **On-ART mortality is now nonzero by default.** Previously, agents on ART had
> zero HIV-specific mortality; on-ART death now follows an age/sex/CD4/adherence-adjusted
> rate anchored to (and always ≤) the off-ART CD4-based hazard. Models calibrated
> against 1.5.10 or earlier will need recalibration.

016's central finding was that **the HIV module under-supplies AIDS deaths** —
~65 % of the UNAIDS peak, ~50 % post-2010 — and its first named candidate cause
was *"nobody on ART can die of HIV (`ti_zero` is cleared at ART initiation,
`hiv.py:616-617`)"*. That is precisely what 1.5.11 changes. So the experiment
016 queued as **018 — the HIV mortality deficit** may be partly resolved by an
upgrade rather than by modelling work.

**This experiment asks: how much does the bump change the model, and how much of
the AIDS-death deficit does 1.5.11's on-ART mortality close?**

The bump has to happen regardless — the question is whether we learn something
from it or just absorb it silently.

## Plan

**Design: three arms at fixed parameters**, the same shape as 015 and 016.
Same parameters, same seeds throughout; **each adjacent pair differs by exactly
one thing.**

| Arm | stisim / starsim | background mortality |
|---|---|---|
| **A** | 1.5.8 / 3.5.0 | `data/eswatini_deaths.csv` as shipped (all-cause) |
| **B** | 1.5.11 / 3.5.2 | as shipped (all-cause) |
| **C** | 1.5.11 / 3.5.2 | HIV-deleted (016's construction) |

- **A → B isolates the version.** Arm A should reproduce 016's arm A exactly —
  that is the control, and the first thing to check.
- **B → C isolates the mortality data**, reproducing 016's A/B on the new stack.

**Why C exists.** 1.5.11 adds on-ART HIV mortality on top of background rates
that *already contain* AIDS deaths, so arm B double-counts harder than arm A
does — and the two death machines compete for the same agents. 016's
observation 5 measured this: deleting background AIDS *raised* HIV-module
deaths from 3 376 to 3 856 at 2005 (+14 %), because PLHIV who would have been
killed by inflated background mortality survived to die of HIV instead.
Measuring 1.5.11's on-ART mortality against double-counted rates therefore
**understates** how much of the deficit it closes. Arm C is the unbiased
measurement, and is the configuration 018 starts from.

**Parameters.** 016's two sets — defaults and high-transmission — 10 seeds
each. 3 arms × 2 sets × 10 seeds = 60 sims, ~13 min locally at 016's rates.
Using 016's exact configuration is deliberate: it makes 017 arm A a
reproduction check, and arms B and C directly comparable to 016's two published
arms.

**The HIV-deleted data is regenerated, not promoted.** Arm C rebuilds it by
importing 016's `mortality_construction.py` and writing to
`outputs/data_hiv_deleted/`, exactly as 016 did. Promoting those rates to a
repo-level input under `data/` remains **018's** decision — 016 asked for that
to be deliberate, and a three-arm comparison does not require it.

**The arms cannot all run simultaneously** — both stisim and starsim are
editable installs from a single checkout each, so the version is a property of
the environment, not of the run. Arm A runs first on the current stack and its
results are written to `outputs/`; then the checkouts move to `v1.5.11` /
`v3.5.2` and arms B and C run together. Recorded here because it means arm A is
not re-runnable after the bump without moving the checkouts back — so `run.py`
resumes from per-sim parquet in `outputs/sims/` and refuses to re-run a cell
whose file already exists.

### Porting the breakage

Known call sites, from a scan of the repo against 1.5.11's breaking changes:

| Site | Break | Expected fix |
|---|---|---|
| [vmmc.py](../../vmmc.py) | `self.circumcised`, `self.ti_circumcised`, `self.pars.eff_circ` moved to `HIV` | Upstream `sti.VMMC` is now itself a prevalence/stock-target intervention with per-stratum top-up — the semantics this subclass existed to provide (it gained them in 1.5.9). `run.py` passes `vmmc_class` per version, so **no repo edit is needed to run**. Deleting the file is an adoption step for 018, once metric 5 confirms equivalence |
| [interventions.py:105](../../interventions.py#L105) | `sti.Prep()` reimplemented (`prep_eff`/`prep_dur`/`prep_adh`, state on `HIV`) | **Leave the bare call as-is** — it constructs in both versions. See "the inherited PrEP default" below |
| [plot_dashboard.py:111](../../plot_dashboard.py#L111) | `sim.interventions['vmmc'].circumcised` | → `sim.diseases.hiv.circumcised`. Not on 017's path (`run.py` uses its own analyzer); fix on adoption |
| [experiments/015](../015_vmmc_prevalence_target/) | own copy of the VMMC subclass | Closed experiment — leave as-is, do not port |

**No repo code changes are required to run this experiment.** `run.py` selects
the VMMC class from the installed version and reads circumcision and PrEP state
through version-agnostic analyzers, so both stacks run off one unmodified
codebase and arm A stays re-runnable after the checkouts move. The actual
porting — deleting `vmmc.py`, fixing the dashboard — happens on adoption in 018.

`never_art` → `art_naive` and the removal of `post_art` do not appear in our
code; the scan found no call sites.

### The inherited PrEP default

`sti.Prep()` is called bare at [interventions.py:105](../../interventions.py#L105).
In **both** 1.5.8 and 1.5.11, `coverage=None` does not mean "off" — it falls
back to a built-in ramp:

```python
coverage = {'year': [2004, 2005, 2015, 2025], 'value': [0, 0.01, 0.5, 0.8]}
```

targeting FSW, which `make_sim` sets at 10 % of women. **Every experiment from
001 onward has run with 50 % of FSW on PrEP by 2015 and 80 % by 2025** — a
parameter nobody in this project chose, inherited from an upstream default.

1.5.11 reimplemented the mechanism (per-agent probability → stock target with
`prep_dur`/`prep_adh`/retention, state on `HIV`), so **A → B moves for two
reasons at once**: on-ART mortality *and* PrEP mechanics.

**Decision: keep the default in all three arms, and measure PrEP explicitly.**
Turning it off would give a cleaner mortality comparison but would break the
arm-A reproduction check against 016 — and a mismatch in arm A would then be
indistinguishable from a porting error, which is the one thing the control
exists to rule out. So PrEP uptake becomes metric 8, reported per arm, and if
the A → B change is large it gets isolated in 018.

`on_prep` lives on the intervention in 1.5.8 and on `HIV` in 1.5.11, and 1.5.8's
`Prep` defines no results at all — so the metric needs an analyzer that reads
whichever location exists, not a results column.

**Specifying PrEP deliberately is deferred, not dropped.** It is half the
decision question in [CLAUDE.md](../../CLAUDE.md), so it deserves its own
treatment rather than an inherited default carried forward by inertia.

**The VMMC deletion is the load-bearing claim here.** 015 established that
prevalence-target vs hazard semantics roughly *doubles* male HIV prevalence, so
if upstream's version differs from ours in any material way, it will show up in
arm B and be easy to mistake for a mortality effect. Circumcision coverage by
age is therefore a required metric, not a nice-to-have.

## Metrics

1. **Arm A vs 016's arm A** — reproduction check. Prevalence 15–49 and
   `hiv.new_deaths` should match to within seed noise. If they don't, something
   other than the version changed and the rest of the experiment is
   uninterpretable. Likewise **arm C vs 016's arm B**, which is the same check
   at the other end.
2. **AIDS deaths over time, all three arms, against UNAIDS.** The headline. 016
   measured 3 376 (arm A, 2005) and 3 856 (arm B, 2005) against UNAIDS' 11 000,
   and reached ~65 % peak-to-peak. How much of that gap does on-ART mortality
   close — biased (A→B) and unbiased (A→C)?
3. **Deaths among PLHIV split by module** (background `ss.Deaths` vs HIV),
   carried forward from 016's metric 3. 016 found the two machines competing
   for the same agents; nonzero on-ART mortality changes that balance, and arm C
   is where the competition is removed.
4. **HIV prevalence 15–49** (all, male, female), all arms. Higher HIV mortality
   removes PLHIV, so B should sit below A — which runs *against* the direction
   coverage needs. C should recover some of it (016: +1.5 points at high
   transmission). Whether C nets out above or below A is the number that matters
   for the next coverage check, and it is not obvious in advance.
5. **Circumcision coverage by age, all arms** — the VMMC equivalence check
   described above, against the SHIMS3 targets 015 used.
6. **Total population**, all arms, against UN estimates. 016 found arm B landed
   within 0.1 % of the 2015 target at default parameters; C should reproduce
   that, and any drift is a version effect on demographics.
7. **Runtime per sim**, all arms. starsim 3.5.2 includes performance work;
   the model-setup sweep in 018 needs the current number.
8. **PrEP uptake — number and share of FSW on PrEP over time**, all arms. The
   confounder check described above. If A and B diverge here, part of any
   prevalence difference is the PrEP rewrite rather than mortality, and the
   headline numbers need that caveat.

## Success criteria

- **Bump is clean and material:** arm A reproduces 016's arm A, arm C
  reproduces 016's arm B on the mortality dimension, VMMC is equivalent, and
  arm C closes a substantial share of the AIDS-death deficit. 018 shrinks from
  "diagnose the deficit" to "promote the corrected mortality and sweep model
  setup", and the bump is adopted.
- **Bump is clean and immaterial:** everything ports, nothing much moves. The
  bump is still adopted — declared and actual should agree, and `vls_coverage`
  and the new PrEP API are needed downstream regardless. 018 keeps its original
  scope, and the deficit is confirmed as structural rather than an artifact of
  ART agents being immortal.
- **Bump changes something we didn't predict:** most likely VMMC, or prevalence
  falling further than expected. This is why arm A exists. Diagnose before
  adopting; do not carry an unexplained change into a calibration.

A note on what would *not* count as success: arm C closing the deficit entirely
would be suspicious rather than reassuring, since 016 established the model is
short ~7 000 deaths a year at the peak and on-ART mortality is bounded above by
the off-ART hazard by construction. If C gets close to UNAIDS, check for
double-counting before celebrating.

## Why this matters beyond the version number

Two things in 1.5.11 are prerequisites for the decision question in
[CLAUDE.md](../../CLAUDE.md), not incidental:

- **`ART(vls_coverage=…)`** now accepts time-varying and age×sex-stratified
  suppression, so PHIA's rise in viral suppression between SHIMS rounds can
  drive the model. A treatment-cascade scenario is not expressible without it.
- **`Prep` reimplemented** with `prep_eff`/`prep_dur`/`prep_adh` and multiple
  mutually-exclusive simultaneous products — i.e. oral and long-acting PrEP as
  distinct arms. That is the other half of the decision question.

Neither is exercised in this experiment. Flagged so the bump is understood as
unblocking work, not just housekeeping.

## Not in scope

- **Promoting** 016's HIV-deleted mortality to a repo-level input under `data/`
  — 018. Arm C regenerates it into `outputs/` instead.
- Improving the mortality construction (016's option 2, apportionment). 016
  judged the crude method good enough; nothing here revisits that.
- Population size, replicate count, establishment threshold — 018.
- Any use of `vls_coverage` or the new PrEP API beyond making the model run.
- Re-running the coverage check — that follows 018.

## Reference

[notes/stisim_1.5.11_review.md](../../notes/stisim_1.5.11_review.md) — our
pre-release review of PR #580, and which points the release addressed.
