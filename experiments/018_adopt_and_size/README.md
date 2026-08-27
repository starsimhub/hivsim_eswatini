# Exp 018 — Adopt the 1.5.11 stack, then size the model

*Combines the adoption work [017](../017_version_bump/SUMMARY.md) cleared with
the population-size and replicate-count sweep [016](../016_double_counted_mortality/SUMMARY.md)
queued. The sweep has been deferred twice — 016 moved it to 017, 017 moved it
here — because each time the model underneath it changed. It is now settled.*

## Question

Two questions, in order.

**1. What does adopting the pending changes actually do?** Four are cleared by
017 and one is a decision from 2026-08-26:

| change | authority | expected effect |
|---|---|---|
| Promote 016's HIV-deleted background mortality to `data/` | 016 accepted; 017 measured it as the only significant prevalence effect (+6.1%, z = 2.61) | known |
| Delete [vmmc.py](../../vmmc.py), use upstream `sti.VMMC` | 017 obs 2: circumcision 0.084 vs 0.084 (2005), 0.474 vs 0.475 (2021) | none |
| Fix [plot_dashboard.py:111](../../plot_dashboard.py#L111) | 017 porting table | none (not on the sim path) |
| Pin `pyproject.toml` to installed versions | declared ≠ actual since PR #2 | none |
| **Remove the inherited PrEP default** | decision, 2026-08-26 | **unmeasured** |

Only the last is unknown. `sti.Prep()` with `coverage=None` falls back to a
built-in ramp reaching 80% of FSW by 2025 — a parameter nobody chose, starting
in **2004**, a decade before PrEP had efficacy evidence. 017 measured realised
protection at ~0.67 of uninfected FSW by 2021. Removing it should *raise* late
incidence and prevalence, and the size of that matters because 2021 is both a
fitting year and the validation hold-out.

**2. How should the model be configured for calibration?** Run time, population
size, replicate count, and where the epidemic reliably establishes.
016 observation 6 and 017 observation 5 both stumbled on the establishment
problem — 2 of 10 seeds produced no epidemic at default parameters, which was
enough to corrupt a mean and send 017's first pass chasing a mechanism for an
effect that was not there. Sizing this properly is the point.

## What is already known, before running anything

**Run time:** ~145 s per sim at N = 10 000, 1985–2026, measured across 80 runs in
017 (146.4 / 147.5 / 142.1 s for its three 10-seed arms). Flat across
stisim 1.5.8 and 1.5.11.

**Population size is close to adequate already.** Expected *infected agent*
counts per PHIA target stratum, computed from 017 arm C's `PopByAgeSex` output
at N = 10 000 (54 target rows in
`calibration_data/prevalence_by_age_sex.csv` — 3 years × 2 sexes × 5-year bands
15–65):

| expected cases | rows | verdict |
|---|---|---|
| < 5 | **0** | none — no stratum is un-calibratable |
| < 10 | **2** | 2011 F 15–19 (5.1), 2007 M 60–64 (5.7) |
| ≥ 20 | 46 | comfortable |

So the answer is likely "10 000 is nearly enough, 20 000 clears the floor" —
but that calculation uses arm C at high transmission (prevalence ~0.22). In a
prior that explores lower prevalence, every count scales down with it. The sweep
tests that rather than assuming it.

**Variance is strongly heteroskedastic**, as expected near an epidemic
threshold. Between-seed SD of trajectory-mean prevalence in 017: **0.0087 at
high transmission (CV ≈ 4%)** vs **0.0520 at default parameters (CV ≈ 47%)** —
a sixfold difference driven by stochastic extinction competing with takeoff.
The replicate count therefore depends on where in parameter space we are, and
sizing for the plausible region is only meaningful once we know where that
region is. Hence part 2c.

## Plan

### Part 1 — adopt, and measure the one unknown

Apply all five changes, then a single A/B:

- **Arm `current`** = 017's arm C. **Already run** (1.5.11, HIV-deleted
  mortality, PrEP on). Re-used from `../017_version_bump/outputs/`, not re-run.
- **Arm `adopted`** = same, with the PrEP default removed.

Because everything else in 017 arm C is identical to the adopted configuration,
this A/B isolates PrEP removal exactly. 2 parameter sets × 10 seeds = 20 sims,
~5 min.

**Mortality promotion:** write `data/eswatini_deaths_hiv_deleted.csv` as a new
file and switch the default `datafolder` to use it. The original
`data/eswatini_deaths.csv` stays, so the all-cause rates remain available and
the provenance of both is legible. `mortality_construction.py` is promoted from
016 to the repo root at the same time — 017 imported it across experiment
folders, which was acceptable once and should not become a pattern.

### Part 2 — size the model

**2a. Population size.** N ∈ {5 000, 10 000, 20 000, 50 000}, 10 seeds, at
high-transmission parameters. Measure per-stratum expected counts, between-seed
CV, and run time. Confirms or refutes the pre-computed sizing above, and
measures how run time scales — the `StructuredSexual` pairing step is the term
to watch, since an O(N²) contact process would rule out N = 50 000.

**2b. Replicate count.** From 2a's CV at each N, and at 3 parameter points
(plausible, low-transmission, high-transmission), read off replicates needed
per the standard bands: CV < 5% → 3–5 replicates, 5–20% → 10–20, > 20% → 50+ or
increase N.

**2c. Establishment threshold.** A grid over the two parameters 014 identified
as controlling establishment — `beta_m2f` × `rel_init_prev` — 5 seeds per cell,
recording the fraction of seeds producing an epidemic (prevalence 15–49 > 5% at
2005). This produces the map that says which parts of the prior are
usable and where the variance explodes. It is a prerequisite for the next
coverage check, not a nice-to-have: 014's coverage failure and 017's phantom
effect both trace to draws that produced no epidemic.

## Success criteria

- **Clean:** PrEP removal moves late prevalence by a knowable amount, the sweep
  confirms N = 10 000–20 000 with a defensible replicate count, and the
  establishment map has a clear usable region. The model is then configured and
  coverage check v3 can open.
- **Awkward but useful:** the sweep says we need N = 50 000, or 20+ replicates in
  the plausible region. That is a real compute finding — it changes what method
  is feasible and hands a concrete number to `method-selection`.
- **Blocking:** run time scales worse than linearly in N, or the establishment
  region is so narrow that the prior cannot be sampled usefully. Either would
  need addressing before any calibration.

## Not in scope

- Specifying PrEP from Eswatini programme data — deferred to the decision
  analysis, where it belongs. 018 removes it and states that as an assumption.
- The residual AIDS-death deficit (~4 000/year at the peak). 017 observation 7
  points at the sex-neutral `ti_zero` clock as the likely site. That is **019**.
- Prior bounds and which parameters to open — `parameter-engineering`, after
  this.
- Any calibration run.

## Notes

- **`rel_death_f` flagged.** 017 arm D showed the model is insensitive to it
  (|z| ≤ 0.9), so it is not urgent, but it is a new-in-1.5.11 parameter with no
  documented provenance that encodes "at the same CD4, women die 26% slower" —
  a claim the seroconverter-cohort literature does not clearly support. For
  `parameter-engineering`, as a candidate for a prior spanning 1.0.
- **Two upstream reports still unfiled:** `on_prep` dead in 1.5.8/1.5.10 (017
  obs 4), and a follow-up to review point #3 on `rel_death_f`'s provenance.
