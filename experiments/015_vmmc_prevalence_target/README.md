# Exp 015 — VMMC prevalence-target fix (in-repo subclass)

*Runs before [../014_prior_expansion/](../014_prior_expansion/) despite the higher
number. Motivated by [../013_matcher_comparison/SUMMARY.md](../013_matcher_comparison/SUMMARY.md),
which cleared the network fix of 012's version-bump prevalence drop and named a
broken VMMC as the leading suspect.*

## Question

The exp-005 VMMC fix gave the stisim `VMMC` class **prevalence-target
semantics**: coverage was read as a cross-sectional circumcision *stock* target
per age bin (matching SHIMS3/PHIA data), topped up each step, never as a hazard.
That fix was a local edit to the editable stisim install
(`star_sim/stisim/stisim/interventions/hiv_interventions.py`). The stisim 1.5.6
rewrite of `VMMC` (#472/#477) **overwrote it**. The current upstream `VMMC`:

1. applies coverage as a **per-step hazard** on the *uncircumcised* pool
   (`n_to_circ = p × n_uncircumcised`), never correcting toward a stock; and
2. computes only an **aggregate** target, ignoring the age stratification in
   `data/vmmc_coverage.csv`.

Result: circumcision coverage overshoots to ~100% for all ages (vs SHIMS3 targets
~20–40%), visible comparing panel F of the 007 vs 012 fit dashboards. Over-
circumcising men (60% acquisition reduction) suppresses male HIV acquisition and
drags overall prevalence down — the leading candidate for the 012 drop.

**This experiment asks:** does re-implementing prevalence-target VMMC as an
in-repo subclass (a) reproduce the age-differentiated SHIMS3 coverage without
overshoot, and (b) how much does the VMMC bug alone move HIV prevalence?

## Plan

**Model base.** starsim 3.5.0 / stisim 1.5.8, `run_sims.make_sim` unchanged.

**The fix.** A `VMMCPrevalenceTarget(sti.VMMC)` subclass (this folder's
`vmmc_prevalence_target.py`) that overrides `step()` to hit an age-stratified
circumcision *prevalence* target: for each age bin, target count =
`p × (all alive males in bin)`; top up the highest-willingness uncircumcised men
to reach it; never remove. The efficacy mechanic (rel_sus reduction) is inherited
from upstream unchanged, so the only behavioural change is *who* is circumcised.
Kept in the repo (not a stisim patch) so the next stisim upgrade cannot wipe it.

**Run.** 10 seeds, broken (upstream `sti.VMMC`) vs fixed (`VMMCPrevalenceTarget`),
everything else identical. Record per arm: circumcision coverage by 5-yr age bin
over time (vs SHIMS3 targets), and HIV prevalence 15–49 (overall + male).

## Success criteria

- **Pass:** the fixed arm reproduces the age-differentiated SHIMS3 coverage
  (like 007's panel F) with no ~100% overshoot, and we get a clean number for the
  prevalence difference the VMMC bug caused.
- **Informative regardless:** even if the prevalence delta is small, it either
  confirms or rules out VMMC as the main driver of 012's drop, narrowing what
  else in the version bump to investigate.
- On success, promote `VMMCPrevalenceTarget` into `interventions.py` so all
  downstream runs (014 onward) use it.
