# Exp 031 — Re-run 026's arms at the corrected testing rates and measure what the cascade arm actually delivers

**Question.** [026](../026_scenarios/SUMMARY.md) reported that pushing ART
coverage to 95% averts **16,525** infections over 2026-2040 — more than any
single PrEP arm — and that adding PrEP on top of the cascade still averts
**4,449**. [029](../029_cascade_audit/SUMMARY.md) then found that the `art_95`
arm **could not reach 0.95 in three of eight strata**, because awareness is a
hard ceiling on ART coverage and the model's awareness sat at 0.897 in men
25-34, 0.900 in men 15-24 and 0.832 in women 15-24. The target simply does not
fill, with no warning and no error, and 026 never measured target against
achieved. So the headline number is for a smaller intervention than its label
claims.

[030](../030_testing_rate/SUMMARY.md) then made the ceiling *lower*. Fitting the
routine testing rate to SHIMS3 gives `k_m = 0.5`, `k_f = 0.6` — the model had
been diagnosing people roughly twice as fast as the survey says. Correcting it
brings awareness down to the observed level, which means the ceiling binds
harder than it did in 029, which was itself a baseline run at the *inflated*
awareness.

This experiment re-runs 026's eight arms at the corrected rates with
`CascadeByAge` attached, and records, for every arm, every stratum and every
year: **the ART target, the awareness ceiling, and the coverage actually
achieved.** Three things follow:

1. **How far short did `art_95` land, and what does that do to 16,525?** This is
   029's outstanding measurement and the reason the experiment exists.
2. **Does the ordering survive?** 026's claim is that the cascade beats every
   single PrEP arm. If the cascade arm delivers materially less than its label,
   the margin over `prep_fsw` (8,771) narrows — possibly to nothing.
3. **Does the PrEP increment at high coverage survive?** `both − art_95` was
   4,449 (10/10 seeds positive). If `art_95` under-delivers, the arm it is
   differenced against is a weaker cascade, and the increment is being measured
   at a lower suppression level than the label implies.

**This is for the CROI abstract** whose deadline 026 recorded as 1 Oct. The
numbers 031 produces are the ones that go in it.

## Decisions taken before this run, and their cost

**The testing rates are adopted as-is, without 030's recommended refinement.**
030's Acceptance section explicitly withheld adoption pending a 2-D refinement
around (0.5, 0.6), on two grounds: one genuine ART ceiling breach (women 15-24,
2031, −1.10 pp of slack), and separability slippage between the 1-D scans and
the joint run (men 15-24 moved from −1.18 to −3.75 pp). Adam's call is to adopt
and proceed, for the deadline.

The cost being accepted: **the women 15-24 stratum will fail to fill its
baseline ART target by 1.3-2.0 pp from 2021 onward**, and it will do so in every
arm including baseline, so it is a common-mode error rather than a contrast
error. 030 obs 7 argues the breach is substantially a fit-quality problem in one
band (the model under-shoots 15-24 awareness by 3.8 pp and over-shoots 50+ by
2.8 pp) rather than a structural impossibility. **That claim is now testable
here**, because this run extends the reachability check to 2040 where 030
stopped at 2031 and the shortfall was still widening (−1.56 in 2021, −2.00 in
2031). If the 2040 shortfall is much worse than the 2031 one, the refinement was
not optional and this experiment says so.

**The model is tagged model-v1.6 at the point this runs**, per CLAUDE.md's
"tag the model before calibrating" rule: v1.5 plus the sex-split testing ramp
with defaults `test_rate_m = 0.5`, `test_rate_f = 0.6`. Changing the *defaults*
rather than passing them per-experiment is the same pattern v1.2 and v1.3 used,
and it means the abstract cites one model version.

**026 ran at model-v1.3, so this is not a clean single-variable comparison.**
Between them sit v1.4 (acute-phase CD4 floor) and v1.5 (non-negative CD4
invariant at both rate lookups). Both tags record those as near-inert — v1.4's
tag states an identical seeded run gives `cum_infections` 82943 either way, and
v1.5 only clips on paths that would otherwise have crashed. So the expected
difference between 026 and 031 is the testing split plus seeds. **That is an
expectation, not a measurement**, and the baseline arm is where it gets checked:
if 031's baseline cumulative infections 2026-2040 differ from 026's 63,162 by
more than seed noise, something other than testing moved and must be found
before any contrast is reported.

## Plan

**Arms.** 026's eight, unchanged — `baseline`, `prep_agyw`, `prep_agyw_risk`,
`prep_fsw`, `art_95`, `both`, `art_95_early`, `prep_at_high_art` — with the same
`ARM_OVERRIDES` early-cascade pair and the same seven named `CONTRASTS`. No new
arms. A testing-scale-up arm is the obvious next thing (029 obs 5: two thirds of
the men 25-34 coverage gap is people who do not know their status, and no arm in
026 moves awareness) but it is a different question and belongs in its own
folder.

**Parameters.** 024 design row 868, unchanged from 026 and 029 — the single draw
that fit 48/48 targets. `PARAM_SETS` keeps its list interface so an NROY sample
can be swapped in later. **One point means differences between arms carry seed
noise only, not parameter uncertainty. The abstract must say so and must not
report a credible interval.** This is 026's constraint and it has not changed.

**Design.** 8 arms x 10 seeds = 80 sims, 1985-2041. 026's equivalent took 1850 s
locally; `CascadeByAge` adds per-band counting on an existing analyzer pass, so
budget the same order of magnitude.

**Run script.** Copied from 026's `run.py`, not re-derived — the same reasoning
that made 026 copy `build_pars` from 025 verbatim. Three changes:

1. `CascadeByAge()` added to the analyzer list, and `cascadebyage.` added to
   `KEEP_PREFIX` so the per-band columns survive the column filter. Per
   CLAUDE.md's workflow note, `n_infected_*` per band must also survive — 016
   and 017 lost their age-stratified fit permanently this way.
2. **`test_rate_m` / `test_rate_f` added to `arm_fingerprint`.** The fingerprint
   is the cache key; 026's docstring records that omitting a spec from it let a
   changed coverage ramp reuse a stale parquet by name. The testing rates are
   now part of what defines a run and must be in the hash.
3. Stop year unchanged at 2041 so 2040 is complete.

**Measurement — the primary output.** `outputs/art_reachability.csv`, extending
030's diagnostic across all eight arms: for each arm, stratum, and year
1985-2040, the interpolated ART target, the awareness ceiling, the coverage
achieved, and the signed shortfall. 030 established the format and the annual
interpolation that matches what stisim does internally.

**Figures.**
- `scenarios.png` — 026's headline figure, regenerated: cumulative infections
  averted by arm with per-seed spread.
- `art_target_vs_achieved.png` — per stratum, the `art_95` arm's target against
  its ceiling against what it delivered, 2026-2040. This is the figure that
  answers question 1.
- `prevalence_fit_vs_phia.png` — the standard figure, from the baseline arm, via
  `standard_figures.plot_prevalence_fit`. Required of every experiment; 030
  skipped it on the argument that testing does not change who is infected, and
  that argument does not hold here because the arms change ART coverage.

## Success criteria

**Good.** The baseline reproduces 026's 63,162 within seed noise; `art_95`'s
realised per-stratum coverage is recorded and the gap to 0.95 is quantified; the
three headline numbers are re-stated as what the model actually delivered rather
than what it was asked for.

**The interesting failure, and the one to watch for.** `art_95` delivers
materially less than 026 credited it with, and the cascade's margin over
`prep_fsw` narrows or reverses. That would not invalidate 029's reasoning — it
would confirm it — but it changes the abstract's headline from "the cascade is
the larger prize" to something conditional on how much of the cascade is
actually reachable. **If that happens, the abstract says so**; the alternative
is reporting an intervention that the model could not deliver.

**The failure that would stop the experiment.** Baseline cumulative infections
differ from 026's by more than seed spread. That means something other than the
testing split changed between v1.3 and v1.6, contrary to what the tags record,
and no contrast should be reported until it is identified.

## Notes

- `run_sims.py` is on CLAUDE.md's shared-code list. Changing the defaults there
  is a behaviour change for Daniel's runs too, and is the reason it gets a tag
  and a `main` merge rather than living on a branch.
- Adopting v1.6 invalidates 028's calibration. 028's marginals were unmoved from
  the prior, so little is lost, but the next wave runs against v1.6.
- Deferred, not dropped: the ~3 pp hard-to-reach awareness floor in the 50+
  bands (030 obs 8). It will not move a PrEP-versus-cascade contrast. It would
  matter to a testing-scale-up arm, which is the next experiment.
