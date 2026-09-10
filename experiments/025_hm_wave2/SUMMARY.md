# Exp 025 — The gate did not report: an unfloored CD4 decline in stisim crashed wave 2, and it was always going to

**Date:** 2026-09-10. **Model:** model-v1.3 (starsim 3.5.2 / stisim 1.5.11).
**Compute:** raccoon, 120 workers, 22 min before the crash.
**Outcome:** **failed run, root cause found and fixed.** No scientific question
this experiment asked has been answered.

**Questions.** (1) Can the pipeline separate `beta_m2f` from `s_f_young` at all?
(2) Does incidence constrain what prevalence cannot? (3) Does the NROY tighten
when σ_disc is halved, or collapse?

**Result. None of the three. The run died in the simulation phase of wave 2 of
the re-identification pre-flight**, after wave 1's 1000 points completed and
committed cleanly:

```
stisim/diseases/hiv.py:540  step_state
stisim/diseases/hiv.py:401  make_p_hiv_death
IndexError: index 6 is out of bounds for axis 0 with size 6
```

This is not a calibration failure, a prior problem, or an engine problem. It is
a latent defect in stisim, and the pre-flight was simply the first run large
enough to hit it.

## The defect

Upstream `sti.HIV` has three CD4-decline functions. Two floor the result at the
stage's end value; `acute_decline` does not:

| function | floored? |
|---|---|
| `acute_decline` | **no** — `cd4 = self.cd4[uids] - per_timestep_decline` |
| `falling_decline` | yes — `np.maximum(cd4_end, ...)` |
| `post_art_decline` | yes — `np.maximum(cd4_end, ...)` |

`acute_decline` sizes a per-timestep decrement to carry CD4 from `cd4_start` to
`cd4_end` in `acute_dur` steps, where `acute_dur = ti_latent - ti_acute` is
rounded to timesteps (upstream's own comment says so). An agent that stays acute
one step longer than that rounding assumed keeps subtracting the same decrement,
and the decrement is large relative to what remains.

Negative CD4 is then **fatal rather than merely wrong**, because the death rate
is looked up positionally:

```python
rate = pars.cd4_death_rates[np.digitize(cd4, pars.cd4_death_bins)]
```

The bins are **descending** — `[1000, 500, 350, 200, 50, 0]` — against a
6-element rate array. For descending bins `np.digitize` returns `len(bins)` only
when the value falls below the last edge. Verified across the range: **index 6
is reachable only for CD4 < 0**; every non-negative CD4 maps to 0–5. `step_state`
guards NaN CD4 ("Invalid entry for CD4") but not negative CD4.

**Demonstrated, not inferred.** Forcing the exact condition the rounding produces
— CD4 already at the latent set-point with one more decrement to apply — on real
acute agents from a live sim, and iterating:

```
upstream sti.HIV    min cd4 by extra step: 69.5 -276.8 -623.1 -969.4 ...
                    first step whose lookup raises IndexError: 2
fixed HIVCD4Floor   min cd4 by extra step: 415.8 415.8 415.8 415.8 ...
                    first step whose lookup raises IndexError: None
```

One extra step overshoots the set-point by ~350 CD4; the second crosses zero.

## Observations

**obs 1 — Not a bad region of parameter space, which is the important negative
result.** None of the seven calibration parameters touch CD4 dynamics —
`s_f_young` is `rel_sus_age`, susceptibility, not survival. The reason wave 2 hit
it and wave 1 did not is exposure: wave 2 draws inside wave 1's NROY, which
concentrates on higher transmission, so more agents are infected and more pass
through the acute phase. Nothing about the prior or the NROY needs revisiting.

**obs 2 — 024 was lucky, not safe.** 024 ran 1000 points at N = 10,000 without
crashing. 025 ran at N = 20,000, doubling agents per sim, and the exp-020 sizing
argument that raised N is what raised exposure. Any future run at this scale
would have hit it.

**obs 3 — It crashes rather than corrupting, so nothing already published is
suspect.** Every negative CD4 raises, on both the off-ART path
(`make_p_hiv_death`) and the on-ART path (`get_art_mortality_hazard`, same
positional lookup). No result produced before this fix can have been silently
distorted. **Exp 026 and the CROI abstract are unaffected.**

**obs 4 — stisim 1.6.1 has the identical defect.** Checked directly against the
1.6.1 wheel: same unfloored `acute_decline`, same 6-element arrays. Upgrading was
never a way around this.

**obs 5 — `uv sync` silently resolved the wrong model, and this nearly went
unnoticed.** `pyproject.toml` constrains `starsim>=3.5.2` / `stisim>=1.5.11` —
greater-than, not equality — so a clean sync on raccoon produced **starsim 3.6.1
/ stisim 1.6.1**, not the model-v1.3 stack that 024 (which supplied the
re-identification truth) and 026 (which consumes the NROY sample) ran on. Pinned
back before launch with `uv pip install ==` and run with `uv run --no-sync`.
Commit `27eb60a` made the repo clone-and-run by dropping `[tool.uv.sources]`;
that traded reproducibility for convenience, because the constraints were never
`==`. **Anyone cloning the repo today gets a different model than produced these
results.** Recorded in `config.yaml`.

**obs 6 — The point cache does not key on model version.** `point_key()` hashes
the parameter values plus N, stop year and seed. It does **not** include the
stisim version or the model tag, so resuming a post-fix run against this cache
would silently serve pre-fix results. Exactly the class of bug the 025 README
already caught for `N_AGENTS`. The 2000 cached points from this run must not be
reused by 027.

**obs 7 — Two documentation defects found while reading the run.**
`run.py:558-559` print `wave 1 emulates: ... / wave 2 emulates: ...`, which is a
stale label from the pre-correction design — `run.py:566` passes all three
features as one list, so all three are emulated every wave, as the README states.
The code is right and the log is misleading. Separately, this folder's
`preflight_*` outputs were **smoke-test artifacts at n = 25 with 12 NROY
samples**, not gate results, and `preflight_recovery.csv` carried an `inside`
column that reads exactly like a verdict — `log_age_gap_sd_mult` shows
`inside=False` there and it means nothing at that sample size. Renamed to
`SMOKE_*_n25.*` on closing so they cannot be mistaken for a result.

## What survives

Pulled off raccoon before the spot machine could be reclaimed:

| | |
|---|---|
| `outputs/ensemble_partial.parquet` | 2000 points consolidated, 84,000 rows — wave 1 complete plus the wave-2 points that finished before the crash |
| `outputs/hm/preflight/checkpoint.pkl` | engine checkpoint, wave 1 committed |
| `outputs/hm/preflight/run_config.json`, `log.txt` | engine-side config and log |
| `outputs/w2_crash.log` | full run log including the traceback |

Wave 1 of the pre-flight is a genuine completed wave at n = 1000. It has not been
analysed here, because a one-wave NROY does not answer the gate's question — the
degeneracy test is about what successive waves cut — and because 027 re-runs it
against a changed model anyway.

## The fix

`hiv_cd4_floor.HIVCD4Floor` — `np.maximum(cd4_end, ...)` in `acute_decline`,
matching the other two decline functions. In this repo as a subclass, not as a
patch to the stisim checkout, per CLAUDE.md; a `git pull` wiped the exp-005 VMMC
patch that way. **On by default** in `run_sims.make_sim`, since it is a fix and
not a knob — `hiv_class=sti.HIV` restores the defect for an A/B.

Tagged **model-v1.4**. Behaviourally a no-op where the floor does not bind: an
identical seeded run gives `cum_infections` 82,943 either way.

Retirement path is the one `vmmc.VMMCPrevalenceTarget` took — upstream the
one-liner, then delete the subclass and the default (exp 018).

## Acceptance

Accepted as a failed run with a diagnosed cause. The experiment's three questions
carry forward to 027 unchanged; nothing learned here bears on their answers.

## Next

**027 re-runs this experiment against model-v1.4**, with the same features, the
same σ, the same design, and the same pre-registered success criteria. Two
things must change mechanically:

1. **A fresh cache.** 027 uses its own `outputs/sims/`; the 2000 points here were
   produced by model-v1.3 and are invalid under v1.4 (obs 6).
2. **The gate reading order gains `age_gap_sd_mult`**, per exp 026's SUMMARY —
   the decision quantity is what PrEP buys AGYW 15–24 after the cascade closes
   while the coverage gap is in men 25–34, so the age-mixing parameters matter
   more to the decision than the feature set informs about them.

**Worth doing separately:** upstream the fix to stisim, and pin `==` in
`pyproject.toml` (or restore `[tool.uv.sources]`) so obs 5 cannot recur.
