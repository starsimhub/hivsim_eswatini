# Exp 014 — Diagnostic prior expansion on the updated stack

*Resumes calibration. This is 009's deferred "010" — the prior expansion —
now run on the updated starsim/stisim stack, after
[../012_version_bump_diagnostics/](../012_version_bump_diagnostics/) documented
the version bump and [../013_matcher_comparison/](../013_matcher_comparison/)
decomposed the network-driven prevalence shift.*

## Question

The prior predictive coverage check in
[../009_coverage_check/SUMMARY.md](../009_coverage_check/SUMMARY.md) failed:
only 34 % of target rows fell inside the 5–95 % envelope (20 % of AIDS-death
rows, 43 % of PHIA-prevalence rows). The diagnosis was a single class of
issue — too little mortality flow, so PLHIV accumulate and late-period
prevalence drifts up — with three suspects, all parameters that were *fixed*
rather than in the prior. Separately,
[../011_network_age_mixing/SUMMARY.md](../011_network_age_mixing/SUMMARY.md)
showed the network's age-mixing was broken (rank-matching ignored
`age_diff_pars`); the fix (#477) now ships as the default matcher in stisim
1.5.6+.

This experiment asks: **with the network fixed and the three suspect
parameters opened into the prior, does the coverage check now pass — or does
the deaths/prevalence miss persist and point to deeper structural work?**

## Plan

**Model base.** starsim 3.5.0 + stisim 1.5.8 (latest releases). The #477
network fix ships as the default `closest_age_tapered_seeking` matcher — no
local patch to carry. `age_diff_pars` set to the DHS Eswatini 2006-07 values,
after 012 confirms the shipped matcher reproduces them (011 validated a
*different* algorithm, `linear_sum_assignment`).

**Prior changes vs. 009** (following 009's three suspects):
1. Add `rel_init_prev` to the prior (currently hard-coded at 0.2). Range
   ~0.05–0.5.
2. Widen `rel_dur_on_art` upper bound (currently 1–20) to ~1–50, or switch
   to log-scale.
3. Expose one HIV-mortality / progression parameter that is currently fixed
   at stisim defaults. Candidate to confirm during config: the CD4-decline /
   time-to-death scale in `stisim/diseases/hiv.py`. Add a single multiplier
   to the prior rather than freeing the whole natural-history block.

All other prior ranges, the target set (frozen in exp 008), and the
observation model are **unchanged** from 009, so the coverage result is
attributable to the network fix plus these three prior changes.

**Run.** Prior predictive coverage check, same machinery as 009: draw from
the prior, simulate, compare to PHIA prevalence + UNAIDS deaths inside the
5–95 % envelope. 50 sims for a first pass on raccoon (009 ran 50 in ~38 s);
scale to ~1000 if a more thorough envelope is wanted. Single replicate.

## Success criteria

- **Pass:** a clear majority of target rows fall inside the envelope, and in
  particular the deaths peak (1995–2005) and post-2010 plateau are covered,
  and 2016 female prevalence no longer overshoots. This would mean 009's
  failure was prior-narrowness plus the network bug — both now fixed — and we
  proceed to `method-selection`.
- **Partial / fail:** coverage improves but deaths still undershoot the
  post-2010 plateau or prevalence still drifts. This would isolate which of
  the three suspects (or the network) mattered, and indicate the model needs
  structural work on HIV mortality before calibration is meaningful.
- A clean failure with a clear per-suspect diagnosis is a valid, useful
  result.
