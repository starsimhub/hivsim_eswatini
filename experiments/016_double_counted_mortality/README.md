# Exp 016 — Double-counted HIV mortality: does the background mortality data already contain AIDS deaths?

*Opened after a code read prompted by a question about on-ART mortality. Not
the experiment originally planned for this slot — the population-size and
replicate sweep moves to 017.*

## Question

`data/eswatini_deaths.csv` supplies the background mortality rates that stisim
hands to `ss.Deaths`, which kills agents regardless of HIV status. Mean annual
all-cause mortality for ages 30–44 in that file:

| Year | Female | Male |
|---|---|---|
| 1985 | 0.0047 | 0.0059 |
| 1995 | 0.0109 | 0.0125 |
| **2005** | **0.0259** | **0.0293** |
| 2015 | 0.0126 | 0.0168 |
| 2025 | 0.0045 | 0.0073 |

Adult mortality rises 5.5× and peaks in 2005, then returns to baseline. **That
hump is the AIDS epidemic.** These are all-cause rates, so AIDS deaths are
already inside the background mortality module — and the HIV module then kills
additional agents on top, via `p_hiv_death` and the `ti_zero` AIDS pathway.

If that reading is right, the model double-counts HIV mortality, with two
consequences that point straight at the failed coverage checks:

1. **Prevalence is suppressed.** PLHIV are subject to AIDS-inflated background
   mortality *and* HIV mortality, so they cannot accumulate.
   [014](../014_prior_expansion/SUMMARY.md) found the ensemble below the
   observation in 85 of 89 target rows.
2. **The deaths target is undercounted.** `hiv.new_deaths` records only deaths
   executed by the HIV module. AIDS deaths delivered by `ss.Deaths` are
   invisible to it. 014's simulated median was ~1 000 deaths/year against an
   observed ~11 000.

**This experiment asks: how large is the effect, and does correcting it move
prevalence and deaths toward the targets?**

014's SUMMARY attributes the coverage failure to prior narrowness and the
epidemic-establishment threshold. If double-counting turns out to be the
dominant driver, that attribution needs correcting — which requires a number,
not just a fix.

## Plan

**Design: A/B at fixed parameters**, the same shape as
[015](../015_vmmc_prevalence_target/SUMMARY.md). Same seeds, same parameters,
only the background mortality data differs. This isolates the effect rather
than confounding it with a prior change.

- **Arm A (current):** `data/eswatini_deaths.csv` as shipped.
- **Arm B (HIV-deleted):** background rates with the AIDS hump removed.

**Parameters.** 009's fixed settings, which are documented and were used for the
first coverage check — so the result also speaks to why 014 sits below 009.
10 seeds per arm.

**Constructing the HIV-deleted rates.** Checked what the repo actually holds
before choosing. All-cause mortality *is* age/sex-disaggregated
(`data/eswatini_deaths.csv`: Time × Sex × AgeStart, decadal). The quantity that
would be subtracted is not: in `data/eswatini_hiv_calib.csv`, AIDS deaths
(`hiv.new_deaths`) is a single national annual column, 1990–2024, both sexes,
no age breakdown — the only undisaggregated quantity in a 43-column file where
prevalence is split by age and sex and PLHIV by age band. So the direct
subtraction `all_cause[age,sex,year] − aids[age,sex,year]` cannot be computed
from what is here.

Three constructions, in increasing order of cost:

1. **Interpolation (this experiment).** Interpolate age/sex rates between 1985
   and 2025, replacing 1995/2005/2015. Assumes both endpoints are AIDS-free
   (2025 still carries residual AIDS mortality) and ignores secular non-AIDS
   improvement. Crude, but sufficient to size the effect, which is all this
   experiment needs.
2. **Apportionment.** Split the national AIDS-death total across strata using
   UNAIDS' own PLHIV distribution as weights — `hiv_epi.n_infected_0_14` and
   `hiv_epi.n_infected_15_100` form a clean child/adult partition, and
   `hiv.n_infected_f` gives a sex split. Assumes AIDS mortality is proportional
   to PLHIV within a stratum, which ignores ART coverage varying by age. Coarse
   (2 age groups × 2 sexes, not 5-year bins) and needs population denominators
   by age/sex/year to convert counts to rates.
3. **Sourced.** HIV-deleted life tables, or age/sex-disaggregated AIDS deaths
   from UNAIDS AIDSinfo. Requires going outside the repo.

**Deliberately using (1) here.** Options 2 and 3 are real data-construction
work, and doing that before knowing whether the effect is 2 % or 50 % is
premature. If this experiment shows the effect is material, the fix gets its own
experiment using (2) or (3).

**Provenance.** No script in this repo generates `data/eswatini_deaths.csv` and
its source is undocumented. Establishing where it came from (UN WPP release and
variant, most likely) is part of this experiment — a calibration input this
consequential should not have unknown provenance.

**Documentation correction.** CLAUDE.md and
[008's README](../008_calibration_targets/README.md) both cite
`data/eswatini_deaths.csv` as the source of the UNAIDS deaths target. It is
actually `data/eswatini_hiv_calib.csv` (see
[008's run.py](../008_calibration_targets/run.py) line 74). CLAUDE.md has been
corrected; 008 is closed and left as written, so the correction is recorded
here.

## Metrics

1. HIV prevalence 15–49 (all, male, female) over time, both arms.
2. `hiv.new_deaths` over time, both arms, against the UNAIDS target.
3. **Deaths among PLHIV attributed to the background module vs the HIV module.**
   This is the direct measurement of the invisible-AIDS-deaths claim, and the
   single most diagnostic number in the experiment.
4. **Total population over time, both arms, against UN estimates.** The key
   sanity check: if the current double-counted model already tracks population,
   removing the inflated mortality will overshoot it — which would mean the
   double-counting is partly compensating for something else and the fix is not
   as simple as deleting the hump.

## Success criteria

- **Confirmed and material:** arm B raises prevalence and deaths substantially
  toward the targets, and metric 3 shows a meaningful number of AIDS deaths
  currently going uncounted. This reframes 009 and 014, and the corrected
  mortality goes into the model before any further coverage check.
- **Confirmed but small:** the effect exists but is a few percent. Worth fixing
  for correctness, but 014's prior-narrowness diagnosis stands, and the route
  through 017/018 is unchanged.
- **Refuted:** prevalence and deaths barely move, or population overshoots badly
  in arm B. Either the double-counting reading is wrong or something else
  compensates for it. A clean refutation is a useful result and cheap to obtain.

## Downstream note

Metric 3 raises a targets question that belongs to `likelihood-design` rather
than here: if AIDS deaths can be delivered by either module, is `hiv.new_deaths`
the right model output to compare against UNAIDS AIDS deaths at all, or should
the target be all deaths among PLHIV? Flagged, not resolved in this experiment.

## Not in scope

Population size, replicate count and the epidemic-establishment threshold — that
sweep moves to **017**, and is better run on a model whose mortality is
understood.
