# Exp 012 — Version-bump diagnostics on the updated stack

*A re-baseline. Before resuming calibration in
[../013_prior_expansion/](../013_prior_expansion/), document how the model
behaves on the updated dependency stack — especially the shipped network
fix — so we know what changed and against which baseline 013 is calibrating.*

## Question

Between experiment 011 and now, both dependencies advanced several releases:
**starsim 3.3.3 → 3.5.0** and **stisim 1.5.5 → 1.5.8**. The headline change is
that the age-mixing fix (#477, confirmed in
[../011_network_age_mixing/SUMMARY.md](../011_network_age_mixing/SUMMARY.md))
now ships as the **default** matcher — but via `closest_age_tapered_seeking`,
a *different* algorithm than the `linear_sum_assignment` we validated in 011.
The stisim 1.5.6 CHANGELOG warns this shifts HIV prevalence by **~59 % in
Eswatini at default settings** and that "downstream calibrated models will
need recalibration."

This experiment asks: **on the updated stack, does the model still run
end-to-end, and how dramatically does behavior change vs. the last known
(009-era) baseline — in HIV prevalence, incidence, deaths, and the realized
partner-age structure?** It is a documentation / re-baseline pass, not a fit.

## Plan

**Model base.** starsim 3.5.0 + stisim 1.5.8, current `run_sims.make_sim`
(unchanged from the 1.5.5 API migration in commit `4bce62a`). No prior
changes, no parameter changes — this characterises the *default* new-stack
behavior.

**Run.** `run.py` runs the model for 10 seeds (1985–2031, 10k agents) and
regenerates the two standard diagnostics via `plot_dashboard.py`:
- **Fit dashboard** — incidence by sex, prevalence by age (F/M) at survey
  years vs PHIA, ART/VMMC coverage.
- **Network dashboard** — debut CDF, risk-group composition, age-pairing
  heatmap, lifetime-partnership distribution.

Collected multi-seed data is saved to `outputs/` so figures can be
regenerated without re-running. Comparison baseline: 009's SUMMARY numbers
and the pre-update dashboards.

## Success criteria

- **Primary:** the model runs clean on the new stack and both dashboards
  render — we can state, with numbers, how prevalence / incidence / deaths /
  age-mixing moved vs. the 009-era baseline.
- **Expected (not a failure):** a large prevalence shift (~59 % order) from
  the network fix. The point is to quantify magnitude and direction, and to
  confirm the shipped `closest_age_tapered_seeking` matcher reproduces the DHS
  partner-age gaps (the check 011 did on the other algorithm).
- **Would flag a problem:** a crash / API break, or a change so extreme
  (e.g. prevalence collapse to ~0 or runaway) that it signals a regression
  rather than the expected recalibration-scale shift.
