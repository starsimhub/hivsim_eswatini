# Experiment 011 — Network age-mixing fix

*Model-development experiment. Pauses calibration to fix a structural
component before continuing.*

## Question

The current stisim `match_pairs` uses rank-based sorting, which (we
suspect) destroys the desired partner-age-gap distribution specified in
`age_diff_pars`. Before any further calibration, we need to:

1. **Confirm** the diagnosis empirically — vary `age_diff_pars` and show
   the realized partner-age distribution doesn't change.
2. **Fix** `match_pairs` so the parameters actually drive the realized
   network.
3. **Set** `age_diff_pars` to values calibrated against DHS Eswatini 2006-07
   (and cross-checked against the Ott/Tanser 2011 KZN paper).
4. **Validate** by simulating and plotting realized partner-age
   distributions vs. data.

## Why pause calibration

Experiment 009 (coverage check) failed: model produces too few deaths
and too much late-period prevalence. Diagnosis pointed to prior-narrowness
(experiments 003-007 hard-coded several params that should be free), but
the network's age-mixing structure is also a candidate. More importantly,
the matching algorithm itself appears broken — meaning even a calibrated
prior over `age_diff_pars` would be uninterpretable. Fix the structure
first, then calibrate.

## Step 1 — Diagnostic: does `age_diff_pars` actually drive realized gaps?

Run 3 sim configurations, 5 seeds each, with sweeping `age_diff_pars` mean values:

| Config | μ for all (woman_age × risk_level) cells | Hypothesis |
|---|---|---|
| **A — zeros** | 0 | If algorithm responds: realized mean ≈ 0 |
| **B — defaults** | stisim defaults (7-8) | Baseline |
| **C — doubled** | doubled (14-16) | If algorithm responds: realized mean ≈ 14-16 |

If realized partner-age distributions are nearly identical across A, B, C
→ confirmed: the matching algorithm doesn't respect the parameter.

Output:
- `outputs/rank_test_realized_gaps.csv` — per-config aggregated gap statistics
- `outputs/rank_test_realized_gaps.png` — overlay plot of all three distributions

Cost: ~15 sims at 10 s each ≈ 2-3 min on local laptop.

## Step 2 — Fix `match_pairs`

Replace the sort-based matching with `scipy.optimize.linear_sum_assignment`
(the original code, removed in stisim commit `be3d20d`). Add a sort-based
fallback if `n > 5000` per call (unlikely to trigger at 10k agents).

Patch in `~/Dropbox/star_sim/stisim/stisim/networks.py` around the
current `match_pairs` (line ~340).

## Step 3 — Set `age_diff_pars` from data

Empirical input:
- **DHS Eswatini 2006-07** (`extract_dhs_partner_age.py`) — see
  [outputs/dhs_partner_age_summary.csv](outputs/dhs_partner_age_summary.csv).
  Mean gap by woman age bin: 8.6 (15-19), 7.2 (20-24), 7.7 (25-49 weighted).
  σ ≈ 5-6 across bins.
- **Ott / Tanser 2011** (rural KZN) — max gap ~4.5 yr, spousal partners
  have larger gaps than casual.

Starting values (updated after Ott reading; level 0 = marital → larger,
level 2 = casual → smaller):

```python
age_diff_pars = dict(
    teens=[(9, 5), (8, 5), (7, 5)],
    young=[(8, 5), (7, 5), (6, 5)],
    adult=[(8, 6), (7, 6), (6, 5)],
)
```

These go in `run_sims.py` via `network_pars`, not in upstream stisim
defaults (we'll PR those separately if appropriate).

## Step 4 — Validate

Run the same diagnostic-style test (5 seeds at the new defaults) with the
fixed match_pairs and the new `age_diff_pars`:

- Plot realized partner-age distribution from network edges, stratified
  by woman age bin
- Overlay DHS Eswatini empirical
- Iterate `(μ, σ)` if realized doesn't match within ~1 yr

## Out of scope

- **Per-relationship-type structure** (different age_diff_pars for
  marital/transitory/informal/commercial) — stisim currently has one
  distribution per (woman_age × risk_level), not per relationship type.
  EMOD's PFA has the richer structure. Defer to experiment 012 if needed.
- **Casual partner age data** — DHS v730 is biased toward marital. Could
  augment with v821a (categorical) or PHIA microdata, but for now we
  use marital data as a starting point.

## Status

`README.md` and DHS extraction script written. DHS summary stats
generated. Step 1 (diagnostic test) to follow.
