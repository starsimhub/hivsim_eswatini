# Exp 013 — Matcher comparison: what drives the network prevalence shift?

*Follows [../012_version_bump_diagnostics/SUMMARY.md](../012_version_bump_diagnostics/SUMMARY.md),
which found the updated stack's fixed network pulls HIV prevalence down but
could not say why.*

## Question

The stisim 1.5.6 network fix (#477) shifts Eswatini HIV prevalence downward,
but the change entangles **two mechanisms** that the CHANGELOG's headline
"~59 %" does not separate:

- **(a) Age structure.** The corrected matcher produces a realistic +6 yr
  male−female age gap (vs ~0 before), which — per the stisim docs —
  concentrates incidence among older men and reduces onward at-risk
  transmissions.
- **(b) Partnership volume.** The corrected matcher also *forms fewer
  partnerships*, because it adds filtering the old rank-zip did not have
  (see below). Fewer partnerships means fewer transmission opportunities.

Both lower prevalence, so the aggregate shift can't tell them apart. **This
experiment asks: how much of the prevalence shift is age structure (a) vs
partnership volume (b)?** The answer determines how we recover fit in
[../014_prior_expansion/](../014_prior_expansion/): if it's volume, we adjust
`p_pair_form`/concurrency; if it's structure, the age gap itself is doing the
epidemiological work.

## The matching methods

Each timestep, the MF network matches *looking* women to *eligible* men. Each
looking woman draws a **desired** partner age (her age + a gap sampled from
`age_diff_pars`); men enter with their **actual** ages. A `match_method`
turns those into pairs. We compare three, chosen to separate (a) from (b):

| `match_method` | How it pairs | Age gap | Partnership volume |
|---|---|---|---|
| **`sort_bisect`** *(old default, pre-1.5.6)* | Sort men by actual age and women by desired age, then **zip by rank**; trim only the non-overlapping age tails. Because both distributions have the same shape, rank-zip pairs by quantile and the age gap **collapses to ~0**. Keeps almost everyone. | ~0 (broken) | **High** |
| **`kdtree_nn`** *(strict, untapered)* | KD-tree nearest-neighbour: each woman matched to the man closest to her *desired* age; contested men go to the nearest woman. Honours the age gap, but **no older-woman taper and no `max_deviation` skip** — so it drops fewer pairs than the new default. | ~correct (+6 yr) | **Medium** |
| **`closest_age_tapered_seeking`** *(new default, 1.5.6+)* | Sorted closest-age two-pointer sweep on desired age, **plus three volume-reducing filters**: (1) older women's chance of looking tapers to zero by age 55; (2) a woman is *skipped* if her nearest free man is more than `max_deviation` (1 yr) off target; (3) pairs with a realized gap beyond `mean_gap + 3·sd` are trimmed. | ~correct (+6 yr) | **Low** |

The contrasts isolate the mechanisms:
- **`sort_bisect` → `kdtree_nn`**: age gap goes from ~0 to +6 yr at *similar*
  volume → isolates mechanism **(a)** age structure.
- **`kdtree_nn` → `closest_age_tapered_seeking`**: age gap stays ~+6 yr while
  volume drops (taper + skip + trim) → isolates mechanism **(b)** volume.

## Plan

**Model base.** starsim 3.5.0 / stisim 1.5.8, `run_sims.make_sim` unchanged.
Only `match_method` varies (passed via `network_pars`); everything else —
parameters, seeds, targets — held fixed. `sort_bisect` is still in the 1.5.8
matcher registry, so this is a same-stisim A/B, no reverting.

**Run.** 10 seeds per matcher (1985–2031, 10k agents). Per run, record:
- HIV prevalence 15-49 and new-infections 15-49 trajectories (from `hiv_epi`).
- Partnership **volume**: active MF edges over time, mean/total lifetime
  partners (total MF pairs ever formed = Σ female lifetime partners).
- Realized male−female **age gap** by female age group (2020 snapshot, via
  `PartnershipSnapshot`).

Scalars → `outputs/results.jsonl` (appended per run); trajectories →
`outputs/`; comparison figure → `figures/`.

## Success criteria

- **Primary:** a clean decomposition — a statement like "of the X-point
  prevalence drop from `sort_bisect` to the new default, ~Y points come from
  the age-gap correction and ~Z from lower partnership volume," backed by the
  `kdtree_nn` midpoint.
- **Confirms your recollection** if `sort_bisect` forms materially more
  partnerships than the new default (higher active-edge count / lifetime
  partners) and that volume gap explains a non-trivial share of the prevalence
  difference.
- **Actionable for 014:** tells us whether to recover prevalence via
  partnership-volume parameters or to treat the age structure as fixed and
  work other levers.
- A null result (volume nearly identical; shift is almost all age structure)
  is equally valid and equally useful.
