# GitHub Issue: starsimhub/stisim

**Title:** `StructuredSexual.match_pairs` only partially honors `age_diff_pars` — per-pair age fidelity is lost within the surviving pool

**Create at:** https://github.com/starsimhub/stisim/issues/new

---

## Summary

`StructuredSexual.match_pairs` has age-feasibility cutoffs at the extremes of the eligible pool (drops young men below the youngest preferred age and old men above the oldest preferred age) — so the configured `age_diff_pars` does have *some* influence on realized partner ages. But within the surviving pool, matching is rank-based on `argsort(desired_ages)` ↔ `argsort(m_ages)`, which doesn't minimize per-pair age mismatch and largely flattens the desired-gap distribution into a function of the marginal age distributions of who's left after the cutoffs.

Net effect: realized partner-age gaps are systematically attenuated relative to the configured distribution, especially in HIV-affected populations where the male age pyramid is uneven. For Eswatini-style demographics, sweeping the configured μ from 0 → 14 yr moves realized mean by only 2.7 yr at sim end (≈19 % transmission of the parameter). For the hivsim default demographics, transmission is better but still incomplete (~57 % across the full sweep, only 27 % at typical operating values).

This makes calibrating or validating against DHS/PHIA partner-age data difficult at the `age_diff_pars` level — the realized network's age structure is dominated by who survives the extreme cutoffs, not by per-pair age preferences.

## What works (the cutoffs)

`match_pairs` lines 356-384 do gate age-feasibility at the global extremes:

```python
youngest_preferred_male_age = desired_ages[ind_f[0]]
youngest_male_age = m_ages[ind_m[0]]
if youngest_male_age < youngest_preferred_male_age:
    # drop young men below the youngest preferred age
    cutoff_index = bisect_left(m_ages[ind_m], youngest_preferred_male_age)
    ind_m = ind_m[cutoff_index:]
# ... symmetric logic at the upper bound
```

This means: as configured μ rises, more young men get dropped from the eligible pool. Useful and biologically sensible — young men without high-desirability matches don't pair as often. The "young men pair less than older men" pattern in real populations is partly captured here.

## What doesn't work (within the surviving pool)

After the cutoffs and a random sub-sampling to equalize pool sizes, matching is just rank-based:

```python
ind_m = np.argsort(m_ages, stable=True)
ind_f = np.argsort(desired_ages, stable=True)
# ... cutoffs ...
p1 = m_eligible.uids[ind_m]
p2 = f_looking[ind_f]
```

Within the surviving pool, a woman with desired age 35 sorted at rank *k* gets paired with whichever man sits at rank *k* of the surviving men — regardless of whether his age is anywhere near 35. The desired age controls her **sort position**, but not which specific man she's matched with. This flattens individual age-gap preferences into a shape dominated by the surviving marginal age distributions.

## Expected behavior

Setting `age_diff_pars=dict(teens=[(μ, σ)]*3, young=[(μ, σ)]*3, adult=[(μ, σ)]*3)` and varying μ should produce realized network partner-age gaps with mean ≈ μ for each woman age bin, within stochastic and demographic-supply tolerance.

## Actual behavior

Running 5 seeds × 3 configurations (μ ∈ {0, 7, 14}) of an Eswatini HIV sim, holding all else equal:

| Woman age bin | μ=0 | μ=7 (default) | μ=14 | DHS Eswatini 2006-07 |
|---|---|---|---|---|
| 15-24 | 1.74 | 2.09 | 5.30 | 8.6 |
| 25-34 | 0.93 | 1.24 | 2.17 | 7.4 |
| 35-49 | **−2.90** | **−2.00** | **−1.54** | 7.7 |

Key observations:

1. **Realized μ ≪ configured μ at typical operating values.** Going from μ=0 to μ=7 only moves the realized mean by 0.42 yr (Eswatini) or 1.9 yr (minimal hivsim). The bulk of the configured signal is lost.
2. **Mean realized gap is *negative* for women 35-49** across all configurations. Older women are paired with men ~2 years their junior on average — biologically implausible. The extreme cutoffs prune well-matched older men in this region (those preferred by mid-range women), leaving older women to be rank-matched against the surviving (younger) male tail.
3. **Sensitivity is nonlinear**: configured μ=0→7 transmits only 6 % through to realized gap in the Eswatini sim; configured μ=7→14 transmits 33 %. The cutoffs are threshold-like — they trigger only when the desired range diverges enough from the male age range to create a young-male surplus.

### Figure

![rank_test_realized_gaps.png](rank_test_realized_gaps.png)

Histograms (left) overlap heavily for μ=0 and μ=7 — the parameter is largely transparent in this regime. The right panel shows realized mean gap by woman age bin vs. DHS Eswatini reference. The configured μ does shift the realized distribution somewhat (especially at extreme values), but the *shape* — particularly the negative gaps for older women — doesn't match data for any setting.

## Proposed fix

Restore `linear_sum_assignment` for the matching step *within the surviving pool* (i.e. after the existing age-feasibility cutoffs run). This minimizes Σ|desired_age − available_age| pair-wise rather than relying on rank order, so each woman is matched to the man with the closest available age to her desired age.

```python
# After existing cutoffs and pool-size equalization, replace:
#   p1 = m_eligible.uids[ind_m]
#   p2 = f_looking[ind_f]
# With:
import scipy.optimize as spo
import scipy.spatial as spsp

threshold = 5000   # fallback to sort-based for very large pools (unlikely)
if len(ind_m) > threshold:
    p1 = m_eligible.uids[ind_m]
    p2 = f_looking[ind_f]
else:
    surviving_m_ages = m_ages[ind_m]
    surviving_desired = desired_ages[ind_f]
    dist_mat = spsp.distance_matrix(
        surviving_m_ages[:, np.newaxis],
        surviving_desired[:, np.newaxis],
    )
    lsa_m, lsa_f = spo.linear_sum_assignment(dist_mat)
    p1 = m_eligible.uids[ind_m[lsa_m]]
    p2 = f_looking[ind_f[lsa_f]]
```

`match_pairs` runs on **new pair formations per timestep only** (not the full population), so the surviving-pool size after cutoffs is typically 100-500 in our sims. `linear_sum_assignment` at O(n³) runs in well under a second per call at these sizes. We can submit a PR if useful.

Note: this preserves the existing age-feasibility cutoffs unchanged. The young-male-supply gatekeeping continues to operate; LSA only changes how the **surviving** pool is paired up.

## A separate (not-blocking) question

Whether the age-feasibility cutoffs should be graded (probabilistic) rather than hard thresholds at the extremes is a separate model-design decision. The current approach drops a 21-yo man if any woman in the pool has desired_age ≥ 22 — but he's a fine match for some women in the pool. A per-pair feasibility check (e.g., reject a match where |desired − actual| > threshold) would be more faithful to the input distribution. Not addressed by this issue.

## Reproduction

The full test (5 seeds × 3 configs of an Eswatini HIV sim) is in [`test_rank_matching.py`](test_rank_matching.py). Minimal stand-alone reproduction using the bundled `hivsim` wrapper (no country data required):

```python
"""Minimal repro: realized partner-age gap is attenuated relative to age_diff_pars."""
import os
os.environ["OMP_NUM_THREADS"] = "1"

import matplotlib.pyplot as plt
import numpy as np
import hivsim
import stisim as sti


def run_one(mu_gap: float, seed: int) -> np.ndarray:
    """Run a small hivsim with the given uniform mean partner-age gap; return realized gaps."""
    sexual = sti.StructuredSexual(
        age_diff_pars=dict(
            teens=[(mu_gap, 3), (mu_gap, 3), (mu_gap, 3)],
            young=[(mu_gap, 3), (mu_gap, 3), (mu_gap, 3)],
            adult=[(mu_gap, 3), (mu_gap, 3), (mu_gap, 3)],
        ),
    )
    sim = hivsim.Sim(
        n_agents=10_000, dur=25, rand_seed=seed, verbose=-1,
        networks=[sexual],
    )
    sim.run()
    e = sim.networks.structuredsexual.edges
    a1, a2 = np.asarray(e.age_p1), np.asarray(e.age_p2)
    ok = ~(np.isnan(a1) | np.isnan(a2))
    return a1[ok] - a2[ok]


configs = {"mu=0": 0, "mu=7": 7, "mu=14": 14}
seeds = range(3)
results = {label: np.concatenate([run_one(mu, s) for s in seeds])
           for label, mu in configs.items()}

for label, gaps in results.items():
    print(f"{label:>6s}: realized mean gap = {gaps.mean():.2f}, std = {gaps.std():.2f}")

bins = np.arange(-15, 31, 1)
fig, ax = plt.subplots(figsize=(8, 5))
for label, gaps in results.items():
    ax.hist(gaps, bins=bins, density=True, alpha=0.5,
            label=f"{label} -> realized {gaps.mean():.1f}")
ax.set_xlabel("Partner age - woman age (years)")
ax.set_ylabel("Density")
ax.legend()
ax.set_title("Realized partner-age gap is attenuated and partially insensitive to age_diff_pars")
plt.tight_layout()
plt.savefig("min_repro.png", dpi=120)
```

Expected output: realized means around 1.0, 2.9, 9.0 for μ ∈ {0, 7, 14}. The parameter has *some* effect but realized μ is ~60 % of configured at the extremes and ~27 % at typical operating values.

## Environment

- stisim: 1.5.4 (with local cherry-picks of `rel_sus_age` and Bellan acute pars from `fix/395-…` and `fix/396-…`)
- starsim: 3.3.3
- Python 3.14.3 / Windows 11 (also reproduces on Linux 24.04 / Python 3.12 on a 120-core Azure VM)

I have a working local fix and partial validation against DHS Eswatini 2006-07 partner-age data; happy to submit a PR if useful.
