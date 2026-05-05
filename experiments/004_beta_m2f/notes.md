# Experiment 004 — Revisit beta_m2f to close the F:M incidence gap

## Motivation

After experiment 003b (debut F=17.5, M=18.5, SD=1, 10 seeds), male incidence tracks PHIA targets reasonably (~0.5 sim vs 0.85 target in 2016, ~0.25 vs 0.2 in 2021). Female incidence under-predicts (~0.7 sim vs 1.7 target in 2016; ~0.4 vs 1.4 in 2021). Young women (15–25) prevalence is also low. The ~2× female-to-male incidence ratio in PHIA is not reproduced.

Debut age alone cannot close this gap — points to **`beta_m2f`** (per-act male-to-female transmission probability) as the next lever.

## Plan

1. Sweep `beta_m2f` (e.g., 1×, 1.5×, 2×, 3× current value), 10 seeds each
2. Run dashboard at most promising value: `python plot_dashboard.py --label 004_beta_m2f`
3. Compare against PHIA female incidence/prevalence targets
4. Promote final figures into `experiments/004_beta_m2f/figures/`
5. Update `experiments/log.md` with rationale, results, and decision
