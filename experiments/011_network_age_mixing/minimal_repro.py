"""Minimal repro: realized partner-age gap is largely insensitive to age_diff_pars.

Posted to the stisim repo as part of the github_issue.md write-up. Uses only
the bundled hivsim wrapper (no country-specific demographics required).
"""
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
    print(f"{label:>6s}: realized mean gap = {gaps.mean():.2f}, std = {gaps.std():.2f}, n={len(gaps)}")

bins = np.arange(-15, 31, 1)
fig, ax = plt.subplots(figsize=(8, 5))
for label, gaps in results.items():
    ax.hist(gaps, bins=bins, density=True, alpha=0.5,
            label=f"{label} -> realized {gaps.mean():.1f}")
ax.set_xlabel("Partner age - woman age (years)")
ax.set_ylabel("Density")
ax.legend()
ax.set_title("Realized partner-age gap is largely insensitive to age_diff_pars")
plt.tight_layout()

from pathlib import Path
out_path = Path(__file__).resolve().parent / "outputs" / "minimal_repro.png"
out_path.parent.mkdir(exist_ok=True)
plt.savefig(out_path, dpi=120)
print(f"\nWrote {out_path}")
