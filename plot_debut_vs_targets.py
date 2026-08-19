"""
Run 10 sims at current debut settings (F=17.5, M=18.5, SD=1) and plot
sim output vs calibration targets:
  1. HIV incidence (15-49) by sex, over time
  2. HIV prevalence by 5-year age bin and sex (at survey years)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_sims import make_sim

FIGURES_DIR = 'figures'
DATA_DIR = 'data'
N_SEEDS = 10

# Age bins available in calibration targets
AGE_BINS_YOUNG = [(15, 20), (20, 25), (25, 30), (30, 35)]   # from hiv module (hiv.prevalence_{sex}_{lo}_{hi})
AGE_BINS_OLD = [(35, 40), (40, 45), (45, 50), (50, 55), (55, 60), (60, 65)]  # from hiv_epi analyzer
ALL_BINS = AGE_BINS_YOUNG + AGE_BINS_OLD


def annual_mean(yearvec, vals):
    years_int = np.floor(yearvec).astype(int)
    unique_years = np.unique(years_int)
    return unique_years, np.array([np.mean(vals[years_int == yr]) for yr in unique_years])


def annual_last(yearvec, vals):
    """Take the last monthly value within each calendar year (snapshot of state)."""
    years_int = np.floor(yearvec).astype(int)
    unique_years = np.unique(years_int)
    out = []
    for yr in unique_years:
        mask = years_int == yr
        out.append(vals[mask][-1])
    return unique_years, np.array(out)


# ── Run sims ─────────────────────────────────────────────────────────────────
print(f'Running {N_SEEDS} sims with debut F=17.5, M=18.5 (SD=1)...')
sims = []
for seed in range(1, N_SEEDS + 1):
    print(f'  seed {seed}...')
    sim = make_sim(seed=seed, verbose=-1)
    sim.run()
    sims.append(sim)
print('Done.')

# ── Collect data ─────────────────────────────────────────────────────────────
yearvec = sims[0].t.yearvec

# Incidence: annual mean per seed
inc_f_seeds = []
inc_m_seeds = []
for sim in sims:
    epi = sim.analyzers['hiv_epi']
    yrs, f_ann = annual_mean(yearvec, np.array(epi.results['incidence_f_15_49']))
    _,   m_ann = annual_mean(yearvec, np.array(epi.results['incidence_m_15_49']))
    inc_f_seeds.append(f_ann)
    inc_m_seeds.append(m_ann)
inc_f = np.array(inc_f_seeds) * 100
inc_m = np.array(inc_m_seeds) * 100

# Prevalence by age bin & sex: annual last-of-year per seed
prev_seeds = {('f', lo, hi): [] for (lo, hi) in ALL_BINS}
prev_seeds.update({('m', lo, hi): [] for (lo, hi) in ALL_BINS})

for sim in sims:
    hiv = sim.results.hiv
    epi = sim.analyzers['hiv_epi']
    for sex in ['f', 'm']:
        for lo, hi in AGE_BINS_YOUNG:
            # These live on sim.results.hiv (not hiv_epi)
            key = f'prevalence_{sex}_{lo}_{hi}'
            vals = np.array(hiv[key])
            _, ann = annual_last(yearvec, vals)
            prev_seeds[(sex, lo, hi)].append(ann)
        for lo, hi in AGE_BINS_OLD:
            key = f'prevalence_{sex}_{lo}_{hi}'
            vals = np.array(epi.results[key])
            _, ann = annual_last(yearvec, vals)
            prev_seeds[(sex, lo, hi)].append(ann)

prev_years = yrs  # same for all
prev = {k: np.array(v) for k, v in prev_seeds.items()}

# ── Load calibration targets ─────────────────────────────────────────────────
targets = pd.read_csv(f'{DATA_DIR}/eswatini_hiv_calib.csv')

# ── Plot 1: Incidence by sex with targets ────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(9, 5))

for inc_arr, color, sex_label, target_col in [
    (inc_f, '#d62728', 'Female', 'hiv_epi.incidence_f_15_49'),
    (inc_m, '#1f77b4', 'Male', 'hiv_epi.incidence_m_15_49'),
]:
    mean = np.mean(inc_arr, axis=0)
    lo = np.percentile(inc_arr, 10, axis=0)
    hi = np.percentile(inc_arr, 90, axis=0)
    ax1.plot(yrs, mean, color=color, lw=2, label=f'{sex_label} — sim')
    ax1.fill_between(yrs, lo, hi, color=color, alpha=0.2)

    # Overlay targets
    tgt = targets[['time', target_col]].dropna()
    ax1.scatter(tgt['time'], tgt[target_col] * 100, color=color, s=50,
                edgecolor='black', linewidth=1, zorder=5, label=f'{sex_label} — target')

ax1.set_xlabel('Year', fontsize=12)
ax1.set_ylabel('HIV incidence per 100 PY (15–49)', fontsize=12)
ax1.set_title(f'HIV Incidence vs Calibration Targets (debut F=17.5, M=18.5, n={N_SEEDS} seeds)',
              fontsize=13)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(1990, 2030)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/incidence_vs_targets.png', dpi=150, bbox_inches='tight')
print(f'Saved {FIGURES_DIR}/incidence_vs_targets.png')
plt.close()

# ── Plot 2: Prevalence by age bin & sex with targets ─────────────────────────
n_bins = len(ALL_BINS)
ncols = 5
nrows = int(np.ceil(n_bins / ncols))

fig2, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows),
                          sharex=True, sharey=True)
axes = axes.flatten()

for i, (lo, hi) in enumerate(ALL_BINS):
    ax = axes[i]

    for sex, color, sex_label in [('f', '#d62728', 'Female'), ('m', '#1f77b4', 'Male')]:
        arr = prev[(sex, lo, hi)] * 100
        mean = np.mean(arr, axis=0)
        p10 = np.percentile(arr, 10, axis=0)
        p90 = np.percentile(arr, 90, axis=0)
        ax.plot(prev_years, mean, color=color, lw=2, label=f'{sex_label} — sim')
        ax.fill_between(prev_years, p10, p90, color=color, alpha=0.2)

        # Overlay targets — check both hiv.* and hiv_epi.* prefixes
        for prefix in ['hiv', 'hiv_epi']:
            col = f'{prefix}.prevalence_{sex}_{lo}_{hi}'
            if col in targets.columns:
                tgt = targets[['time', col]].dropna()
                if len(tgt):
                    ax.scatter(tgt['time'], tgt[col] * 100, color=color, s=40,
                               edgecolor='black', linewidth=1, zorder=5,
                               label=f'{sex_label} — target')
                break

    ax.set_title(f'Ages {lo}–{hi}', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1990, 2030)
    if i % ncols == 0:
        ax.set_ylabel('Prevalence (%)', fontsize=11)
    if i >= n_bins - ncols:
        ax.set_xlabel('Year', fontsize=11)

# Hide unused axes
for j in range(n_bins, len(axes)):
    axes[j].axis('off')

# Single legend
handles, labels = axes[0].get_legend_handles_labels()
fig2.legend(handles, labels, loc='upper center', ncol=4, fontsize=11,
            bbox_to_anchor=(0.5, 1.02))

fig2.suptitle(f'HIV Prevalence by Age & Sex vs Targets (debut F=17.5, M=18.5, n={N_SEEDS} seeds)',
              fontsize=14, y=1.05)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/prevalence_by_age_sex_vs_targets.png', dpi=150, bbox_inches='tight')
print(f'Saved {FIGURES_DIR}/prevalence_by_age_sex_vs_targets.png')
plt.close()
