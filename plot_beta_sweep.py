"""
Sweep rel_beta_f2m (female-to-male transmission relative to male-to-female)
to explore the F:M incidence gap.

  rel_beta_f2m = beta_f2m / beta_m2f
  default 0.5 -> M2F is 2x F2M per partnership-month
  lower values -> wider F>M incidence gap (closer to PHIA's ~2x ratio)

Usage:
    python plot_beta_sweep.py                      # 3 seeds, default values
    python plot_beta_sweep.py --n_seeds 5          # 5 seeds per value
    python plot_beta_sweep.py --values 0.3,0.5,1.0 # custom sweep
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from run_sims import make_sim

FIGURES_DIR = 'figures'
CALIB_DIR = 'calibration_data'

DEFAULT_VALUES = [0.25, 0.4, 0.5, 0.75, 1.0]
DEFAULT_N_SEEDS = 3
F_COLOR = '#d62728'
M_COLOR = '#1f77b4'
BAND_ALPHA = 0.15


def annual_mean(yearvec, vals):
    yi = np.floor(yearvec).astype(int)
    uy = np.unique(yi)
    return uy, np.array([np.nanmean(vals[yi == yr]) for yr in uy])


def run_sweep(values, n_seeds):
    results = {}  # rel -> {'inc_f': (n_seeds, n_years), 'inc_m': ..., 'years': ...}
    for rel in values:
        print(f'\nrel_beta_f2m = {rel}')
        inc_f_seeds, inc_m_seeds = [], []
        for seed in range(1, n_seeds + 1):
            print(f'  seed {seed}/{n_seeds}...')
            sim = make_sim(seed=seed, verbose=-1,
                           hiv_pars={'rel_beta_f2m': rel})
            sim.run()
            yv = sim.t.yearvec
            epi = sim.analyzers['hiv_epi']
            yrs, ann_f = annual_mean(yv, np.array(epi.results['incidence_f_15_49']))
            _,   ann_m = annual_mean(yv, np.array(epi.results['incidence_m_15_49']))
            inc_f_seeds.append(ann_f)
            inc_m_seeds.append(ann_m)
        results[rel] = {
            'years': yrs,
            'inc_f': np.array(inc_f_seeds) * 100,
            'inc_m': np.array(inc_m_seeds) * 100,
        }
    return results


def plot_sweep(results, inc_calib, n_seeds, out_path):
    values = sorted(results.keys())
    cmap = plt.get_cmap('viridis')
    color_map = {rel: cmap(i / max(len(values) - 1, 1))
                 for i, rel in enumerate(values)}

    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(2, 2, left=0.07, right=0.97, bottom=0.08, top=0.92,
                  wspace=0.22, hspace=0.32)

    for ax_idx, sex_key, sex_int, sex_label, color_target in [
        (0, 'inc_f', 1, 'Female', F_COLOR),
        (1, 'inc_m', 0, 'Male',   M_COLOR),
    ]:
        ax = fig.add_subplot(gs[0, ax_idx])
        for rel in values:
            r = results[rel]
            mean = np.nanmean(r[sex_key], axis=0)
            lo = np.nanpercentile(r[sex_key], 10, axis=0)
            hi = np.nanpercentile(r[sex_key], 90, axis=0)
            mask = (r['years'] >= 1990) & (r['years'] <= 2030)
            ax.plot(r['years'][mask], mean[mask], color=color_map[rel], lw=2,
                    label=f'rel_beta_f2m = {rel}')
            ax.fill_between(r['years'][mask], lo[mask], hi[mask],
                            color=color_map[rel], alpha=BAND_ALPHA, linewidth=0)
        sub = inc_calib[(inc_calib['Gender'] == sex_int) & (inc_calib['Endage'] == 49)]
        if len(sub):
            ax.errorbar(sub['Year'], sub['Incidence'],
                        yerr=[sub['Incidence'] - sub['lb'], sub['ub'] - sub['Incidence']],
                        fmt='o', color='black', capsize=4, markersize=7,
                        zorder=5, label='PHIA target')
        ax.set_title(f'{sex_label} incidence (15–49)', fontsize=13)
        ax.set_ylabel('Incidence per 100 PY', fontsize=11)
        ax.set_xlabel('Year', fontsize=11)
        ax.set_xlim(1990, 2030)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc='upper right')

    # F:M incidence ratio over time
    ax = fig.add_subplot(gs[1, 0])
    for rel in values:
        r = results[rel]
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.nanmean(r['inc_f'], axis=0) / np.nanmean(r['inc_m'], axis=0)
        mask = (r['years'] >= 1990) & (r['years'] <= 2030)
        ax.plot(r['years'][mask], ratio[mask], color=color_map[rel], lw=2,
                label=f'rel_beta_f2m = {rel}')
    # PHIA F:M ratio at survey years
    f_sub = inc_calib[(inc_calib['Gender'] == 1) & (inc_calib['Endage'] == 49)]
    m_sub = inc_calib[(inc_calib['Gender'] == 0) & (inc_calib['Endage'] == 49)]
    merged = f_sub[['Year', 'Incidence']].merge(
        m_sub[['Year', 'Incidence']], on='Year', suffixes=('_f', '_m'))
    if len(merged):
        merged['ratio'] = merged['Incidence_f'] / merged['Incidence_m']
        ax.scatter(merged['Year'], merged['ratio'], color='black', s=70,
                   zorder=5, label='PHIA F:M target')
    ax.axhline(1.0, color='grey', ls='--', lw=0.8, alpha=0.6)
    ax.set_title('F:M incidence ratio', fontsize=13)
    ax.set_ylabel('Ratio (F / M)', fontsize=11)
    ax.set_xlabel('Year', fontsize=11)
    ax.set_xlim(1990, 2030)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='upper right')

    # Summary: F & M incidence at 2016 vs target, by sweep value
    ax = fig.add_subplot(gs[1, 1])
    target_year = 2016
    f_target = inc_calib[(inc_calib['Gender'] == 1) & (inc_calib['Endage'] == 49) &
                         (inc_calib['Year'] == target_year)]
    m_target = inc_calib[(inc_calib['Gender'] == 0) & (inc_calib['Endage'] == 49) &
                         (inc_calib['Year'] == target_year)]
    f_sims, m_sims = [], []
    for rel in values:
        r = results[rel]
        idx = np.where(r['years'] == target_year)[0]
        if len(idx):
            f_sims.append(np.nanmean(r['inc_f'][:, idx[0]]))
            m_sims.append(np.nanmean(r['inc_m'][:, idx[0]]))
        else:
            f_sims.append(np.nan)
            m_sims.append(np.nan)
    x = np.arange(len(values))
    width = 0.35
    ax.bar(x - width/2, f_sims, width, color=F_COLOR, label='Female sim')
    ax.bar(x + width/2, m_sims, width, color=M_COLOR, label='Male sim')
    if len(f_target):
        ax.axhline(f_target['Incidence'].iloc[0], color=F_COLOR, ls='--',
                   lw=1.5, label=f'F target ({target_year})')
    if len(m_target):
        ax.axhline(m_target['Incidence'].iloc[0], color=M_COLOR, ls='--',
                   lw=1.5, label=f'M target ({target_year})')
    ax.set_xticks(x)
    ax.set_xticklabels([str(v) for v in values])
    ax.set_xlabel('rel_beta_f2m', fontsize=11)
    ax.set_ylabel(f'{target_year} incidence per 100 PY', fontsize=11)
    ax.set_title(f'Incidence at {target_year} vs target', fontsize=13)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(fontsize=8, loc='upper right')

    fig.suptitle(f'rel_beta_f2m sweep (n={n_seeds} seeds per value)', fontsize=15)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'\nSaved {out_path}')
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_seeds', type=int, default=DEFAULT_N_SEEDS)
    parser.add_argument('--values', type=str,
                        default=','.join(str(v) for v in DEFAULT_VALUES),
                        help='comma-separated rel_beta_f2m values')
    parser.add_argument('--outdir', type=str, default=FIGURES_DIR,
                        help='where to save the sweep figure (default: figures/ scratch)')
    args = parser.parse_args()

    values = [float(v) for v in args.values.split(',')]
    os.makedirs(args.outdir, exist_ok=True)
    inc_calib = pd.read_csv(f'{CALIB_DIR}/incidence_by_sex.csv')

    print(f'Sweeping rel_beta_f2m = {values} with {args.n_seeds} seeds each '
          f'({len(values) * args.n_seeds} sims total)')
    results = run_sweep(values, args.n_seeds)
    plot_sweep(results, inc_calib, args.n_seeds,
               f'{args.outdir}/beta_f2m_sweep.png')
