"""
Sensitivity sweep: scale condom_data values down by various factors.

Holds rel_beta_f2m=0.25 (chosen in 004) fixed.
LL pairing kept near zero (marital — broadly low across scenarios).
FSW-Client scaled with the rest (could argue for holding higher; sweep is the test).

Each scenario corresponds to a defensible interpretation of the gap between
DHS-reported condom use ("at last sex with non-marital partner") and
act-level usage during all sex acts.

Usage:
    python plot_condom_sweep.py
    python plot_condom_sweep.py --n_seeds 5
    python plot_condom_sweep.py --outdir experiments/004_beta_m2f/figures
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
DATA_DIR = 'data'
CALIB_DIR = 'calibration_data'
DEFAULT_N_SEEDS = 3
REL_BETA_F2M = 0.25  # held fixed from experiment 004 finding

F_COLOR = '#d62728'
M_COLOR = '#1f77b4'
BAND_ALPHA = 0.15

# Scenario name -> (multiplier on non-LL rows, description)
SCENARIOS = [
    ('A_default',    1.0,  'Stisim default (DHS-reported)'),
    ('B_scale_0.7',  0.7,  '0.7x — DHS overestimation correction'),
    ('C_scale_0.5',  0.5,  '0.5x — act-level vs ever-use gap'),
    ('D_scale_0.3',  0.3,  '0.3x — aggressive (typical SSA)'),
]


def scale_condom_data(base_df, factor):
    """Scale all rows except (0,0) by factor; clip LL row at floor of 0.01."""
    df = base_df.copy()
    year_cols = [c for c in df.columns if c != 'partnership']
    for i, row in df.iterrows():
        if row['partnership'] == '(0,0)':
            continue  # leave Low-Low alone (marital baseline)
        df.loc[i, year_cols] = row[year_cols].astype(float) * factor
    return df


def annual_mean(yearvec, vals):
    yi = np.floor(yearvec).astype(int)
    uy = np.unique(yi)
    return uy, np.array([np.nanmean(vals[yi == yr]) for yr in uy])


def run_scenario(scenario_label, condom_df, n_seeds):
    print(f'\nScenario: {scenario_label}')
    inc_f_seeds, inc_m_seeds = [], []
    yrs_ref = None
    for seed in range(1, n_seeds + 1):
        print(f'  seed {seed}/{n_seeds}...')
        sim = make_sim(
            seed=seed, verbose=-1,
            hiv_pars={'rel_beta_f2m': REL_BETA_F2M},
            network_pars={'condom_data': condom_df},
        )
        sim.run()
        yv = sim.t.yearvec
        epi = sim.analyzers['hiv_epi']
        yrs, ann_f = annual_mean(yv, np.array(epi.results['incidence_f_15_49']))
        _,   ann_m = annual_mean(yv, np.array(epi.results['incidence_m_15_49']))
        inc_f_seeds.append(ann_f)
        inc_m_seeds.append(ann_m)
        yrs_ref = yrs
    return {
        'years': yrs_ref,
        'inc_f': np.array(inc_f_seeds) * 100,
        'inc_m': np.array(inc_m_seeds) * 100,
    }


def plot_sweep(results, scenario_meta, condom_dfs, inc_calib, n_seeds, out_path):
    cmap = plt.get_cmap('plasma')
    n = len(scenario_meta)
    color_map = {label: cmap(0.1 + 0.75 * (i / max(n - 1, 1)))
                 for i, (label, _, _) in enumerate(scenario_meta)}

    fig = plt.figure(figsize=(16, 11))
    gs = GridSpec(3, 2, left=0.06, right=0.97, bottom=0.06, top=0.93,
                  wspace=0.22, hspace=0.42)

    # (A, B) Incidence over time, by sex
    for ax_idx, sex_key, sex_int, sex_label in [
        (0, 'inc_f', 1, 'Female'),
        (1, 'inc_m', 0, 'Male'),
    ]:
        ax = fig.add_subplot(gs[0, ax_idx])
        for label, _, desc in scenario_meta:
            r = results[label]
            mean = np.nanmean(r[sex_key], axis=0)
            lo = np.nanpercentile(r[sex_key], 10, axis=0)
            hi = np.nanpercentile(r[sex_key], 90, axis=0)
            mask = (r['years'] >= 1990) & (r['years'] <= 2030)
            ax.plot(r['years'][mask], mean[mask], color=color_map[label], lw=2,
                    label=desc)
            ax.fill_between(r['years'][mask], lo[mask], hi[mask],
                            color=color_map[label], alpha=BAND_ALPHA, linewidth=0)
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

    # (C) F:M ratio
    ax = fig.add_subplot(gs[1, 0])
    for label, _, desc in scenario_meta:
        r = results[label]
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.nanmean(r['inc_f'], axis=0) / np.nanmean(r['inc_m'], axis=0)
        mask = (r['years'] >= 1990) & (r['years'] <= 2030)
        ax.plot(r['years'][mask], ratio[mask], color=color_map[label], lw=2, label=desc)
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

    # (D) Incidence at 2016 vs target
    ax = fig.add_subplot(gs[1, 1])
    target_year = 2016
    f_target = inc_calib[(inc_calib['Gender'] == 1) & (inc_calib['Endage'] == 49) &
                         (inc_calib['Year'] == target_year)]
    m_target = inc_calib[(inc_calib['Gender'] == 0) & (inc_calib['Endage'] == 49) &
                         (inc_calib['Year'] == target_year)]
    f_sims, m_sims = [], []
    labels = []
    for label, _, desc in scenario_meta:
        r = results[label]
        idx = np.where(r['years'] == target_year)[0]
        if len(idx):
            f_sims.append(np.nanmean(r['inc_f'][:, idx[0]]))
            m_sims.append(np.nanmean(r['inc_m'][:, idx[0]]))
        else:
            f_sims.append(np.nan)
            m_sims.append(np.nan)
        labels.append(label.split('_', 1)[0])
    x = np.arange(len(labels))
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
    ax.set_xticklabels(labels)
    ax.set_xlabel('Scenario', fontsize=11)
    ax.set_ylabel(f'{target_year} incidence per 100 PY', fontsize=11)
    ax.set_title(f'Incidence at {target_year} vs target', fontsize=13)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(fontsize=8, loc='upper right')

    # (E, F) Show condom_data for a representative pairing per scenario
    # E: Mixed-risk (1,0) pairings as a proxy for "general population"
    ax = fig.add_subplot(gs[2, 0])
    for label, factor, desc in scenario_meta:
        df = condom_dfs[label]
        row = df.loc[df['partnership'] == '(1,0)']
        years = [int(c) for c in df.columns if c != 'partnership']
        vals = row.iloc[0][[str(c) for c in years]].astype(float).values
        ax.plot(years, vals, color=color_map[label], lw=2, marker='o', label=desc)
    ax.set_title('Condom use input — mixed (Med-M, Low-F) pairings', fontsize=13)
    ax.set_xlabel('Year', fontsize=11)
    ax.set_ylabel('Per-act condom use probability', fontsize=11)
    ax.set_ylim(0, 1.0)
    ax.set_xlim(1985, 2025)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='lower right')

    # F: FSW-Client
    ax = fig.add_subplot(gs[2, 1])
    for label, factor, desc in scenario_meta:
        df = condom_dfs[label]
        row = df.loc[df['partnership'] == '(fsw,client)']
        years = [int(c) for c in df.columns if c != 'partnership']
        vals = row.iloc[0][[str(c) for c in years]].astype(float).values
        ax.plot(years, vals, color=color_map[label], lw=2, marker='o', label=desc)
    ax.set_title('Condom use input — FSW–Client pairings', fontsize=13)
    ax.set_xlabel('Year', fontsize=11)
    ax.set_ylabel('Per-act condom use probability', fontsize=11)
    ax.set_ylim(0, 1.0)
    ax.set_xlim(1985, 2025)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='lower right')

    fig.suptitle(
        f'Condom-use sensitivity (rel_beta_f2m={REL_BETA_F2M}, n={n_seeds} seeds)',
        fontsize=15)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'\nSaved {out_path}')
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_seeds', type=int, default=DEFAULT_N_SEEDS)
    parser.add_argument('--outdir', type=str, default=FIGURES_DIR)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    base_df = pd.read_csv(f'{DATA_DIR}/condom_use.csv')
    inc_calib = pd.read_csv(f'{CALIB_DIR}/incidence_by_sex.csv')

    print(f'Running {len(SCENARIOS)} scenarios x {args.n_seeds} seeds '
          f'= {len(SCENARIOS) * args.n_seeds} sims')
    print(f'rel_beta_f2m = {REL_BETA_F2M} (held fixed)')

    condom_dfs = {label: scale_condom_data(base_df, factor)
                  for label, factor, _ in SCENARIOS}
    results = {}
    for label, factor, desc in SCENARIOS:
        results[label] = run_scenario(label, condom_dfs[label], args.n_seeds)

    plot_sweep(results, SCENARIOS, condom_dfs, inc_calib, args.n_seeds,
               f'{args.outdir}/condom_use_sweep.png')
