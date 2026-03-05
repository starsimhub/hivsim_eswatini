"""
HIV calibration figure: model vs UNAIDS/survey data.

Supports both single-sim (line) and multi-sim calibration stats (median + band).

Panels:
    A: HIV prevalence (15-49)
    B: New HIV infections
    C: People on ART
    D: People living with HIV
    E: HIV-related deaths
    F: Total population
"""

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from utils import set_font

LOCATION = 'eswatini'
RESULTS_DIR = 'results'
DATA_DIR = 'data'
FIGURES_DIR = 'figures'
START_YEAR = 1990
END_YEAR = 2025

# Colors
HIV_COLOR = '#2171b5'
DATA_COLOR = 'k'
BAND_ALPHA = 0.2


def plot_panel(ax, results, col, data_df=None, data_col=None, title='',
               color=HIV_COLOR, ylabel=None, si_ticks=False):
    """
    Plot a single calibration panel with model band/line and data points.

    Args:
        results: either a calib_stats DataFrame (multi-index columns with percentiles)
                 or a single-sim DataFrame (flat columns, 'timevec' for time)
    """
    multi = isinstance(results.columns, pd.MultiIndex)  # Multi-sim vs single-sim

    if multi:
        years = results.index
        mask = (years >= START_YEAR) & (years <= END_YEAR)
        yrs = years[mask]
        med = results.loc[yrs, (col, '50%')]
        lo = results.loc[yrs, (col, '10%')]
        hi = results.loc[yrs, (col, '90%')]
        ax.fill_between(yrs, lo, hi, alpha=BAND_ALPHA, color=color, linewidth=0)
        ax.plot(yrs, med, color=color, linewidth=1.5, label='Model (median)')
    else:
        years = results['timevec']
        mask = (years >= START_YEAR) & (years <= END_YEAR)
        ax.plot(years[mask], results[col][mask], color=color, linewidth=1.5, label='Model')

    if data_df is not None and data_col is not None and data_col in data_df.columns:
        d = data_df[['time', data_col]].dropna()
        d = d[(d.time >= START_YEAR) & (d.time <= END_YEAR)]
        if len(d):
            ax.scatter(d.time, d[data_col], color=DATA_COLOR, s=15, zorder=5,
                       label='Data', edgecolors='none')

    ax.set_title(title, fontsize=18, pad=6)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xlim(START_YEAR, END_YEAR)
    ax.set_ylim(bottom=0)
    if si_ticks:
        sc.SIticks(ax, axis='y')
    ax.tick_params(labelsize=13)


if __name__ == '__main__':

    # Load model results — either calibration stats (multi-sim) or single sim
    use_calib = True
    if use_calib:
        results = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_calib_stats.df')  # Multi-sim: median + bands
    else:
        results = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_sim.df')  # Single sim: line only

    data = pd.read_csv(f'{DATA_DIR}/eswatini_hiv_calib.csv')
    art_data = pd.read_csv(f'{DATA_DIR}/n_art.csv').rename(columns={'year': 'time'})

    set_font(size=18)
    fig = plt.figure(figsize=(22, 14))
    gs = GridSpec(2, 3, left=0.06, right=0.98, bottom=0.06, top=0.93,
                  wspace=0.28, hspace=0.35)

    # (A) HIV prevalence 15-49
    ax = fig.add_subplot(gs[0, 0])
    plot_panel(ax, results, 'hiv.prevalence_15_49', data, 'hiv.prevalence_15_49',
               '(A) HIV prevalence (15\u201349)', ylabel='Prevalence')
    ax.legend(fontsize=12, loc='upper right', frameon=False)

    # (B) New infections
    ax = fig.add_subplot(gs[0, 1])
    plot_panel(ax, results, 'hiv.new_infections', data, 'hiv.new_infections',
               '(B) New HIV infections', si_ticks=True)

    # (C) People on ART
    ax = fig.add_subplot(gs[0, 2])
    plot_panel(ax, results, 'hiv.n_on_art', art_data, 'n_art',
               '(C) People on ART', si_ticks=True)

    # (D) People living with HIV
    ax = fig.add_subplot(gs[1, 0])
    plot_panel(ax, results, 'hiv.n_infected', data, 'hiv.n_infected',
               '(D) People living with HIV', si_ticks=True)

    # (E) HIV-related deaths
    ax = fig.add_subplot(gs[1, 1])
    plot_panel(ax, results, 'hiv.new_deaths', data, 'hiv.new_deaths',
               '(E) HIV-related deaths', si_ticks=True)

    # (F) Total population
    ax = fig.add_subplot(gs[1, 2])
    plot_panel(ax, results, 'n_alive', data, 'n_alive',
               '(F) Total population', color='#555555', si_ticks=True)

    plt.savefig(f'{FIGURES_DIR}/fig1_hiv_calibration.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig1_hiv_calibration.png')
