"""
HIV prevalence by age and sex, compared to PHIA survey data.

Four stacked panels, one per survey year (2007, 2011, 2016, 2021).
Bars = model prevalence, diamonds = PHIA data with 95% CIs.
"""

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from utils import set_font

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

# Colors
F_COLOR = '#d46e9c'
M_COLOR = '#4a90d9'

# Sim age bins and the PHIA 5-year bins they map to
SIM_BINS = [
    (15, 20, [15]),
    (20, 25, [20]),
    (25, 30, [25]),
    (30, 35, [30]),
    (35, 50, [35, 40, 45]),
    (50, 65, [50, 55, 60]),
]


def aggregate_phia(phia, year, gender, sim_lo, sim_hi, phia_starts):
    """Average PHIA 5-year bins into a sim bin, weighted by Count."""
    rows = phia[(phia.Year == year) & (phia.Gender == gender) &
                (phia['start age'].isin(phia_starts))]
    if len(rows) == 0:
        return np.nan, np.nan, np.nan

    counts = rows['Count'].values
    prev = rows['NationalPrevalence'].values
    lb = rows['lb'].values
    ub = rows['ub'].values

    if np.any(np.isnan(counts)) or np.sum(counts) == 0:
        # Simple average if counts missing
        return np.nanmean(prev), np.nanmean(lb), np.nanmean(ub)

    w = counts / counts.sum()
    return np.average(prev, weights=w), np.average(lb, weights=w), np.average(ub, weights=w)


def plot_hiv_by_age():
    sim_df = sc.loadobj(f'{RESULTS_DIR}/eswatini_sim.df')
    phia = pd.read_csv('raw_data/SWAZILAND_nationalprevalence_all_updatedPHIA3.csv')

    survey_years = [2007, 2011, 2016, 2021]
    bin_labels = [f'{lo}\u2013{hi}' for lo, hi, _ in SIM_BINS]
    x = np.arange(len(SIM_BINS))
    width = 0.35

    set_font(size=18)
    fig = plt.figure(figsize=(10, 16))
    gs = GridSpec(len(survey_years), 1, left=0.10, right=0.95, bottom=0.05, top=0.95, hspace=0.25)

    ymax = 70

    for idx, year in enumerate(survey_years):
        ax = fig.add_subplot(gs[idx])

        # Find closest sim year
        sim_years = sim_df['timevec'].values
        sim_year = sim_years[np.argmin(np.abs(sim_years - year))]

        sim_row = sim_df[sim_df['timevec'] == sim_year].iloc[0]

        # Model bars
        f_vals, m_vals = [], []
        for lo, hi, _ in SIM_BINS:
            f_col = f'hiv.prevalence_f_{lo}_{hi}'
            m_col = f'hiv.prevalence_m_{lo}_{hi}'
            f_vals.append(sim_row.get(f_col, 0) * 100)
            m_vals.append(sim_row.get(m_col, 0) * 100)

        f_label = 'Female' if idx == 0 else None
        m_label = 'Male' if idx == 0 else None
        ax.bar(x - width/2, f_vals, width, color=F_COLOR, alpha=0.85, label=f_label)
        ax.bar(x + width/2, m_vals, width, color=M_COLOR, alpha=0.85, label=m_label)

        # PHIA data overlay
        kw = dict(marker='D', s=50, zorder=5, edgecolors='k', linewidths=0.5)
        first = True
        for j, (lo, hi, phia_starts) in enumerate(SIM_BINS):
            for gender, sex_color, offset in [(1, F_COLOR, -width/2), (0, M_COLOR, width/2)]:
                val, lb, ub = aggregate_phia(phia, year, gender, lo, hi, phia_starts)
                if np.isnan(val):
                    continue
                label = f'PHIA {year}' if first else None
                ax.scatter(x[j] + offset, val * 100, color=sex_color, label=label, **kw)
                ax.errorbar(x[j] + offset, val * 100,
                            yerr=[[(val - lb) * 100], [(ub - val) * 100]],
                            fmt='none', ecolor=sex_color, capsize=2, linewidth=1, alpha=0.7)
                first = False

        ax.set_ylim(0, ymax)
        ax.text(0.98, 0.88, str(year), transform=ax.transAxes, ha='right', va='top',
                fontsize=16, fontweight='bold')

        if idx == len(survey_years) - 1:
            ax.set_xticks(x)
            ax.set_xticklabels(bin_labels)
            ax.set_xlabel('Age group')
        else:
            ax.set_xticks(x)
            ax.set_xticklabels([])

        if idx == 0:
            ax.set_title('HIV prevalence by age and sex')
            ax.legend(frameon=False, fontsize=14, loc='upper left', ncol=3)

        if idx == len(survey_years) // 2:
            ax.set_ylabel('Prevalence (%)')

    plt.savefig(f'{FIGURES_DIR}/fig2_hiv_by_age.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig2_hiv_by_age.png')


if __name__ == '__main__':
    plot_hiv_by_age()
