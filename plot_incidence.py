"""
Plot HIV incidence (ages 15-49) by sex with uncertainty across 10 seeds.
Reference data: UNAIDS/PHIA point estimates for 2011 and 2016 only.
Note: 2021 data is held out for validation and NOT plotted here.
"""

from pathlib import Path
import numpy as np
import sciris as sc
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from utils import run_one, summarize

REPO_DIR = Path(__file__).parent
RESULTS_DIR = REPO_DIR / 'results'
FIGURES_DIR = REPO_DIR / 'figures'
N_SEEDS = 10

# PHIA-derived HIV incidence per person-year, by sex, ages 15-49 (15-49 in 2016, 18-49 in 2011).
# Source: PHIA "Swaziland_Incidence_Data2.csv" central estimates (raw_data/, not redistributable).
# TODO(cliff): replace with the canonical PHIA report citation (table/page) once confirmed.
# Withheld for validation: 2021.
PHIA_INCIDENCE = {
    2011: {'f': 0.0314, 'm': 0.0165, 'age_note': '18–49'},
    2016: {'f': 0.0173, 'm': 0.0085, 'age_note': '15–49'},
}


def main():
    """ Run N_SEEDS sims and plot HIV incidence (15-49) by sex against PHIA reference points. """
    print(f'Running {N_SEEDS} simulations...')
    columns = ['timevec', 'hiv_epi.incidence_f_15_49', 'hiv_epi.incidence_m_15_49']
    dfs = sc.parallelize(run_one, range(1, N_SEEDS + 1), columns=columns)

    years = dfs[0]['timevec'].values
    inc_f = np.array([df['hiv_epi.incidence_f_15_49'].values for df in dfs]) * 100
    inc_m = np.array([df['hiv_epi.incidence_m_15_49'].values for df in dfs]) * 100

    f_stats = summarize(inc_f)
    m_stats = summarize(inc_m)

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.fill_between(years, f_stats['lo'], f_stats['hi'], color='#d62728', alpha=0.2, label='_nolegend_')
    ax.plot(years, f_stats['median'], color='#d62728', lw=2, label='Female (model median)')
    ax.fill_between(years, m_stats['lo'], m_stats['hi'], color='#1f77b4', alpha=0.2, label='_nolegend_')
    ax.plot(years, m_stats['median'], color='#1f77b4', lw=2, label='Male (model median)')

    for year, vals in PHIA_INCIDENCE.items():
        note = vals['age_note']
        ax.scatter(year, vals['f'] * 100, color='#d62728', s=60, zorder=5,
                   marker='o', label=f'Female PHIA ({note})' if year == 2016 else '_nolegend_')
        ax.scatter(year, vals['m'] * 100, color='#1f77b4', s=60, zorder=5,
                   marker='o', label=f'Male PHIA ({note})' if year == 2016 else '_nolegend_')

    # 2011 used a different age band (18-49 vs 15-49); flag this for the reader.
    ax.annotate('2011\n(18–49)', xy=(2011, PHIA_INCIDENCE[2011]['f'] * 100),
                xytext=(2013, PHIA_INCIDENCE[2011]['f'] * 100 + 0.3),
                fontsize=8, color='grey',
                arrowprops=dict(arrowstyle='->', color='grey', lw=0.8))

    ax.set_xlim(1990, 2031)
    ax.set_ylim(0, None)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('HIV incidence (% per year)', fontsize=12)
    ax.set_title('HIV Incidence (ages 15–49), Eswatini\nModel: 10 seeds (median ± 10th–90th percentile)', fontsize=13)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x:.1f}%'))
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = FIGURES_DIR / 'incidence_15_49_by_sex.png'
    plt.savefig(out, dpi=150)
    print(f'Saved {out}')
    plt.show()


if __name__ == '__main__':
    main()
