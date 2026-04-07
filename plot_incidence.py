"""
Plot HIV incidence (ages 15-49) by sex with uncertainty across 10 seeds.
Reference data: UNAIDS/PHIA point estimates for 2011 and 2016 only.
Note: 2021 data is held out for validation and NOT plotted here.
"""

import numpy as np
import sciris as sc
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from run_sims import make_sim

RESULTS_DIR = 'results'
N_SEEDS = 10

# Reference data: [year, female_incidence, male_incidence]
# 2011: ages 18-49 (PHIA); 2016: ages 15-49 (PHIA)
# 2021 withheld for validation
REF_DATA = {
    2011: {'f': 0.0314, 'm': 0.0165, 'age_note': '18–49'},
    2016: {'f': 0.0173, 'm': 0.0085, 'age_note': '15–49'},
}


def run_one(seed):
    sim = make_sim(seed=seed, verbose=-1)
    sim.run()
    df = sim.to_df(resample='year', use_years=True, sep='.')
    return df[['timevec', 'hiv_epi.incidence_f_15_49', 'hiv_epi.incidence_m_15_49']]


def main():
    print(f'Running {N_SEEDS} simulations...')
    dfs = sc.parallelize(run_one, range(1, N_SEEDS + 1))

    # Stack into arrays indexed by year
    years = dfs[0]['timevec'].values
    inc_f = np.array([df['hiv_epi.incidence_f_15_49'].values for df in dfs])
    inc_m = np.array([df['hiv_epi.incidence_m_15_49'].values for df in dfs])

    # Convert to % (multiply by 100)
    inc_f *= 100
    inc_m *= 100

    # Summary statistics across seeds
    def summarize(arr):
        return {
            'median': np.median(arr, axis=0),
            'lo': np.percentile(arr, 10, axis=0),
            'hi': np.percentile(arr, 90, axis=0),
        }

    f_stats = summarize(inc_f)
    m_stats = summarize(inc_m)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(10, 5))

    # Female model lines
    ax.fill_between(years, f_stats['lo'], f_stats['hi'],
                    color='#d62728', alpha=0.2, label='_nolegend_')
    ax.plot(years, f_stats['median'], color='#d62728', lw=2, label='Female (model median)')

    # Male model lines
    ax.fill_between(years, m_stats['lo'], m_stats['hi'],
                    color='#1f77b4', alpha=0.2, label='_nolegend_')
    ax.plot(years, m_stats['median'], color='#1f77b4', lw=2, label='Male (model median)')

    # Reference data points
    for year, vals in REF_DATA.items():
        note = vals['age_note']
        ax.scatter(year, vals['f'] * 100, color='#d62728', s=60, zorder=5,
                   marker='o', label=f'Female PHIA ({note})' if year == 2016 else '_nolegend_')
        ax.scatter(year, vals['m'] * 100, color='#1f77b4', s=60, zorder=5,
                   marker='o', label=f'Male PHIA ({note})' if year == 2016 else '_nolegend_')

    # Annotate the 2011 point to flag the age-bin difference
    ax.annotate('2011\n(18–49)', xy=(2011, REF_DATA[2011]['f'] * 100),
                xytext=(2013, REF_DATA[2011]['f'] * 100 + 0.3),
                fontsize=8, color='grey',
                arrowprops=dict(arrowstyle='->', color='grey', lw=0.8))

    # Axes formatting
    ax.set_xlim(1990, 2031)
    ax.set_ylim(0, None)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('HIV incidence (% per year)', fontsize=12)
    ax.set_title('HIV Incidence (ages 15–49), Eswatini\nModel: 10 seeds (median ± 10th–90th percentile)', fontsize=13)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x:.1f}%'))
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('figures/incidence_15_49_by_sex.png', dpi=150)
    print('Saved figures/incidence_15_49_by_sex.png')
    plt.show()


if __name__ == '__main__':
    main()
