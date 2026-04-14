"""
Plot ART coverage by age group and sex: model vs PHIA data.

Runs 10 simulations with a custom analyzer that tracks the proportion of
PLHIV on ART within each (age bin x sex) stratum at every timestep.
"""

import numpy as np
import sciris as sc
import pandas as pd
import matplotlib.pyplot as plt
import starsim as ss
from run_sims import make_sim

N_SEEDS = 10
CALIB_DIR = 'calibration_data'
FIGURES_DIR = 'figures'

# Age bins matching the PHIA/ART coverage data
ART_BINS = [(15, 25), (25, 35), (35, 45)]


class ARTbyAgeSex(ss.Analyzer):
    """Track proportion of PLHIV on ART by age bin and sex."""

    def init_results(self):
        super().init_results()
        results = []
        for sex in ['f', 'm']:
            for lo, hi in ART_BINS:
                results.append(ss.Result(f'p_art_{sex}_{lo}_{hi}', dtype=float, scale=False))
        self.define_results(*results)

    def step(self):
        sim = self.sim
        ti = self.ti
        hiv = sim.diseases.hiv
        ppl = sim.people
        alive = ppl.alive

        for sex_key, sex_bool in [('f', ppl.female), ('m', ppl.male)]:
            for lo, hi in ART_BINS:
                in_bin = alive & sex_bool & (ppl.age >= lo) & (ppl.age < hi)
                n_infected = (in_bin & hiv.infected).count()
                if n_infected > 0:
                    n_on_art = (in_bin & hiv.on_art).count()
                    self.results[f'p_art_{sex_key}_{lo}_{hi}'][ti] = n_on_art / n_infected
                else:
                    self.results[f'p_art_{sex_key}_{lo}_{hi}'][ti] = np.nan


# ── Run simulations ──────────────────────────────────────────────────────────

def run_one(seed):
    analyzer = ARTbyAgeSex()
    sim = make_sim(seed=seed, verbose=-1, analyzers=[analyzer])
    sim.run()
    return sim.to_df(resample='year', use_years=True, sep='.')


print(f'Running {N_SEEDS} simulations with ART-by-age-sex analyzer...')
dfs = [run_one(seed) for seed in range(1, N_SEEDS + 1)]
years = dfs[0]['timevec'].values


def stack(col):
    return np.array([
        df[col].values if col in df.columns else np.full(len(years), np.nan)
        for df in dfs
    ])


def summarize(arr):
    return {
        'median': np.nanmedian(arr, axis=0),
        'lo': np.nanpercentile(arr, 10, axis=0),
        'hi': np.nanpercentile(arr, 90, axis=0),
    }


# ── Load PHIA data ───────────────────────────────────────────────────────────

phia = pd.read_csv(f'{CALIB_DIR}/art_coverage_by_age_sex.csv')

# ── Plot: one panel per age bin, lines for M/F ──────────────────────────────

F_COLOR = '#d62728'
M_COLOR = '#1f77b4'

fig, axes = plt.subplots(1, len(ART_BINS), figsize=(5 * len(ART_BINS), 5), sharey=True)
if len(ART_BINS) == 1:
    axes = [axes]

for i, (lo, hi) in enumerate(ART_BINS):
    ax = axes[i]
    bin_label = f'[{lo},{hi})'

    for sex_key, sex_int, color, label in [('f', 1, F_COLOR, 'Female'), ('m', 0, M_COLOR, 'Male')]:
        col = f'artbyagesex.p_art_{sex_key}_{lo}_{hi}'
        arr = stack(col) * 100
        stats = summarize(arr)

        # Model band + median
        mask = years >= 2004
        ax.fill_between(years[mask], stats['lo'][mask], stats['hi'][mask],
                         color=color, alpha=0.15)
        ax.plot(years[mask], stats['median'][mask], color=color, lw=2, label=label)

        # PHIA data points
        sub = phia[(phia['Gender'] == sex_int) & (phia['AgeBin'] == bin_label)]
        for _, row in sub.iterrows():
            val = row['NationalARTPrevalence'] * 100
            lb = row['lb'] * 100
            ub = row['ub'] * 100
            marker = 'o' if row['Year'] != 2021 else 'D'
            alpha = 1.0 if row['Year'] != 2021 else 0.5
            ax.errorbar(row['Year'], val,
                        yerr=[[val - lb], [ub - val]],
                        fmt=marker, color=color, capsize=4, markersize=6,
                        zorder=5, alpha=alpha)

    ax.set_title(f'Ages {lo}\u2013{hi}', fontsize=14)
    ax.set_xlim(2004, 2028)
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('Year', fontsize=11)
    if i == 0:
        ax.set_ylabel('% of PLHIV on ART', fontsize=12)
        ax.legend(fontsize=10, loc='upper left')

# Add note about 2021 diamonds
axes[-1].annotate('Diamonds = 2021\n(validation, not calibrated)',
                  xy=(0.98, 0.02), xycoords='axes fraction',
                  fontsize=8, ha='right', va='bottom', color='grey')

fig.suptitle('ART Coverage by Age Group and Sex, Eswatini\n'
             'Model (10 seeds, median + 10\u201390th pct) vs PHIA survey data',
             fontsize=14, y=1.02)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/art_coverage_by_age_sex.png', dpi=200, bbox_inches='tight')
print(f'Saved {FIGURES_DIR}/art_coverage_by_age_sex.png')
plt.close()
print('Done!')
