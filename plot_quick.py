"""
Quick single-sim plotting: incidence by sex, ART by age/sex, prevalence by age/sex.
Runs 1 sim (~40s) instead of 10 (~7min). Use for fast iteration.
"""

import numpy as np
import sciris as sc
import pandas as pd
import starsim as ss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_sims import make_sim

CALIB_DIR = 'calibration_data'
FIGURES_DIR = 'figures'
F_COLOR = '#d62728'
M_COLOR = '#1f77b4'

# ── Age bins for ART analyzer ────────────────────────────────────────────────
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


# ── Run single sim ───────────────────────────────────────────────────────────
print('Running 1 simulation...')
sim = make_sim(seed=1, verbose=-1, analyzers=[ARTbyAgeSex()])
sim.run()
df = sim.to_df(resample='year', use_years=True, sep='.')
years = df['timevec'].values
print('Sim done.')


# ── Plot 1: Incidence by sex ────────────────────────────────────────────────
inc_data = pd.read_csv(f'{CALIB_DIR}/incidence_by_sex.csv')

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(years, df['hiv_epi.incidence_f_15_49'] * 100, color=F_COLOR, lw=2, label='Female (model)')
ax.plot(years, df['hiv_epi.incidence_m_15_49'] * 100, color=M_COLOR, lw=2, label='Male (model)')

for _, row in inc_data.iterrows():
    if row['Startage'] not in [15, 18] or row['Endage'] != 49:
        continue
    color = F_COLOR if row['Gender'] == 1 else M_COLOR
    inc = row['Incidence']
    label = f"{'Female' if row['Gender'] == 1 else 'Male'} PHIA"
    ax.errorbar(row['Year'], inc,
                yerr=[[inc - row['lb']], [row['ub'] - inc]],
                fmt='o', color=color, capsize=5, markersize=7, zorder=5)

ax.set_xlim(1990, 2031)
ax.set_ylim(0, None)
ax.set_xlabel('Year')
ax.set_ylabel('HIV incidence (% per year)')
ax.set_title('HIV Incidence (15-49) — single sim')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/incidence_15_49_by_sex.png', dpi=150)
print('Saved incidence_15_49_by_sex.png')
plt.close()


# ── Plot 2: ART coverage by age/sex ─────────────────────────────────────────
phia = pd.read_csv(f'{CALIB_DIR}/art_coverage_by_age_sex.csv')

fig, axes = plt.subplots(1, len(ART_BINS), figsize=(5 * len(ART_BINS), 5), sharey=True)

for i, (lo, hi) in enumerate(ART_BINS):
    ax = axes[i]
    bin_label = f'[{lo},{hi})'

    for sex_key, sex_int, color, label in [('f', 1, F_COLOR, 'Female'), ('m', 0, M_COLOR, 'Male')]:
        col = f'artbyagesex.p_art_{sex_key}_{lo}_{hi}'
        vals = df[col].values * 100
        mask = years >= 2004
        ax.plot(years[mask], vals[mask], color=color, lw=2, label=label)

        sub = phia[(phia['Gender'] == sex_int) & (phia['AgeBin'] == bin_label)]
        for _, row in sub.iterrows():
            val = row['NationalARTPrevalence'] * 100
            lb = row['lb'] * 100
            ub = row['ub'] * 100
            marker = 'o' if row['Year'] != 2021 else 'D'
            alpha = 1.0 if row['Year'] != 2021 else 0.5
            ax.errorbar(row['Year'], val, yerr=[[val - lb], [ub - val]],
                        fmt=marker, color=color, capsize=4, markersize=6, zorder=5, alpha=alpha)

    ax.set_title(f'Ages {lo}\u2013{hi}', fontsize=14)
    ax.set_xlim(2004, 2028)
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('Year')
    if i == 0:
        ax.set_ylabel('% of PLHIV on ART')
        ax.legend(fontsize=10, loc='upper left')

fig.suptitle('ART Coverage by Age/Sex — single sim', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/art_coverage_by_age_sex.png', dpi=200, bbox_inches='tight')
print('Saved art_coverage_by_age_sex.png')
plt.close()


# ── Plot 3: Prevalence by age/sex ───────────────────────────────────────────
prev_data = pd.read_csv(f'{CALIB_DIR}/prevalence_by_age_sex.csv')
age_bins = [(15,20),(20,25),(25,30),(30,35),(35,40),(40,45),(45,50),(50,55),(55,60),(60,65)]

col_map = {}
for sex_key, prefix_core, prefix_epi in [('f','hiv.prevalence_f_','hiv_epi.prevalence_f_'),
                                          ('m','hiv.prevalence_m_','hiv_epi.prevalence_m_')]:
    for lo, hi in age_bins:
        if lo < 35:
            col_map[(sex_key, lo, hi)] = f'{prefix_core}{lo}_{hi}'
        else:
            col_map[(sex_key, lo, hi)] = f'{prefix_epi}{lo}_{hi}'

survey_markers = {2007: 's', 2011: 'o', 2016: '^'}
survey_colors  = {2007: '#555555', 2011: '#e07b00', 2016: '#007700'}

fig, axes = plt.subplots(len(age_bins), 2, figsize=(11, 22), sharex=True)

for col_idx, (sex_key, sex_label, model_color) in enumerate([
        ('f', 'Women', F_COLOR), ('m', 'Men', M_COLOR)]):
    gender_int = 1 if sex_key == 'f' else 0

    for row_idx, (lo, hi) in enumerate(age_bins):
        ax = axes[row_idx, col_idx]
        model_col = col_map[(sex_key, lo, hi)]
        ax.plot(years, df[model_col] * 100, color=model_color, lw=1.5)

        sub = prev_data[(prev_data['Gender'] == gender_int) & (prev_data['start age'] == lo)]
        for _, srow in sub.iterrows():
            yr = int(srow['Year'])
            p = srow['NationalPrevalence'] * 100
            lb = srow['lb'] * 100
            ub = srow['ub'] * 100
            ax.errorbar(yr, p, yerr=[[p - lb], [ub - p]],
                        fmt=survey_markers[yr], color=survey_colors[yr],
                        capsize=3, markersize=5, zorder=5,
                        label=str(yr) if row_idx == 0 and col_idx == 0 else '_nolegend_')

        ax.set_ylim(0, 75)
        ax.set_xlim(1990, 2031)
        if row_idx == 0:
            ax.set_title(sex_label, fontsize=13, color=model_color, fontweight='bold')
        ax.set_ylabel(f'{lo}\u2013{hi}', fontsize=8, rotation=0, labelpad=30, va='center')
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)

    axes[-1, col_idx].set_xlabel('Year', fontsize=11)

from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0],[0], marker=survey_markers[yr], color=survey_colors[yr],
           linestyle='None', markersize=6, label=str(yr))
    for yr in [2007, 2011, 2016]
]
fig.legend(handles=legend_handles, loc='upper center', ncol=3, fontsize=10, bbox_to_anchor=(0.5, 1.0))
fig.suptitle('HIV Prevalence by Age/Sex — single sim', fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/prevalence_by_age_sex.png', dpi=150, bbox_inches='tight')
print('Saved prevalence_by_age_sex.png')
plt.close()

print('All done!')
