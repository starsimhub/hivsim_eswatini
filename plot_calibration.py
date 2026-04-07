"""
plot_calibration.py
Runs 10 simulations and generates three calibration plots:
  1. HIV incidence (15-49) by sex — model vs PHIA 2011 & 2016 with 95% CIs
  2. ART coverage — model vs PHIA survey points by sex
  3. HIV prevalence by 5-year age bin and sex — model vs SDHS 2007, SHIMS 2011, PHIA 2016

2021 data is held out for validation and NOT plotted here.
"""

import numpy as np
import sciris as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from run_sims import make_sim

N_SEEDS = 10
CALIB_DIR = 'calibration_data'
FIGURES_DIR = 'figures'


# ── Simulation ────────────────────────────────────────────────────────────────

def run_one(seed):
    sim = make_sim(seed=seed, verbose=-1)
    sim.run()
    return sim.to_df(resample='year', use_years=True, sep='.')


print(f'Running {N_SEEDS} simulations...')
dfs = [run_one(seed) for seed in range(1, N_SEEDS + 1)]
years = dfs[0]['timevec'].values


def stack(col):
    """Stack a single column across all sims into a 2D array (n_sims × n_years)."""
    return np.array([
        df[col].values if col in df.columns else np.full(len(years), np.nan)
        for df in dfs
    ])


def summarize(arr):
    return {
        'median': np.nanmedian(arr, axis=0),
        'lo':     np.nanpercentile(arr, 10, axis=0),
        'hi':     np.nanpercentile(arr, 90, axis=0),
    }


# ── Plot 1: Incidence ─────────────────────────────────────────────────────────

inc_data = pd.read_csv(f'{CALIB_DIR}/incidence_by_sex.csv')

inc_f_arr = stack('hiv_epi.incidence_f_15_49') * 100
inc_m_arr = stack('hiv_epi.incidence_m_15_49') * 100
f_stats   = summarize(inc_f_arr)
m_stats   = summarize(inc_m_arr)

fig, ax = plt.subplots(figsize=(10, 5))

# Model uncertainty bands + medians
ax.fill_between(years, f_stats['lo'], f_stats['hi'], color='#d62728', alpha=0.2)
ax.plot(years, f_stats['median'], color='#d62728', lw=2, label='Female (model median)')
ax.fill_between(years, m_stats['lo'], m_stats['hi'], color='#1f77b4', alpha=0.2)
ax.plot(years, m_stats['median'], color='#1f77b4', lw=2, label='Male (model median)')

# Survey reference points with 95% CIs
plotted_f = plotted_m = False
for _, row in inc_data.iterrows():
    if row['Startage'] not in [15, 18] or row['Endage'] != 49:
        continue
    if row['Gender'] == 1:   # Female
        color = '#d62728'
        do_label = not plotted_f
        plotted_f = True
        age_note = f"{int(row['Startage'])}–{int(row['Endage'])}"
        lbl = f'Female PHIA ({age_note})' if do_label else '_nolegend_'
    elif row['Gender'] == 0:  # Male
        color = '#1f77b4'
        do_label = not plotted_m
        plotted_m = True
        age_note = f"{int(row['Startage'])}–{int(row['Endage'])}"
        lbl = f'Male PHIA ({age_note})' if do_label else '_nolegend_'
    else:
        continue

    inc = row['Incidence']
    ax.errorbar(row['Year'], inc,
                yerr=[[inc - row['lb']], [row['ub'] - inc]],
                fmt='o', color=color, capsize=5, markersize=7, zorder=5, label=lbl)

    if row['Startage'] == 18:
        ax.annotate('18–49', xy=(row['Year'], inc),
                    xytext=(row['Year'] + 1.5, inc + 0.25),
                    fontsize=8, color='grey',
                    arrowprops=dict(arrowstyle='->', color='grey', lw=0.8))

ax.set_xlim(1990, 2031)
ax.set_ylim(0, None)
ax.set_xlabel('Year', fontsize=12)
ax.set_ylabel('HIV incidence (% per year)', fontsize=12)
ax.set_title('HIV Incidence (ages 15–49), Eswatini\n'
             'Model: 10 seeds (median ± 10th–90th pct) | PHIA 2011, 2016 with 95% CI', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/incidence_15_49_by_sex.png', dpi=150)
print('Saved: figures/incidence_15_49_by_sex.png')
plt.close()


# ── Plot 2: ART Coverage ──────────────────────────────────────────────────────

art_data = pd.read_csv(f'{CALIB_DIR}/art_coverage_by_age_sex.csv')
n_art_df = pd.read_csv('data/n_art.csv').dropna()
unaids   = pd.read_csv('data/eswatini_hiv_calib.csv')[['time', 'hiv.n_infected']].dropna()
merged   = unaids.merge(n_art_df, left_on='time', right_on='year')
merged['pct_on_art'] = merged['n_art'] / merged['hiv.n_infected'] * 100

art_model_arr = stack('hiv.p_on_art') * 100
art_stats     = summarize(art_model_arr)

fig, ax = plt.subplots(figsize=(10, 5))

ax.fill_between(years, art_stats['lo'], art_stats['hi'], color='purple', alpha=0.2)
ax.plot(years, art_stats['median'], color='purple', lw=2, label='Model overall (median)')
ax.plot(merged['time'], merged['pct_on_art'], 'k--', lw=1.5,
        label='Input: n_art.csv ÷ UNAIDS n_infected')

# PHIA survey points: age-weighted average by sex (2011, 2016 only)
calib_art = art_data[art_data['Year'] != 2021]
for gender_int, color, label in [(1, '#d62728', 'Female'), (0, '#1f77b4', 'Male')]:
    sub = calib_art[calib_art['Gender'] == gender_int]
    for year in [2011, 2016]:
        yr = sub[sub['Year'] == year]
        if len(yr) == 0:
            continue
        w = yr['Count']
        wt_mean = (yr['NationalARTPrevalence'] * w).sum() / w.sum() * 100
        wt_lb   = (yr['lb'] * w).sum() / w.sum() * 100
        wt_ub   = (yr['ub'] * w).sum() / w.sum() * 100
        offset  = 0.3 if gender_int == 1 else -0.3
        lbl = f'{label} (PHIA, age-wtd)' if year == 2011 else '_nolegend_'
        ax.errorbar(year + offset, wt_mean,
                    yerr=[[wt_mean - wt_lb], [wt_ub - wt_mean]],
                    fmt='o', color=color, capsize=5, markersize=7, zorder=5, label=lbl)

ax.set_xlim(2003, 2031)
ax.set_ylim(0, 105)
ax.set_xlabel('Year', fontsize=12)
ax.set_ylabel('% of HIV+ on ART', fontsize=12)
ax.set_title('ART Coverage, Eswatini\n'
             'Model vs. data inputs and PHIA survey points (2011, 2016)', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/art_coverage.png', dpi=150)
print('Saved: figures/art_coverage.png')
plt.close()


# ── Plot 3: Prevalence by Age Bin and Sex ─────────────────────────────────────

prev_data = pd.read_csv(f'{CALIB_DIR}/prevalence_by_age_sex.csv')

age_bins = [(15,20),(20,25),(25,30),(30,35),(35,40),(40,45),(45,50),(50,55),(55,60),(60,65)]

# Map (sex_key, lo, hi) → model output column name
col_map = {}
for sex_key, prefix_core, prefix_epi in [('f','hiv.prevalence_f_','hiv_epi.prevalence_f_'),
                                          ('m','hiv.prevalence_m_','hiv_epi.prevalence_m_')]:
    for lo, hi in age_bins:
        if lo < 35:
            col_map[(sex_key, lo, hi)] = f'{prefix_core}{lo}_{hi}'
        else:
            col_map[(sex_key, lo, hi)] = f'{prefix_epi}{lo}_{hi}'

n_bins = len(age_bins)
fig, axes = plt.subplots(n_bins, 2, figsize=(11, 22), sharex=True)

survey_markers = {2007: 's', 2011: 'o', 2016: '^'}
survey_colors  = {2007: '#555555', 2011: '#e07b00', 2016: '#007700'}

for col_idx, (sex_key, sex_label, model_color) in enumerate([
        ('f', 'Women', '#d62728'), ('m', 'Men', '#1f77b4')]):

    gender_int = 1 if sex_key == 'f' else 0

    for row_idx, (lo, hi) in enumerate(age_bins):
        ax = axes[row_idx, col_idx]

        # Model trajectory
        model_col = col_map[(sex_key, lo, hi)]
        arr   = stack(model_col) * 100
        stats = summarize(arr)
        ax.fill_between(years, stats['lo'], stats['hi'], color=model_color, alpha=0.2)
        ax.plot(years, stats['median'], color=model_color, lw=1.5)

        # Survey points with 95% CIs
        sub = prev_data[(prev_data['Gender'] == gender_int) &
                        (prev_data['start age'] == lo)]
        for _, srow in sub.iterrows():
            yr  = int(srow['Year'])
            p   = srow['NationalPrevalence'] * 100
            lb  = srow['lb'] * 100
            ub  = srow['ub'] * 100
            ax.errorbar(yr, p, yerr=[[p - lb], [ub - p]],
                        fmt=survey_markers[yr], color=survey_colors[yr],
                        capsize=3, markersize=5, zorder=5,
                        label=str(yr) if row_idx == 0 and col_idx == 0 else '_nolegend_')

        ax.set_ylim(0, 75)
        ax.set_xlim(1990, 2031)
        if row_idx == 0:
            ax.set_title(sex_label, fontsize=13, color=model_color, fontweight='bold')
        ax.set_ylabel(f'{lo}–{hi}', fontsize=8, rotation=0, labelpad=30, va='center')
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)

    axes[-1, col_idx].set_xlabel('Year', fontsize=11)

# Legend for survey years
from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0],[0], marker=survey_markers[yr], color=survey_colors[yr],
           linestyle='None', markersize=6, label=str(yr))
    for yr in [2007, 2011, 2016]
]
legend_handles += [
    plt.Rectangle((0,0),1,1, fc='grey', alpha=0.3, label='Model 10–90th pct')
]
fig.legend(handles=legend_handles, loc='upper center', ncol=4,
           fontsize=10, bbox_to_anchor=(0.5, 1.0))

fig.suptitle('HIV Prevalence by Age Bin and Sex, Eswatini\n'
             'Model (10 seeds) vs. SDHS 2007, SHIMS 2011, PHIA 2016',
             fontsize=13, y=1.02)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/prevalence_by_age_sex.png', dpi=150, bbox_inches='tight')
print('Saved: figures/prevalence_by_age_sex.png')
plt.close()

print('All done!')
