"""
Experiment 03: Coverage check against PHIA survey data.

Changes from experiment 02:
  - p_pair_form added to prior (0.30–0.70) to probe timing–magnitude tradeoff
  - Coverage plots now show PHIA survey measurements (orange diamonds) as the
    PRIMARY calibration targets, with UNAIDS modeled estimates (grey dots) as
    secondary context
  - New panels: female/male prevalence 20–25, female/male incidence 15–49

Run from the repo root:
    uv run python experiments/03_coverage_phia/run.py
"""

import os
import sys

os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                  NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1')

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import sciris as sc
import starsim as ss
import stisim as sti

from interventions import make_interventions
from analyzers import hiv_epi

OUTDIR = os.path.join(os.path.dirname(__file__), 'results')
FIGDIR = os.path.join(os.path.dirname(__file__), 'figures')
DATA_DIR = os.path.join(ROOT, 'data')

# --- Prior --------------------------------------------------------------------

PRIOR = dict(
    beta_m2f       = (0.002, 0.020),
    eff_condom     = (0.50,  0.90),
    rel_init_prev  = (0.05,  0.50),
    rel_dur_on_art = (1.0,   20.0),
    prop_f0        = (0.55,  0.90),
    prop_m0        = (0.55,  0.80),
    m1_conc        = (0.05,  0.30),
    p_pair_form    = (0.30,  0.70),  # new: test timing–magnitude tradeoff
)

N_DRAWS = 100
N_CPUS  = 10


# --- Sampling -----------------------------------------------------------------

def sample_prior(n, seed=42):
    from scipy.stats import qmc
    keys  = list(PRIOR.keys())
    lows  = np.array([PRIOR[k][0] for k in keys])
    highs = np.array([PRIOR[k][1] for k in keys])
    sampler = qmc.LatinHypercube(d=len(PRIOR), seed=seed)
    scaled  = qmc.scale(sampler.random(n), lows, highs)
    return [dict(zip(keys, row)) for row in scaled]


# --- Sim building -------------------------------------------------------------

def make_sim_with_pars(pars, seed=1):
    sexual = sti.StructuredSexual(
        prop_f0=pars['prop_f0'],
        prop_f2=0.1,
        prop_m0=pars['prop_m0'],
        prop_m2=0.1,
        f1_conc=0.15,
        f2_conc=0.25,
        m1_conc=pars['m1_conc'],
        m2_conc=0.5,
        p_pair_form=ss.bernoulli(p=pars['p_pair_form']),  # now calibrated
        condom_data=pd.read_csv(os.path.join(DATA_DIR, 'condom_use.csv')),
        fsw_shares=ss.bernoulli(p=0.10),
        client_shares=ss.bernoulli(p=0.20),
    )

    hiv = sti.HIV(
        beta_m2f=pars['beta_m2f'],
        eff_condom=pars['eff_condom'],
        rel_init_prev=pars['rel_init_prev'],
        rel_dur_on_art=pars['rel_dur_on_art'],
        init_prev_data=pd.read_csv(os.path.join(DATA_DIR, 'init_prev_hiv.csv')),
    )

    sim = sti.Sim(
        pars=dict(use_migration=True, rand_seed=seed, n_agents=10_000,
                  start=1985, stop=2031, verbose=-1),
        datafolder=DATA_DIR + '/',
        demographics='eswatini',
        diseases=[hiv],
        networks=[sexual, ss.MaternalNet()],
        interventions=make_interventions(),
        analyzers=[hiv_epi()],
    )
    return sim


# --- Runner -------------------------------------------------------------------

def run_one(idx, pars):
    try:
        sim = make_sim_with_pars(pars, seed=idx)
        sim.run()
        df = sim.to_df(resample='year', use_years=True, sep='.')
        df['draw'] = idx
        for k, v in pars.items():
            df[f'par_{k}'] = v
        return df
    except Exception as e:
        print(f'  draw {idx} failed: {e}')
        return None


# --- Plotting -----------------------------------------------------------------

# Columns that come from PHIA household surveys (sparse — only at survey years).
# Everything else is drawn from UNAIDS modeled estimates (annual series).
PHIA_COLS = {
    'hiv.prevalence_f_15_20', 'hiv.prevalence_f_20_25',
    'hiv.prevalence_f_25_30', 'hiv.prevalence_f_30_35',
    'hiv.prevalence_m_15_20', 'hiv.prevalence_m_20_25',
    'hiv.prevalence_m_25_30', 'hiv.prevalence_m_30_35',
    'hiv_epi.prevalence_f_35_40', 'hiv_epi.prevalence_m_35_40',
    'hiv_epi.incidence_f_15_49', 'hiv_epi.incidence_m_15_49',
    'hiv_epi.incidence_f_18_49', 'hiv_epi.incidence_m_18_49',
}

# (col, title, is_percent)
TARGETS = [
    ('hiv.prevalence_15_49',      'HIV prevalence 15-49\n(UNAIDS annual + PHIA at surveys)',   True),
    ('hiv.prevalence_f_20_25',    'Female prevalence 20-25\n(PHIA surveys: 2007/2011/2016/2021)', True),
    ('hiv.prevalence_m_20_25',    'Male prevalence 20-25\n(PHIA surveys: 2007/2011/2016/2021)',   True),
    ('hiv.n_infected',            'PLHIV (total)\n(UNAIDS annual)',                              False),
    ('hiv_epi.incidence_f_15_49', 'Female incidence 15-49\n(PHIA surveys: 2011, 2016)',         True),
    ('hiv_epi.incidence_m_15_49', 'Male incidence 15-49\n(PHIA surveys: 2011, 2016)',           True),
    ('hiv.n_on_art',              'People on ART\n(programme data)',                             False),
    ('n_alive',                   'Total population\n(census/UN)',                               False),
]


def plot_coverage(all_df, calib_data, art_data, n_draws):
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    axes = axes.flatten()

    # Align prevalence_15_49 — this column appears in both PHIA (2007/2011/2016/2021)
    # and UNAIDS annual series. We show annual series as grey dots and PHIA as orange.
    # For simplicity we just show all non-NaN values; PHIA years will overlap UNAIDS.

    for ax, (col, title, is_pct) in zip(axes, TARGETS):
        if col not in all_df.columns:
            ax.text(0.5, 0.5, f'{col}\nnot in model output', ha='center', va='center',
                    transform=ax.transAxes, color='red', fontsize=9)
            ax.set_title(title, fontsize=9)
            continue

        pivot = all_df.pivot_table(index='timevec', columns='draw', values=col)
        lo  = pivot.quantile(0.10, axis=1)
        mid = pivot.quantile(0.50, axis=1)
        hi  = pivot.quantile(0.90, axis=1)

        ax.fill_between(pivot.index, lo, hi, alpha=0.20, color='steelblue')
        ax.plot(pivot.index, hi,  color='steelblue', lw=0.6, ls='--', alpha=0.5)
        ax.plot(pivot.index, lo,  color='steelblue', lw=0.6, ls='--', alpha=0.5)
        ax.plot(pivot.index, mid, color='steelblue', lw=1.5, label='Prior median')
        ax.fill_between([], [], [], alpha=0.20, color='steelblue', label='Prior 10–90%')

        # Overlay data — PHIA (orange diamonds) vs UNAIDS (grey circles)
        if col == 'hiv.n_on_art':
            # Use the separate n_art.csv (keyed 'year' and 'n_art')
            dvals = art_data.set_index('year')['n_art'].dropna()
            ax.scatter(dvals.index, dvals.values, color='dimgrey', s=18,
                       zorder=5, alpha=0.8, label='Programme data')
        elif col in calib_data.columns:
            dvals = calib_data.set_index('time')[col].dropna()
            is_phia = col in PHIA_COLS
            color  = 'darkorange'    if is_phia else 'dimgrey'
            marker = 'D'             if is_phia else 'o'
            size   = 55              if is_phia else 18
            lw     = 0.5             if is_phia else 0.0
            label  = 'PHIA survey'   if is_phia else 'UNAIDS est.'
            ax.scatter(dvals.index, dvals.values,
                       color=color, marker=marker, s=size, zorder=6,
                       edgecolors='black' if is_phia else 'none', linewidths=lw,
                       label=label)

        ax.set_title(title, fontsize=8)
        ax.set_xlabel('Year', fontsize=8)
        ax.tick_params(labelsize=7)
        if is_pct:
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(
                lambda y, _: f'{y*100:.0f}%'))
        ax.legend(fontsize=7, loc='upper left')

    fig.suptitle(
        f'Prior predictive check v3  (n={n_draws} draws) — PHIA survey data as primary targets\n'
        'Orange diamonds = PHIA survey   ·   Grey dots = UNAIDS modelled estimates\n'
        'Does the prior 10–90% band cover the PHIA survey points?',
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    return fig


# --- Diagnostic: p_pair_form vs peak year ------------------------------------

def plot_diagnostic(all_df):
    """Show how p_pair_form affects epidemic timing vs magnitude."""
    # Find peak prevalence year for each draw
    prev_col = 'hiv.prevalence_15_49'
    if prev_col not in all_df.columns:
        return None

    rows = []
    for draw_id, grp in all_df.groupby('draw'):
        grp_sorted = grp.sort_values('timevec')
        peak_idx   = grp_sorted[prev_col].idxmax()
        peak_year  = grp_sorted.loc[peak_idx, 'timevec']
        peak_prev  = grp_sorted.loc[peak_idx, prev_col]
        peak_plhiv = grp_sorted.loc[peak_idx, 'hiv.n_infected'] if 'hiv.n_infected' in grp_sorted.columns else np.nan
        rows.append({
            'draw':      draw_id,
            'peak_year': peak_year,
            'peak_prev': peak_prev,
            'peak_plhiv': peak_plhiv,
            'p_pair_form': grp_sorted['par_p_pair_form'].iloc[0] if 'par_p_pair_form' in grp_sorted.columns else np.nan,
            'beta_m2f':    grp_sorted['par_beta_m2f'].iloc[0]    if 'par_beta_m2f'    in grp_sorted.columns else np.nan,
        })
    diag = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    sc1 = axes[0].scatter(diag['beta_m2f'], diag['peak_year'],
                          c=diag['p_pair_form'], cmap='viridis', s=40, alpha=0.8)
    axes[0].axhline(2007, color='red', lw=1, ls='--', label='Target peak ~2007')
    axes[0].set_xlabel('beta_m2f')
    axes[0].set_ylabel('Peak prevalence year')
    axes[0].set_title('Timing vs beta_m2f\n(colour = p_pair_form)')
    axes[0].legend(fontsize=8)
    plt.colorbar(sc1, ax=axes[0], label='p_pair_form')

    sc2 = axes[1].scatter(diag['beta_m2f'], diag['peak_plhiv'] / 1000,
                          c=diag['p_pair_form'], cmap='viridis', s=40, alpha=0.8)
    axes[1].axhline(200, color='red', lw=1, ls='--', label='Target ~200K PLHIV')
    axes[1].set_xlabel('beta_m2f')
    axes[1].set_ylabel('Peak PLHIV (K)')
    axes[1].set_title('Magnitude vs beta_m2f\n(colour = p_pair_form)')
    axes[1].legend(fontsize=8)
    plt.colorbar(sc2, ax=axes[1], label='p_pair_form')

    sc3 = axes[2].scatter(diag['peak_year'], diag['peak_plhiv'] / 1000,
                          c=diag['p_pair_form'], cmap='viridis', s=40, alpha=0.8)
    axes[2].axvline(2007, color='red', lw=1, ls='--', label='Target year')
    axes[2].axhline(200,  color='blue', lw=1, ls='--', label='Target PLHIV')
    axes[2].set_xlabel('Peak prevalence year')
    axes[2].set_ylabel('Peak PLHIV (K)')
    axes[2].set_title('Timing vs magnitude\n(colour = p_pair_form)')
    axes[2].legend(fontsize=8)
    plt.colorbar(sc3, ax=axes[2], label='p_pair_form')

    fig.suptitle('Diagnostic: does p_pair_form resolve the timing–magnitude tradeoff?', fontsize=12)
    fig.tight_layout()
    return fig


# --- Main ---------------------------------------------------------------------

if __name__ == '__main__':
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(FIGDIR, exist_ok=True)

    draws = sample_prior(N_DRAWS)
    print(f'Running {N_DRAWS} prior draws ({N_CPUS} CPUs)...')
    t0 = sc.tic()
    results = sc.parallelize(run_one, list(enumerate(draws)), ncpus=N_CPUS, die=False)
    sc.toc(t0)

    results = [r for r in results if r is not None]
    print(f'{len(results)}/{N_DRAWS} draws completed successfully.')

    all_df = pd.concat(results, ignore_index=True)
    sc.saveobj(os.path.join(OUTDIR, 'coverage_results.obj'), all_df)

    calib_data = pd.read_csv(os.path.join(DATA_DIR, 'eswatini_hiv_calib.csv'))
    art_data   = pd.read_csv(os.path.join(DATA_DIR, 'n_art.csv'))

    fig = plot_coverage(all_df, calib_data, art_data, len(results))
    out = os.path.join(FIGDIR, '03_coverage_phia.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved: {out}')

    fig_diag = plot_diagnostic(all_df)
    if fig_diag is not None:
        out_diag = os.path.join(FIGDIR, '03_diagnostic_ppairform.png')
        fig_diag.savefig(out_diag, dpi=150, bbox_inches='tight')
        print(f'Saved: {out_diag}')
