"""
Experiment 02: Prior predictive check v2 — corrected prior and cause-deleted deaths.

Changes from experiment 01:
  - beta_m2f upper bound widened 0.014 → 0.020 (best draws in 01 were 25% below data PLHIV)
  - rel_init_prev narrowed 0.05–5.0 → 0.05–0.50 (upper bound of 5.0 caused epidemic 15yr too early)
  - eswatini_deaths.csv now contains cause-deleted rates (AIDS-era spike removed)

Run from the repo root:
    uv run python experiments/02_coverage_check_v2/run.py
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
    beta_m2f       = (0.002, 0.020),  # widened from 0.014
    eff_condom     = (0.50,  0.90),
    rel_init_prev  = (0.05,  0.50),   # narrowed from 5.0
    rel_dur_on_art = (1.0,   20.0),
    prop_f0        = (0.55,  0.90),
    prop_m0        = (0.55,  0.80),
    m1_conc        = (0.05,  0.30),
)

N_DRAWS = 100  # more draws than 01 for a sharper picture
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
        p_pair_form=0.5,
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

TARGETS = [
    ('hiv.prevalence_15_49', 'HIV prevalence 15–49', True),
    ('hiv.new_infections',   'New HIV infections',   False),
    ('hiv.n_on_art',         'People on ART',        False),
    ('hiv.n_infected',       'PLHIV',                False),
    ('hiv.new_deaths',       'HIV deaths',           False),
    ('n_alive',              'Total population',     False),
]


def plot_coverage(all_df, data, n_draws):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    for ax, (col, title, is_pct) in zip(axes, TARGETS):
        if col not in all_df.columns:
            ax.text(0.5, 0.5, f'{col}\nnot in results', ha='center', va='center',
                    transform=ax.transAxes, color='red')
            ax.set_title(title)
            continue

        pivot = all_df.pivot_table(index='timevec', columns='draw', values=col)
        lo  = pivot.quantile(0.10, axis=1)
        mid = pivot.quantile(0.50, axis=1)
        hi  = pivot.quantile(0.90, axis=1)

        ax.fill_between(pivot.index, lo, hi, alpha=0.25, color='steelblue')
        ax.plot(pivot.index, hi,  color='steelblue', lw=0.8, ls='--', alpha=0.6)
        ax.plot(pivot.index, lo,  color='steelblue', lw=0.8, ls='--', alpha=0.6)
        ax.plot(pivot.index, mid, color='steelblue', lw=1.5, label='Prior median')
        ax.fill_between([], [], [], alpha=0.25, color='steelblue', label='Prior 10–90%')

        if col in data.columns:
            dvals = data.set_index('time')[col].dropna()
            ax.scatter(dvals.index, dvals.values, color='black', s=20, zorder=5, label='Data')

        ax.set_title(title, fontsize=11)
        ax.set_xlabel('Year')
        if is_pct:
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(
                lambda y, _: f'{y*100:.0f}%'))
        ax.legend(fontsize=8)

    fig.suptitle(
        f'Prior predictive check v2  (n={n_draws} draws) — corrected prior + cause-deleted deaths\n'
        'Does the data fall inside the prior envelope?',
        fontsize=12, y=1.01,
    )
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

    data = pd.read_csv(os.path.join(DATA_DIR, 'eswatini_hiv_calib.csv'))
    fig  = plot_coverage(all_df, data, len(results))
    out  = os.path.join(FIGDIR, '02_coverage_check_v2.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved: {out}')
