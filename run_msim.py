"""
Run multi-sim with top calibrated parameter sets.

This is the main analysis pipeline — separate from calibration.
Calibration finds the best parameters; this script runs them and
generates all the results needed for figures and tables.

Usage:
    python run_msim.py              # Run top 200 pars, generate stats
    python run_msim.py --n_pars 50  # Quick run with fewer pars
"""

import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import numpy as np
import sciris as sc
import stisim as sti
import pandas as pd
from run_sims import make_sim

LOCATION = 'eswatini'
RESULTS_DIR = 'results'

# Percentiles for summary statistics
percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]
percentiles = [.5] + [p for pair in percentile_pairs for p in pair]


def check_hiv_alive(sim):
    """Check that HIV didn't die out."""
    hiv_ni = sim.results.hiv.new_infections[-60:]  # Last 5 years
    return float(np.sum(hiv_ni)) > 0


def run_msim(n_pars=200, start=1985, stop=2026):
    """
    Run top n_pars calibrated parameter sets.
    No seed variation — each par set is a genuinely distinct fit.
    """
    pars_df = sc.loadobj(f'{RESULTS_DIR}/{LOCATION}_pars.df')
    base = make_sim(start=start, stop=stop, verbose=-1)
    msim = sti.make_calib_sims(
        calib_pars=pars_df, sim=base, n_parsets=n_pars, check_fn=check_hiv_alive,
    )
    print(f'Kept {len(msim.sims)} sims')
    return msim.sims


def prune_columns(df):
    """
    Drop columns we don't need for any figure or analysis.
    Safe to extend KEEP_PREFIXES if new analyzers are added.
    """
    KEEP_PREFIXES = [
        # Core epi
        'timevec', 'n_alive',
        'hiv.prevalence', 'hiv.n_infected', 'hiv.new_infections', 'hiv.new_deaths',
        'hiv.n_on_art', 'hiv.p_on_art', 'hiv.n_diagnosed', 'hiv.new_diagnoses',
        'hiv.new_agents_on_art', 'hiv.incidence',

        # hiv_epi analyzer (UNAIDS age ranges, PHIA fine bins, incidence)
        'hiv_epi.',
    ]

    keep = [c for c in df.columns if any(c.startswith(p) or c == p for p in KEEP_PREFIXES)]
    dropped = len(df.columns) - len(keep)
    if dropped:
        print(f'  Pruned {dropped} columns → {len(keep)} kept')
    return df[keep]


def sim_to_df(sim):
    """Convert a single sim to a pruned dataframe with par_idx."""
    df = sim.to_df(resample='year', use_years=True, sep='.')
    df['par_idx'] = sim.par_idx
    return df


def save_results(sims):
    """
    Generate percentile statistics and save.
    Parallelizes to_df() across sims, then prunes and aggregates.
    """
    print(f'Generating results from {len(sims)} sims...')

    dfs = sc.parallelize(sim_to_df, sims)
    resdf = prune_columns(pd.concat(dfs))
    print(f'  Combined DataFrame: {len(resdf)} rows, {len(resdf.columns)} columns')

    # Generate percentile statistics grouped by year
    cs = resdf.groupby(resdf.timevec).describe(percentiles=percentiles)
    sc.saveobj(f'{RESULTS_DIR}/{LOCATION}_calib_stats.df', cs)
    print(f'Saved {RESULTS_DIR}/{LOCATION}_calib_stats.df')

    return cs


if __name__ == '__main__':

    n_pars = 200
    stop = 2026

    sims = run_msim(n_pars=n_pars, stop=stop)

    if len(sims) > 0:
        cs = save_results(sims)
    else:
        print('No surviving sims — cannot generate stats')

    print('Done!')
