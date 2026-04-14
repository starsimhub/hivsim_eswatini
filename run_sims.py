"""
Run HIV Eswatini model
"""

# %% Imports and settings
import numpy as np
import sciris as sc
import pandas as pd
import starsim as ss
import stisim as sti

# From this repo
from interventions import make_interventions
from analyzers import hiv_epi

# Constants
LOCATION = 'eswatini'
DATA_DIR = 'data'
RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'


def make_sim(seed=1, start=1985, stop=2031, verbose=1/12, analyzers=None):

    # Network
    sexual = sti.StructuredSexual(
        prop_f0=0.6,
        prop_f2=0.1,  # 60% LR, 30% MR, 10% HR
        prop_m0=0.5,
        prop_m2=0.1,  # 50% LR, 50% MR, 10% HR
        debut_pars_f=[17.5, 2.5],  # DHS eSwatini: median ~17-18yr for women
        debut_pars_m=[18.5, 2.5],  # DHS eSwatini: median ~18-19yr for men
        f1_conc=0.15,
        f2_conc=0.25,
        m1_conc=0.15,
        m2_conc=0.5,
        p_pair_form=0.5,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
        fsw_shares=ss.bernoulli(p=0.10),
        client_shares=ss.bernoulli(p=0.20),
    )
    maternal = ss.MaternalNet()
    networks = [sexual, maternal]

    # Diseases
    hiv = sti.HIV(
        beta_m2f=0.01,
        eff_condom=0.85,
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=.1,
    )

    # Interventions
    interventions = make_interventions()

    # Default analyzers
    default_analyzers = [hiv_epi()]
    if analyzers is not None:
        default_analyzers += list(analyzers)
    analyzers = default_analyzers

    simpars = dict(
        use_migration=True, rand_seed=seed, n_agents=10e3, start=start, stop=stop, verbose=verbose,
    )
    hiv_data = pd.read_csv(f'{DATA_DIR}/eswatini_hiv_calib.csv')

    sim = sti.Sim(
        pars=simpars,
        datafolder='data/',
        demographics=LOCATION.lower(),
        diseases=[hiv],
        networks=networks,
        interventions=interventions,
        analyzers=analyzers,
        data=hiv_data,
    )

    print('Created sim')
    return sim


if __name__ == '__main__':

    seed = 1

    sim = make_sim(stop=2031, seed=seed)
    sim.run()
    sim.plot('hiv', annualize=True)

    df = sim.to_df(resample='year', use_years=True, sep='.')
    sc.saveobj(f'results/{LOCATION}_sim.df', df)
