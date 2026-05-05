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


def make_sim(seed=1, start=1985, stop=2031, verbose=1/12, analyzers=None,
             hiv_pars=None, network_pars=None):

    # Condom data: stisim defaults x0.5 for non-marital pairings (act-level usage
    # is roughly half of DHS-reported "ever-used at last sex"). LL pairing kept
    # at marital baseline (~1%). See experiment 004 for sensitivity analysis.
    condom_df = pd.read_csv(f'data/condom_use.csv')
    year_cols = [c for c in condom_df.columns if c != 'partnership']
    non_ll = condom_df['partnership'] != '(0,0)'
    condom_df.loc[non_ll, year_cols] = condom_df.loc[non_ll, year_cols].astype(float) * 0.5

    # Network
    network_kwargs = dict(
        prop_f0=0.6,
        prop_f2=0.1,  # 60% LR, 30% MR, 10% HR
        prop_m0=0.5,
        prop_m2=0.1,  # 50% LR, 50% MR, 10% HR
        debut_pars_f=[17.5, 1],  # DHS eSwatini: median ~17-18yr for women
        debut_pars_m=[18.5, 1],  # DHS eSwatini: median ~18-19yr for men
        f1_conc=0.15,
        f2_conc=0.25,
        m1_conc=0.15,
        m2_conc=0.5,
        p_pair_form=0.5,
        condom_data=condom_df,
        fsw_shares=ss.bernoulli(p=0.10),
        client_shares=ss.bernoulli(p=0.20),
    )
    if network_pars:
        network_kwargs.update(network_pars)
    sexual = sti.StructuredSexual(**network_kwargs)
    maternal = ss.MaternalNet()
    networks = [sexual, maternal]

    # Diseases. rel_beta_f2m=0.25 (M->F per-act risk is 4x F->M) chosen in
    # experiment 004 to reproduce the PHIA F:M incidence ratio of ~2x.
    # rel_sus_age boosts susceptibility for women 15-24 — biologically motivated
    # by mucosal/cervical immaturity and consistent with PHIA's elevated young-women
    # incidence vs older women.
    hiv_kwargs = dict(
        beta_m2f=0.01,
        rel_beta_f2m=0.25,
        rel_sus_age=[(15, 25, 'f', 1.7), (25, 50, 'f', 1.0), (15, 50, 'm', 1.0)],
        eff_condom=0.85,
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=.1,
    )
    if hiv_pars:
        hiv_kwargs.update(hiv_pars)
    hiv = sti.HIV(**hiv_kwargs)

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
