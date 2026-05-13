"""
Run HIV Eswatini model
"""

# %% Imports and settings
from pathlib import Path
import numpy as np
import sciris as sc
import pandas as pd
import starsim as ss
import stisim as sti

# From this repo
from interventions import make_interventions
from analyzers import hiv_epi

# Constants — paths anchored to the repo root so scripts work regardless of CWD.
LOCATION = 'eswatini'
REPO_DIR = Path(__file__).parent
DATA_DIR = REPO_DIR / 'data'
RESULTS_DIR = REPO_DIR / 'results'
FIGURES_DIR = REPO_DIR / 'figures'


def make_sim(seed=1, start=1985, stop=2031, verbose=1/12, analyzers=None):
    """ Build the central Eswatini HIV simulation.

    Constructs the structured sexual network, the HIV disease module, default
    interventions (testing, ART, PrEP), and the ``hiv_epi`` analyzer. Additional
    analyzers can be appended via ``analyzers``. Many model assumptions (network
    proportions, beta, condom efficacy) are intentionally hardcoded here and are
    only swapped out by the calibrator.

    Args:
        seed (int):  Random seed (``rand_seed``) for the sim.
        start (int|float):  Sim start year.
        stop (int|float):  Sim stop year.
        verbose (float|int):  Verbosity passed to ``sti.Sim``; ``-1`` suppresses output.
        analyzers (list|None):  Extra analyzers to append to the default ``[hiv_epi()]``.

    Returns:
        sti.Sim: An uninitialized sim ready to ``.run()``.

    **Example**:

        sim = make_sim(seed=1, stop=2026)
        sim.run()
    """

    # Network
    sexual = sti.StructuredSexual(
        prop_f0=0.6,
        prop_f2=0.1,  # 60% LR, 30% MR, 10% HR
        prop_m0=0.5,
        prop_m2=0.1,  # 50% LR, 50% MR, 10% HR
        f1_conc=0.15,
        f2_conc=0.25,
        m1_conc=0.15,
        m2_conc=0.5,
        p_pair_form=0.5,
        condom_data=pd.read_csv(DATA_DIR / 'condom_use.csv'),
        fsw_shares=ss.bernoulli(p=0.10),
        client_shares=ss.bernoulli(p=0.20),
    )
    maternal = ss.MaternalNet()
    networks = [sexual, maternal]

    # Diseases
    hiv = sti.HIV(
        beta_m2f=0.01,
        eff_condom=0.85,
        init_prev_data=pd.read_csv(DATA_DIR / 'init_prev_hiv.csv'),
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
        use_migration=True, rand_seed=seed, n_agents=int(10e3), start=start, stop=stop, verbose=verbose,
    )
    hiv_data = pd.read_csv(DATA_DIR / 'eswatini_hiv_calib.csv')

    sim = sti.Sim(
        pars=simpars,
        datafolder=str(DATA_DIR) + '/',
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
    sc.saveobj(RESULTS_DIR / f'{LOCATION}_sim.df', df)
