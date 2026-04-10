"""
Define Eswatini-specific interventions, including HIV testing algorithms and syphilis testing/treatment interventions.
"""

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti



def get_testing_products():
    """
    Define HIV products and testing interventions
    """

    scaleup_years = np.arange(1990, 2021)  # Years for testing
    years = np.arange(1990, 2041)  # Years for simulation
    n_years = len(scaleup_years)
    fsw_prob = np.concatenate([np.linspace(0, 0.75, n_years), np.linspace(0.75, 0.85, len(years) - n_years)])
    low_cd4_prob = np.concatenate([np.linspace(0, 0.85, n_years), np.linspace(0.85, 0.95, len(years) - n_years)])
    gp_prob = np.concatenate([np.linspace(0, 0.5, n_years), np.linspace(0.5, 0.6, len(years) - n_years)])

    # FSW agents who haven't been diagnosed or treated yet
    def fsw_eligibility(sim):
        return sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    fsw_testing = sti.HIVTest(
        years=years,
        test_prob_data=fsw_prob,
        name='fsw_testing',
        eligibility=fsw_eligibility,
        label='fsw_testing',
    )

    # Non-FSW agents who haven't been diagnosed or treated yet
    def other_eligibility(sim):
        return ~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    other_testing = sti.HIVTest(
        years=years,
        test_prob_data=gp_prob,
        name='other_testing',
        eligibility=other_eligibility,
        label='other_testing',
    )

    # Agents whose CD4 count is below 200.
    def low_cd4_eligibility(sim):
        return (sim.diseases.hiv.cd4 < 200) & ~sim.diseases.hiv.diagnosed

    low_cd4_testing = sti.HIVTest(
        years=years,
        test_prob_data=low_cd4_prob,
        name='low_cd4_testing',
        eligibility=low_cd4_eligibility,
        label='low_cd4_testing',
    )
    tests = [fsw_testing, other_testing, low_cd4_testing]

    return tests


def make_interventions():

    n_art = pd.read_csv(f'data/n_art.csv').set_index('year')
    n_art['p_art'] = np.nan
    n_art.loc[2025:, 'p_art'] = 0.90  # Switch to proportion target after historical data ends
    # n_vmmc = pd.read_csv(f'data/n_vmmc.csv').set_index('year')
    tests = get_testing_products()
    art = sti.ART(coverage=n_art)
    # vmmc = sti.VMMC(coverage=n_vmmc)
    prep = sti.Prep()

    interventions = tests + [
        art,
        # vmmc,
        prep,
    ]

    return interventions


