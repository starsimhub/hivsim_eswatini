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

    # ANC testing: test undiagnosed pregnant women in first trimester
    def anc_eligibility(sim):
        return sim.demographics.pregnancy.tri1_uids[
            ~sim.diseases.hiv.diagnosed[sim.demographics.pregnancy.tri1_uids]
        ]

    anc_testing = sti.HIVTest(
        test_prob_data=0.9,
        dt_scale=False,
        name='anc_testing',
        eligibility=anc_eligibility,
        label='anc_testing',
    )

    tests = [fsw_testing, other_testing, low_cd4_testing, anc_testing]

    return tests


def _normalize_age_bin_format(df):
    # Convert "[15,25)" / "[15:25)" interval notation to "15-25" (the format
    # expected by upstream stisim's ss.parse_age_range as of v1.5.5).
    if 'AgeBin' in df.columns:
        df = df.copy()
        df['AgeBin'] = df['AgeBin'].astype(str).str.replace(r'^\[(\d+)[,:](\d+)\)$', r'\1-\2', regex=True)
    return df


def make_interventions(vmmc_class=None):
    # Upstream sti.VMMC gained prevalence/stock-target semantics in stisim 1.5.9
    # -- the behaviour the in-repo VMMCPrevalenceTarget subclass existed to
    # supply. Exp 017 confirmed the two are behaviourally identical (circumcision
    # 15-49: 0.084 vs 0.084 at 2005, 0.474 vs 0.475 at 2021), so vmmc.py was
    # deleted in exp 018. vmmc_class is kept as an injection point for A/B tests.
    vmmc_class = vmmc_class or sti.VMMC

    art_data = _normalize_age_bin_format(pd.read_csv('data/art_coverage.csv'))
    vmmc_data = _normalize_age_bin_format(pd.read_csv('data/vmmc_coverage.csv'))
    tests = get_testing_products()
    art = sti.ART(coverage=art_data)
    vmmc = vmmc_class(coverage=vmmc_data)

    # NO PrEP. `sti.Prep()` with coverage=None does not mean "off" -- both 1.5.8
    # and 1.5.11 fall back to a built-in ramp reaching 80% of FSW by 2025,
    # starting in 2004, a decade before PrEP had efficacy evidence. Every
    # experiment from 001 to 017 ran with that undeclared default; exp 017
    # measured realised protection at ~0.67 of uninfected FSW by 2021.
    # Removed in exp 018 (decision 2026-08-26): a fabricated intervention in the
    # calibration window biases the transmission parameters that absorb it.
    # PrEP returns, deliberately specified from programme data, for the
    # decision analysis -- it is half the question in CLAUDE.md.
    interventions = tests + [
        art,
        vmmc,
    ]

    return interventions


