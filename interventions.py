"""
Define Eswatini-specific interventions, including HIV testing algorithms and syphilis testing/treatment interventions.
"""

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti



# Adopted as model-v1.6 (exp 030). These are FITTED values, not assumptions:
# scanned against the eight SHIMS3 awareness strata, with the two informative
# male bands (25-34 and 35-49) picking out k_m = 0.5 independently of each
# other. k_f rests on the 15-24 band alone -- the other three female bands are
# saturated across the whole scanned range and carry almost no information.
#
# Pass 1.0 / 1.0 explicitly to reproduce experiments 001-029, which all ran with
# a single sex-neutral rate.
TEST_RATE_M = 0.5
TEST_RATE_F = 0.6


def get_testing_products(test_rate_m=None, test_rate_f=None,
                         test_boost=1.0, test_boost_start=2026,
                         test_boost_reach=2030):
    """
    Define HIV products and testing interventions

    Args:
        test_rate_m, test_rate_f: multipliers on the general-population testing
            ramp, by sex. Default to the fitted TEST_RATE_M / TEST_RATE_F above.

    Why the general-population ramp is split by sex (exp 030)
    ---------------------------------------------------------
    Exp 029 measured the model's awareness against SHIMS3 2021 and found it too
    high by up to 14.6 pp (men 25-34: 0.894 against 0.748), while ART coverage
    among PLHIV matched to a mean of 1.4 pp -- because ART coverage is an input
    and the linkage step absorbs the error silently.

    Most of the sex difference in *awareness* is not behavioural: men acquire
    HIV later, so their PLHIV stock at 25-34 is epidemiologically young (new
    infections are 6.45% of the stock per year against women's 1.76%) and a
    larger share has not yet had time to test. Correcting for that, the residual
    difference in mean time-unaware is roughly 30%, not the 3.6x the unaware
    fractions suggest. A modest sex split in the routine testing rate is
    therefore the right correction; a never-testing subgroup is not supported by
    the data and was considered and rejected in 029.

    Because the ramp is linear from zero, scaling it by k is the same as fitting
    its 2020 plateau: linspace(0, 0.5) * k == linspace(0, 0.5 * k). Early years
    scale too, but in absolute terms barely (1995: 0.08 -> 0.04 at k = 0.5).

    ANC testing is deliberately NOT scaled -- it is a real sex-specific route to
    diagnosis, not a modelling artefact, and women should keep it.
    """
    if test_rate_m is None:
        test_rate_m = TEST_RATE_M
    if test_rate_f is None:
        test_rate_f = TEST_RATE_F

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

    # Non-FSW agents who haven't been diagnosed or treated yet, split by sex so
    # the two rates can differ. At test_rate_m == test_rate_f == 1.0 the union of
    # these two is exactly the single `other_testing` that preceded them.
    def other_eligibility_m(sim):
        return (~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed
                & ~sim.diseases.hiv.on_art & sim.people.male)

    def other_eligibility_f(sim):
        return (~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed
                & ~sim.diseases.hiv.on_art & sim.people.female)

    # test_boost is a FORWARD-LOOKING scenario lever, and is deliberately
    # separate from test_rate_m/f.
    #
    # test_rate_m/f are the FITTED historical rates (exp 030) and must apply
    # across the whole 1990-2040 series, because that is what was fitted.
    # Multiplying them to represent a scenario silently rewrites testing back
    # to 1990 and produces a different epidemic BEFORE the scenario starts --
    # exp 034's first run did exactly that, and its cascade arms entered 2026
    # with up to 13% lower incidence than status quo purely from the
    # retroactive change. The viral-load-prevalence figure is what exposed it:
    # the curves separated before the scenario-start line.
    #
    # So the boost ramps linearly from 1.0 at test_boost_start to test_boost at
    # test_boost_reach, and is flat before and after -- matching how the ART
    # and VLS coverage targets ramp in scenarios.py.
    boost = np.ones_like(gp_prob, dtype=float)
    if test_boost != 1.0:
        span = max(test_boost_reach - test_boost_start, 1e-9)
        w = np.clip((years - test_boost_start) / span, 0.0, 1.0)
        boost = 1.0 + (test_boost - 1.0) * w

    other_testing_m = sti.HIVTest(
        years=years,
        test_prob_data=np.clip(gp_prob * test_rate_m * boost, 0, 1),
        name='other_testing_m',
        eligibility=other_eligibility_m,
        label='other_testing_m',
    )

    other_testing_f = sti.HIVTest(
        years=years,
        test_prob_data=np.clip(gp_prob * test_rate_f * boost, 0, 1),
        name='other_testing_f',
        eligibility=other_eligibility_f,
        label='other_testing_f',
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

    tests = [fsw_testing, other_testing_m, other_testing_f, low_cd4_testing,
             anc_testing]

    return tests


def _normalize_age_bin_format(df):
    # Convert "[15,25)" / "[15:25)" interval notation to "15-25" (the format
    # expected by upstream stisim's ss.parse_age_range as of v1.5.5).
    if 'AgeBin' in df.columns:
        df = df.copy()
        df['AgeBin'] = df['AgeBin'].astype(str).str.replace(r'^\[(\d+)[,:](\d+)\)$', r'\1-\2', regex=True)
    return df


def make_interventions(vmmc_class=None, art_vls_coverage='phia',
                       vls_stock_target=True, art_coverage=None,
                       test_rate_m=None, test_rate_f=None, test_boost=1.0,
                       test_boost_start=2026, test_boost_reach=2030):
    # Upstream sti.VMMC gained prevalence/stock-target semantics in stisim 1.5.9
    # -- the behaviour the in-repo VMMCPrevalenceTarget subclass existed to
    # supply. Exp 017 confirmed the two are behaviourally identical (circumcision
    # 15-49: 0.084 vs 0.084 at 2005, 0.474 vs 0.475 at 2021), so vmmc.py was
    # deleted in exp 018. vmmc_class is kept as an injection point for A/B tests.
    vmmc_class = vmmc_class or sti.VMMC

    # art_coverage: override the measured ART coverage table. Needed for the
    # decision analysis -- the cascade axis IS this table, and Eswatini's
    # headroom is in coverage (0.650-0.970 by stratum) rather than in
    # suppression among the treated, which is already 0.96+. Without an
    # override there is no way to express a treatment scale-up scenario.
    # Defaults to the measured series, so calibration behaviour is unchanged.
    art_data = _normalize_age_bin_format(
        pd.read_csv('data/art_coverage.csv') if art_coverage is None
        else art_coverage)
    vmmc_data = _normalize_age_bin_format(pd.read_csv('data/vmmc_coverage.csv'))
    tests = get_testing_products(test_rate_m=test_rate_m,
                                 test_rate_f=test_rate_f,
                                 test_boost=test_boost,
                                 test_boost_start=test_boost_start,
                                 test_boost_reach=test_boost_reach)

    # art_vls_coverage: fraction of ART initiators achieving viral suppression.
    # Defaults to 'phia' -- the measured series from vls_construction.py (SHIMS2
    # Table 9.3.A, SHIMS3 Table 8.1), adopted as model-v1.2 by exp 021.
    #
    # Passing None does NOT mean "no VLS adjustment": stisim then defaults to
    # 1.0, i.e. every treated agent is virally suppressed, transmitting at
    # effective_art_efficacy = 0.99 instead of nonsupp_art_efficacy = 0.35 -- a
    # 65x difference in residual transmission. Exp 021 measured that default as
    # overstating population viral suppression by 8.8 percentage points in 2016.
    # None is retained only so 021's control arm stays reproducible.
    if isinstance(art_vls_coverage, str):
        if art_vls_coverage != 'phia':
            raise ValueError(f"art_vls_coverage must be 'phia', None, or a "
                             f"coverage object; got {art_vls_coverage!r}")
        # Imported lazily: only needed when the default is used, and it keeps
        # interventions.py importable without the survey transcription.
        from vls_construction import build, to_vls_coverage
        # fill_back_to: no measurement exists before SHIMS2 (2016), so
        # suppression among the treated is held flat back to the model start.
        # Early-ART-era suppression was plausibly worse, so this understates the
        # correction rather than overstating it. Recorded in 021's config.
        art_vls_coverage = to_vls_coverage(build(), fill_back_to=1985)

    art = sti.ART(coverage=art_data, vls_coverage=art_vls_coverage)
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

    # Viral suppression as a stock target, adopted as model-v1.3 by exp 022.
    #
    # sti.ART applies vls_coverage only at initiation, so without this an agent
    # keeps whatever suppression status they were assigned when they started
    # treatment -- for life. Better regimens or adherence support (Eswatini's
    # TLD transition from ~2019) then cannot reach existing patients, who by
    # 2021 are most of the treated population. 021 measured the symptom:
    # realized suppression lagged its own input by 1.8 points in 2021.
    #
    # This runs after ART, so it re-targets the stock each step. Exp 022 arm C
    # measured it as tracking the input to four decimals at 2021, at no cost to
    # the prevalence fit (MAE 0.0584 -> 0.0586).
    #
    # Pass vls_stock_target=False for flow-only semantics (022 arms A and B).
    # Skipped when there is no coverage table to target, since stisim's default
    # of 1.0 leaves nothing to re-target.
    if vls_stock_target and art_vls_coverage is not None:
        from vls_stock_target import VLSStockTarget
        interventions = interventions + [
            VLSStockTarget(vls_coverage=art_vls_coverage)]

    return interventions


