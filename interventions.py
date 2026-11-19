"""
Define Eswatini-specific interventions, including HIV testing algorithms and syphilis testing/treatment interventions.
"""

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti



def anc_test_prob(years, start=1995, p_lo=0.70, p_hi=0.90, end=2015):
    """Per-TIMESTEP ANC testing probability giving a per-PREGNANCY coverage ramp.

    Exp 027: zero before `start`, then a linear ramp from `p_lo` of pregnant
    women tested at `start` to `p_hi` at `end`, flat after. The defaults are
    the researcher-specified starting point (DTC_check/context_dtc.md: 1995,
    70% -> 90%, midpoint 2005), left over from the EMOD calibration; `start`,
    `p_lo` and `end` are calibrated in exp 027 (calib_testing.py).

    `anc_testing` runs with dt_scale=False on women in their first trimester,
    which is 13 weeks, i.e. ~3 monthly timesteps. A per-step probability p
    gives a per-pregnancy probability 1-(1-p)^3, so the coverage P is
    converted as p = 1-(1-P)^(1/3). (The pre-027 value, 0.9 per step from the
    sim start in 1985, meant ~99.9% of pregnant women were tested from 1985.)
    """
    years = np.asarray(years, dtype=float)
    P = np.interp(years, [start, max(end, start + 1e-6)], [p_lo, p_hi])
    P[years < start] = 0.0
    return 1 - (1 - P) ** (1 / 3)


# Exp 027 step 2c (DTC_check/context_dtc.md, 2026-10-28 revision): ANC testing
# "begin[s] in 1992, with midpoint in 2003 and max at 2010". Per-pregnancy
# coverage, interpolated linearly between the knots, 0 before and flat after.
# The maximum (90%) is the researcher's earlier ANC specification.
ANC_KNOTS_2C = ((1992, 0.0), (2003, 0.45), (2010, 0.90))


def anc_test_prob_knots(years, knots=ANC_KNOTS_2C):
    """Per-TIMESTEP ANC testing probability from per-PREGNANCY coverage knots.

    Same per-pregnancy -> per-step conversion as `anc_test_prob`.
    """
    k_years, k_cov = zip(*knots)
    P = np.interp(np.asarray(years, dtype=float), k_years, k_cov)
    return 1 - (1 - P) ** (1 / 3)


def anc_test_prob_logistic(years, start=1992.5, mid=2003.0, end=2010.0, p_max=0.90,
                           frac_at_end=0.95):
    """Per-TIMESTEP ANC testing probability from a logistic PER-PREGNANCY ramp.

    Exp 027 step 2d (DTC_check/context_dtc.md, 2026-11-19): "ANC testing
    starts in 1992.5, midpoint at 2003, max by 2010", and diagnosis should
    start lower than under 2c's linear knots. Coverage is 0 before `start`,
    then p_max / (1 + exp(-k (t - mid))), with k set so the curve reaches
    `frac_at_end` of p_max at `end` (95% by 2010 -> k = ln(19)/7 = 0.42/yr).
    It is then ~1% of p_max at 1992.5, rather than rising linearly from it:
    1998 coverage is ~10% here vs ~23% under ANC_KNOTS_2C.
    `p_max` (90%) is the researcher's earlier ANC level. Same per-pregnancy ->
    per-step conversion as `anc_test_prob`.
    """
    years = np.asarray(years, dtype=float)
    k = np.log(frac_at_end / (1 - frac_at_end)) / (end - mid)
    P = p_max / (1 + np.exp(-k * (years - mid)))
    P[years < start] = 0.0
    return 1 - (1 - P) ** (1 / 3)


# Exp 027 step 2d: general-population testing cohorts (DTC_check/context_dtc.md,
# 2026-11-19). Under-15s share the 0-24 rate: no 90-90-90 data exist for them,
# and the researcher's assumption is that they match 15-24 year olds.
AGE_COHORTS = ((0, 25), (25, 35), (35, 50), (50, 200))


def get_testing_products(testing_pars=None):
    """
    Define HIV products and testing interventions

    testing_pars (dict, optional): exp 027's cascade inputs. None (default)
    reproduces the testing used by experiments up to 026 exactly. Otherwise:
        rel_test_m, rel_test_f: multipliers on the general-population testing
            ramp for men and women, which is then split into two streams
            (other_testing_m / other_testing_f). Tuned so % diagnosed matches
            SHIMS2 by sex -- see experiments/027_cascade_prep_recal.
        gp_start (optional, default 1990): year general-population testing
            begins; it ramps linearly from 0 then to rel x 0.5 in 2020 and
            rel x 0.6 in 2040 (the original ramp's 2020/2040 anchors).
            gp_start_m / gp_start_f (optional) set it per sex.
        anc (bool): use an ANC ramp (anc_test_prob) instead of a flat 0.9 per
            step from the sim start. anc_start / anc_lo / anc_end (optional)
            override anc_test_prob's defaults.
        anc_knots (optional): ((year, coverage), ...) per-pregnancy ANC
            coverage knots (anc_test_prob_knots); takes precedence over
            anc_start / anc_lo / anc_end. Exp 027 step 2c uses ANC_KNOTS_2C.
        anc_logistic (optional): True, or a dict of anc_test_prob_logistic
            arguments; takes precedence over anc_knots. Exp 027 step 2d.
        rel_test_age (optional): {'m': [4 values], 'f': [4 values]}, testing
            levels per AGE_COHORTS; replaces rel_test_m/f and splits the
            general-population streams into sex x cohort (other_testing_m_0_25,
            ...). Exp 027 step 2d.
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

    if testing_pars is None:
        gp_tests = [sti.HIVTest(
            years=years,
            test_prob_data=gp_prob,
            name='other_testing',
            eligibility=other_eligibility,
            label='other_testing',
        )]
    else:
        # Exp 027: the original ramp's 2020/2040 anchors, a calibrated start
        # year (0 before it), scaled separately for men and women.
        def gp_shape(s):
            start = testing_pars.get(f'gp_start_{s}', testing_pars.get('gp_start', 1990))
            return np.interp(years, [start, 2020, 2040], [0, 0.5, 0.6])

        def male_elig(sim):
            return other_eligibility(sim) & sim.people.male

        def female_elig(sim):
            return other_eligibility(sim) & sim.people.female

        if 'rel_test_age' in testing_pars:
            # Exp 027 step 2d: one stream per sex x age cohort, same ramp shape.
            def cohort_elig(sex_elig, lo, hi):
                def elig(sim):
                    age = sim.people.age
                    return sex_elig(sim) & (age >= lo) & (age < hi)
                return elig

            gp_tests = [
                sti.HIVTest(years=years,
                            test_prob_data=np.clip(gp_shape(s) * rel, 0, 1),
                            name=f'other_testing_{s}_{lo}_{hi}',
                            eligibility=cohort_elig(elig, lo, hi),
                            label=f'other_testing_{s}_{lo}_{hi}')
                for s, elig in (('m', male_elig), ('f', female_elig))
                for (lo, hi), rel in zip(AGE_COHORTS, testing_pars['rel_test_age'][s])]
        else:
            gp_tests = [
                sti.HIVTest(years=years,
                            test_prob_data=np.clip(gp_shape(s) * testing_pars[f'rel_test_{s}'], 0, 1),
                            name=f'other_testing_{s}', eligibility=elig,
                            label=f'other_testing_{s}')
                for s, elig in (('m', male_elig), ('f', female_elig))]

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

    if testing_pars is not None and testing_pars.get('anc', False):
        anc_years = np.arange(1985, 2041)
        if testing_pars.get('anc_logistic'):
            # Quarter-year grid so the 1992.5 start is represented exactly.
            anc_years = np.arange(1985, 2041, 0.25)
            kw = testing_pars['anc_logistic']
            anc_prob = anc_test_prob_logistic(anc_years, **(kw if isinstance(kw, dict) else {}))
        elif 'anc_knots' in testing_pars:
            anc_prob = anc_test_prob_knots(anc_years, testing_pars['anc_knots'])
        else:
            anc_prob = anc_test_prob(
                anc_years, **{arg: testing_pars[k] for k, arg in
                              (('anc_start', 'start'), ('anc_lo', 'p_lo'), ('anc_end', 'end'))
                              if k in testing_pars})
        anc_testing = sti.HIVTest(
            years=anc_years,
            test_prob_data=anc_prob,
            dt_scale=False,
            name='anc_testing',
            eligibility=anc_eligibility,
            label='anc_testing',
        )
    else:
        anc_testing = sti.HIVTest(
            test_prob_data=0.9,
            dt_scale=False,
            name='anc_testing',
            eligibility=anc_eligibility,
            label='anc_testing',
        )

    tests = [fsw_testing] + gp_tests + [low_cd4_testing, anc_testing]

    return tests


def _normalize_age_bin_format(df):
    # Convert "[15,25)" / "[15:25)" interval notation to "15-25" (the format
    # expected by upstream stisim's ss.parse_age_range as of v1.5.5).
    if 'AgeBin' in df.columns:
        df = df.copy()
        df['AgeBin'] = df['AgeBin'].astype(str).str.replace(r'^\[(\d+)[,:](\d+)\)$', r'\1-\2', regex=True)
    return df


# SHIMS2 Table 15.5.A (SHIMS2 Final Report, p. 117): % of males circumcised,
# medical + non-medical, by age. Exp 027 step 2d replaces the 2016 rows of
# data/vmmc_coverage.csv with these, matching how its 2021 rows were built from
# SHIMS3 (medical + non-medical). The old 2016 rows matched SHIMS2's medical-only
# column at 15-39 but were 5-8 points low at 40-49 and 60-64 (EMOD carry-overs).
# [10,15) has no survey value and keeps its extrapolated 0.44.
SHIMS2_VMMC = {'[15,20)': 0.383, '[20,25)': 0.316, '[25,30)': 0.290, '[30,35)': 0.258,
               '[35,40)': 0.279, '[40,45)': 0.230, '[45,50)': 0.184, '[50,55)': 0.171,
               '[55,60)': 0.146, '[60,65)': 0.199}


def make_interventions(vmmc_class=None, art_vls_coverage='phia',
                       vls_stock_target=True, art_coverage=None,
                       testing_pars=None, prep=False, shims2_year=2016, pmtct=False,
                       vmmc_shims2=False):
    # art_coverage='cohorts': exp 027 step 2d's forced ART on the testing age
    # cohorts, data/art_coverage_cohorts.csv (art_coverage_cohort_construction.py).
    # vmmc_shims2: exp 027 step 2d -- SHIMS2_VMMC totals in the 2016 rows, and
    # those rows moved to shims2_year. Off by default.
    # testing_pars / prep / pmtct: exp 027's cascade, PrEP and PMTCT inputs.
    # All default to off, so every earlier experiment is unchanged. See
    # get_testing_products, prep_initiations.py and pmtct.py.
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
    if isinstance(art_coverage, str):
        if art_coverage != 'cohorts':
            raise ValueError(f"art_coverage must be None, 'cohorts' or a DataFrame; got {art_coverage!r}")
        art_coverage = pd.read_csv('data/art_coverage_cohorts.csv')
    art_data = _normalize_age_bin_format(
        pd.read_csv('data/art_coverage.csv') if art_coverage is None
        else art_coverage)
    # shims2_year: the year SHIMS2-derived inputs are pinned to. The ART and
    # VLS tables label SHIMS2 (fieldwork 2016-17) as 2016; exp 027 compares
    # the model to SHIMS2 at 2017.0, and with the 2016 label the forced ART
    # coverage has already moved 1/5 of the way to the 2021 values by then,
    # which inflated the second 90 by ~4 points. 2016 (default) leaves the
    # inputs as 024-026 used them.
    if shims2_year != 2016:
        art_data = art_data.copy()
        art_data.loc[art_data.Year == 2016, 'Year'] = shims2_year
    vmmc_raw = pd.read_csv('data/vmmc_coverage.csv')
    if vmmc_shims2:
        is16 = vmmc_raw.Year == 2016
        for ab, p in SHIMS2_VMMC.items():
            vmmc_raw.loc[is16 & (vmmc_raw.AgeBin == ab), 'p_vmmc'] = p
        vmmc_raw.loc[is16, 'Year'] = shims2_year
    vmmc_data = _normalize_age_bin_format(vmmc_raw)
    tests = get_testing_products(testing_pars)

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
    #
    # 'phia_logistic' (exp 027 step 2c): the same two survey points joined by a
    # logistic curve per sex (vls_construction.to_vls_logistic), rising from an
    # assumed 80% in the early ART era to a 98% ceiling, instead of flat to
    # SHIMS2 then linear to SHIMS3.
    if isinstance(art_vls_coverage, str):
        if art_vls_coverage not in ('phia', 'phia_logistic'):
            raise ValueError(f"art_vls_coverage must be 'phia', 'phia_logistic', None, "
                             f"or a coverage object; got {art_vls_coverage!r}")
        # Imported lazily: only needed when the default is used, and it keeps
        # interventions.py importable without the survey transcription.
        from vls_construction import build, to_vls_coverage, to_vls_logistic
        if art_vls_coverage == 'phia_logistic':
            art_vls_coverage = to_vls_logistic(build(), shims2_year=shims2_year)
        else:
            # fill_back_to: no measurement exists before SHIMS2 (2016), so
            # suppression among the treated is held flat back to the model start.
            # Early-ART-era suppression was plausibly worse, so this understates the
            # correction rather than overstating it. Recorded in 021's config.
            art_vls_coverage = to_vls_coverage(build(), fill_back_to=1985)
            if shims2_year != 2016:
                art_vls_coverage.loc[art_vls_coverage.Year == 2016, 'Year'] = shims2_year

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

    # Exp 027 step 2c: PMTCT prophylaxis for diagnosed pregnant women not on
    # ART, 2003-2010 ramp. After ART, which owns the on-ART mothers.
    if pmtct:
        from pmtct import PMTCT
        interventions = interventions + [PMTCT()]

    # Exp 027: PrEP is back, this time from programme initiation counts
    # (data/prep_initiations.csv) rather than stisim's undeclared default.
    if prep:
        from prep_initiations import make_prep
        interventions = interventions + make_prep()

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


