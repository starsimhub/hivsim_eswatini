"""
Smoke test: a tiny sim should build, run, and produce a sane results DataFrame.

This test does NOT validate scientific correctness — it only catches "the sim
won't even start" regressions (broken imports, missing data files, API drift in
starsim/stisim, etc.).
"""

import pytest


@pytest.fixture(scope='module')
def tiny_sim_df():
    """ Run one very small, fast sim and return its annualized DataFrame.

    TODO: tune ``n_agents`` and ``stop`` if this runs too slowly on CI — the
    minimum that exercises the full code path is what we want.
    """
    from run_sims import make_sim
    sim = make_sim(seed=1, start=2000, stop=2005, verbose=-1)
    # Shrink the population to keep the test fast. Tweak if needed.
    sim.pars['n_agents'] = 500
    sim.run()
    return sim.to_df(resample='year', use_years=True, sep='.')


def test_smoke_runs_and_produces_results(tiny_sim_df):
    """ The sim runs without error and emits at least one year of output. """
    assert len(tiny_sim_df) > 0
    assert 'timevec' in tiny_sim_df.columns


def test_smoke_prevalence_in_unit_interval(tiny_sim_df):
    """ HIV prevalence (overall and by analyzer age band) must lie in [0, 1]. """
    prev_cols = [c for c in tiny_sim_df.columns if 'prevalence' in c]
    assert len(prev_cols) > 0, 'expected at least one prevalence column'
    for col in prev_cols:
        vals = tiny_sim_df[col].dropna()
        assert (vals >= 0).all(), f'{col} has negative values'
        assert (vals <= 1).all(), f'{col} exceeds 1.0'


def test_smoke_non_negative_counts(tiny_sim_df):
    """ Count-style results (n_infected, new_infections, n_alive) must be non-negative. """
    count_cols = [c for c in tiny_sim_df.columns
                  if any(k in c for k in ['n_infected', 'new_infections', 'n_alive'])]
    for col in count_cols:
        vals = tiny_sim_df[col].dropna()
        assert (vals >= 0).all(), f'{col} has negative values'


if __name__ == '__main__':
    pytest.main(['-x', '-v', __file__])
