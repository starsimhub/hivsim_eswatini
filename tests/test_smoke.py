"""
Smoke test: the model builds and runs without crashing.

Runs a short window so this stays fast (~seconds) -- not a calibration
check, just confirms make_sim()/interventions/analyzers wire up correctly
against the installed stisim/starsim versions.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from run_sims import make_sim


def test_smoke():
    sim = make_sim(seed=1, start=1985, stop=1988, verbose=-1)
    sim.run()
    assert sim.results.hiv.n_infected[-1] > 0, 'HIV died out immediately -- check init_prev/beta'


if __name__ == '__main__':
    test_smoke()
    print('Smoke test passed')
