"""
Run sims and extract network structure data for the network figure.
"""

import sciris as sc
import numpy as np
import starsim as ss
import stisim as sti
from run_sims import make_sim
from analyzers import NetworkSnapshot

RESULTS_DIR = 'results'


def run_network_sims(n_sims=3):
    """Run sims with network analyzers."""
    sims = sc.autolist()
    for seed in range(1, n_sims + 1):
        network_analyzers = [
            sti.partner_age_diff(year=2020),
            NetworkSnapshot(year=2020),
        ]
        sim = make_sim(seed=seed, stop=2026, analyzers=network_analyzers)
        sims += sim

    sims = ss.parallel(sims).sims
    return sims


def extract_network_data(sims):
    """Extract network data from completed sims and save."""
    data = sc.objdict()

    # Lifetime partner distributions (aggregate across sims)
    lp_f, lp_m = [], []
    for sim in sims:
        snap = sim.analyzers['network_snapshot']
        lp_f.extend(list(snap.lifetime_partners_data['Female']))
        lp_m.extend(list(snap.lifetime_partners_data['Male']))
    data.lifetime_partners_f = np.array(lp_f)
    data.lifetime_partners_m = np.array(lp_m)

    # Age differences (aggregate across sims)
    age_diffs = {}
    for sim in sims:
        pad = sim.analyzers['partner_age_diff']
        for key, vals in pad.age_diffs.items():
            age_diffs.setdefault(key, []).extend(list(vals))
    data.age_diffs = {k: np.array(v) for k, v in age_diffs.items()}

    # Risk groups, debut, partnership-by-age from first sim
    snap = sims[0].analyzers['network_snapshot']
    data.risk_group_data = snap.risk_group_data
    data.debut_data = snap.debut_data
    data.rel_dur_data = snap.rel_dur_data

    # Average partnership-by-age across sims
    all_stable, all_casual = [], []
    for sim in sims:
        s = sim.analyzers['network_snapshot']
        all_stable.append(s.partnership_by_age['prop_stable'])
        all_casual.append(s.partnership_by_age['prop_casual'])
    data.partnership_by_age = dict(
        age_bins=snap.partnership_by_age['age_bins'],
        prop_stable=np.nanmean(all_stable, axis=0),
        prop_casual=np.nanmean(all_casual, axis=0),
    )

    sc.saveobj(f'{RESULTS_DIR}/network_data.obj', data)
    print(f'Saved network data to {RESULTS_DIR}/network_data.obj')
    return data


if __name__ == '__main__':
    sims = run_network_sims(n_sims=3)
    extract_network_data(sims)
