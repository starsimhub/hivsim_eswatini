"""
Analyzers for the Eswatini HIV model
"""

import numpy as np
import starsim as ss


class NetworkSnapshot(ss.Analyzer):
    """
    Capture a snapshot of network properties at a specified year.
    Used for the network structure figure.
    """
    def __init__(self, year=2020, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.year = year
        self.name = 'network_snapshot'
        self.risk_group_data = None
        self.debut_data = None
        self.lifetime_partners_data = None
        self.partnership_by_age = None
        self.rel_dur_data = None

    def step(self):
        if self.sim.t.yearvec[self.ti] == self.year:
            self._capture_snapshot()

    def _capture_snapshot(self):
        sim = self.sim
        nw = sim.networks.structuredsexual
        ppl = sim.people
        active = nw.participant & ppl.alive

        # Risk group composition
        rg_data = {}
        for sex_label, sex_bool in [('Female', ppl.female), ('Male', ppl.male)]:
            for rg in [0, 1, 2]:
                rg_data[(sex_label, rg)] = int(((nw.risk_group == rg) & sex_bool & active).count())
            rg_data[(sex_label, 'total')] = int((sex_bool & active).count())
        rg_data[('Female', 'fsw')] = int((nw.fsw & ppl.female & active).count())
        rg_data[('Male', 'client')] = int((nw.client & ppl.male & active).count())
        self.risk_group_data = rg_data

        # Lifetime partners for debuted agents
        debuted = nw.participant & ppl.alive & (ppl.age >= nw.debut)
        lp_data = {}
        for sex_label, sex_bool in [('Female', ppl.female), ('Male', ppl.male)]:
            mask = sex_bool & debuted
            lp_data[sex_label] = np.array(nw.lifetime_partners[mask])
        self.lifetime_partners_data = lp_data

        # Sexual debut age by sex
        debut_data = {}
        for sex_label, sex_bool in [('Female', ppl.female), ('Male', ppl.male)]:
            mask = sex_bool & debuted
            debut_data[sex_label] = np.array(nw.debut[mask])
        self.debut_data = debut_data

        # Female partnership status by age
        age_bins = np.arange(15, 51)
        pba = dict(age_bins=age_bins, prop_stable=[], prop_casual=[])
        for age in age_bins:
            in_age = ppl.female & ppl.alive & (ppl.age >= age) & (ppl.age < age + 1)
            n_total = int(in_age.count())
            if n_total > 0:
                n_stable = int((in_age & (nw.stable_partners >= 1)).count())
                n_casual = int((in_age & (nw.casual_partners >= 1)).count())
                pba['prop_stable'].append(n_stable / n_total)
                pba['prop_casual'].append(n_casual / n_total)
            else:
                pba['prop_stable'].append(np.nan)
                pba['prop_casual'].append(np.nan)
        pba['prop_stable'] = np.array(pba['prop_stable'])
        pba['prop_casual'] = np.array(pba['prop_casual'])
        self.partnership_by_age = pba

    def finalize(self):
        super().finalize()
        nw = self.sim.networks.structuredsexual
        dur_by_type = {0: [], 1: []}
        dt_year = self.sim.t.dt_year
        for _, rels in nw.relationship_durs.items():
            for rel in rels:
                etype = rel.get('edge_type', -1)
                if etype in dur_by_type:
                    dur_by_type[etype].append(rel['dur'] * dt_year)
        self.rel_dur_data = dur_by_type
