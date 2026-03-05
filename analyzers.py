"""
Analyzers for the Eswatini HIV model
"""

import numpy as np
import starsim as ss


def count(arr): return np.count_nonzero(arr)


class hiv_epi(ss.Analyzer):
    """
    Track HIV indicators by custom age ranges to match UNAIDS and PHIA calibration targets.

    Produces results that map directly to columns in eswatini_hiv_calib.csv:
        - Prevalence, incidence, n_infected, new_infections for broad UNAIDS age ranges
          (0-14, 10-19, 15-24, 15-49, 15-100)
        - Prevalence by 5-year age/sex bins above 35 to match PHIA surveys
          (35-40, 40-45, 45-50, 50-55, 55-60, 60-65)
        - Sex-specific incidence for 15-49 and 18-49
    """

    # Broad age ranges from UNAIDS data
    unaids_bins = [(0, 14), (10, 19), (15, 24), (15, 49), (15, 100)]

    # 5-year PHIA bins above 35 (below 35 the sim already has these natively)
    phia_bins = [(35, 40), (40, 45), (45, 50), (50, 55), (55, 60), (60, 65)]

    # Sex-specific incidence ranges from PHIA surveys
    incidence_bins = [(15, 49), (18, 49)]

    def init_results(self):
        super().init_results()
        results = []

        # UNAIDS broad age ranges: prevalence, n_infected, new_infections
        for lo, hi in self.unaids_bins:
            results.append(ss.Result(f'prevalence_{lo}_{hi}', dtype=float, scale=False))
            results.append(ss.Result(f'n_infected_{lo}_{hi}', dtype=int, scale=True))
            results.append(ss.Result(f'new_infections_{lo}_{hi}', dtype=int, scale=True))

        # PHIA 5-year age/sex bins above 35: prevalence only
        for sex in ['f', 'm']:
            for lo, hi in self.phia_bins:
                results.append(ss.Result(f'prevalence_{sex}_{lo}_{hi}', dtype=float, scale=False))

        # Sex-specific incidence for PHIA ranges
        for sex in ['f', 'm']:
            for lo, hi in self.incidence_bins:
                results.append(ss.Result(f'incidence_{sex}_{lo}_{hi}', dtype=float, scale=False))

        # Overall incidence for broad range
        results.append(ss.Result('incidence_15_100', dtype=float, scale=False))

        self.define_results(*results)

    def step(self):
        sim = self.sim
        ti = self.ti
        hiv = sim.diseases.hiv
        ppl = sim.people
        alive = ppl.alive
        dt = sim.t.dt_year

        # UNAIDS broad age ranges (both sexes)
        for lo, hi in self.unaids_bins:
            in_bin = alive & (ppl.age >= lo) & (ppl.age < hi)
            n_bin = in_bin.count()
            if n_bin > 0:
                self.results[f'prevalence_{lo}_{hi}'][ti] = float(np.mean(hiv.infected[in_bin]))
                self.results[f'n_infected_{lo}_{hi}'][ti] = count(hiv.infected[in_bin])
                self.results[f'new_infections_{lo}_{hi}'][ti] = count(hiv.ti_infected[in_bin] == ti)

        # PHIA 5-year bins above 35 by sex
        for sex_key, sex_bool in [('f', ppl.female), ('m', ppl.male)]:
            for lo, hi in self.phia_bins:
                in_bin = alive & sex_bool & (ppl.age >= lo) & (ppl.age < hi)
                n_bin = in_bin.count()
                if n_bin > 0:
                    self.results[f'prevalence_{sex_key}_{lo}_{hi}'][ti] = float(np.mean(hiv.infected[in_bin]))

        # Sex-specific incidence for PHIA ranges
        for sex_key, sex_bool in [('f', ppl.female), ('m', ppl.male)]:
            for lo, hi in self.incidence_bins:
                in_bin = alive & sex_bool & (ppl.age >= lo) & (ppl.age < hi)
                n_susceptible = count(~hiv.infected[in_bin])
                if n_susceptible > 0:
                    new_inf = count(hiv.ti_infected[in_bin] == ti)
                    self.results[f'incidence_{sex_key}_{lo}_{hi}'][ti] = new_inf / (n_susceptible * dt)

        # Overall incidence 15-100
        in_bin = alive & (ppl.age >= 15) & (ppl.age < 100)
        n_susceptible = count(~hiv.infected[in_bin])
        if n_susceptible > 0:
            new_inf = count(hiv.ti_infected[in_bin] == ti)
            self.results['incidence_15_100'][ti] = new_inf / (n_susceptible * dt)


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
