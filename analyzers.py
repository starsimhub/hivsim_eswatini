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


class PopByAgeSex(ss.Analyzer):
    """Population and infection counts by 5-year age band and sex, 0-100.

    Complements `hiv_epi`, which carries the calibration-target age ranges but
    not the alive-counts needed as stratified denominators. Promoted here by exp
    018 after being written locally three times (016, 017, 018): 016 needed it
    for the mortality attribution, 018 needs it to convert population-scaled
    results back into *agent* counts for the rare-event floor check.
    """
    AGE_BINS = [(a, a + 5) for a in range(0, 100, 5)]

    def __init__(self, *args, name='popagesex', **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name

    def init_results(self):
        super().init_results()
        res = []
        for sex in ('f', 'm'):
            for lo, hi in self.AGE_BINS:
                res.append(ss.Result(f'n_alive_{sex}_{lo}_{hi}', dtype=int, scale=True))
                res.append(ss.Result(f'n_infected_{sex}_{lo}_{hi}', dtype=int, scale=True))
                # Added by exp 023: incidence by age band needs a numerator.
                # hiv_epi carries incidence only as 15-49 and 18-49 aggregates
                # by sex, but the SHIMS2 incidence targets are 5- and 10-year
                # bands. Denominator is n_alive - n_infected (HIV has no
                # recovered class, so that is the susceptible stock).
                res.append(ss.Result(f'new_infections_{sex}_{lo}_{hi}',
                                     dtype=int, scale=True))
            res.append(ss.Result(f'prevalence_{sex}_15_49', dtype=float, scale=False))
        res.append(ss.Result('n_alive_total', dtype=int, scale=True))
        res.append(ss.Result('n_infected_total', dtype=int, scale=True))
        res.append(ss.Result('new_infections_total', dtype=int, scale=True))
        self.define_results(*res)

    def step(self):
        sim, ti = self.sim, self.ti
        ppl, hiv = sim.people, sim.diseases.hiv
        alive = ppl.alive

        for sex, sex_bool in (('f', ppl.female), ('m', ppl.male)):
            for lo, hi in self.AGE_BINS:
                in_bin = alive & sex_bool & (ppl.age >= lo) & (ppl.age < hi)
                self.results[f'n_alive_{sex}_{lo}_{hi}'][ti] = in_bin.count()
                self.results[f'n_infected_{sex}_{lo}_{hi}'][ti] = (in_bin & hiv.infected).count()
                self.results[f'new_infections_{sex}_{lo}_{hi}'][ti] = (
                    in_bin & (hiv.ti_infected == ti)).count()
            adults = alive & sex_bool & (ppl.age >= 15) & (ppl.age < 50)
            if adults.count() > 0:
                self.results[f'prevalence_{sex}_15_49'][ti] = float(np.mean(hiv.infected[adults]))

        self.results['n_alive_total'][ti] = alive.count()
        self.results['n_infected_total'][ti] = (alive & hiv.infected).count()
        self.results['new_infections_total'][ti] = (alive & (hiv.ti_infected == ti)).count()


class Cascade(ss.Analyzer):
    """The 95-95-95 cascade by sex: diagnosed, on ART, virally suppressed.

    Exists because `hiv.p_on_art` pools the sexes and stisim carries no result
    at all for the `on_effective_art` / `on_nonsuppressive_art` split, so the
    model's *population* viral suppression -- the headline PHIA indicator -- is
    not observable from the standard outputs.

    That split matters more than it looks: suppressed agents transmit at
    `effective_art_efficacy` = 0.99 and unsuppressed ones at
    `nonsupp_art_efficacy` = 0.35, a 65x difference in residual transmission.

    Two age ranges: 15-49 to match PHIA's headline tables, and 15+ for the
    all-adult figures those tables also report.

    `p_vls` is suppression among all PLHIV (the PHIA indicator, an *outcome* of
    coverage x suppression). `p_vls_given_art` is suppression among the treated
    (the `sti.ART(vls_coverage=...)` *input*). Conflating them double-counts the
    coverage ramp -- see vls_construction.py.

    Added by exp 021.
    """

    def __init__(self, *args, name='cascade', **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name

    RANGES = {'15_49': (15, 50), '15plus': (15, 200)}

    def init_results(self):
        super().init_results()
        res = []
        for rng in self.RANGES:
            for sex in ('f', 'm', 'all'):
                for stem, dtype, scale in (
                        ('n_infected', int, True), ('n_on_art', int, True),
                        ('n_effective_art', int, True),
                        ('n_nonsupp_art', int, True),
                        ('p_on_art', float, False), ('p_vls', float, False),
                        ('p_vls_given_art', float, False)):
                    res.append(ss.Result(f'{stem}_{sex}_{rng}', dtype=dtype,
                                         scale=scale))
        self.define_results(*res)
        return

    def step(self):
        sim, ti = self.sim, self.ti
        ppl, hiv = sim.people, sim.diseases.hiv
        for rng, (lo, hi) in self.RANGES.items():
            in_age = ppl.alive & (ppl.age >= lo) & (ppl.age < hi)
            for sex, sex_bool in (('f', ppl.female), ('m', ppl.male),
                                  ('all', np.ones(len(ppl), dtype=bool))):
                base = in_age & sex_bool
                inf = (base & hiv.infected).count()
                art = (base & hiv.on_art).count()
                eff = (base & hiv.on_effective_art).count()
                nsp = (base & hiv.on_nonsuppressive_art).count()
                r = self.results
                r[f'n_infected_{sex}_{rng}'][ti] = inf
                r[f'n_on_art_{sex}_{rng}'][ti] = art
                r[f'n_effective_art_{sex}_{rng}'][ti] = eff
                r[f'n_nonsupp_art_{sex}_{rng}'][ti] = nsp
                r[f'p_on_art_{sex}_{rng}'][ti] = art / inf if inf else np.nan
                r[f'p_vls_{sex}_{rng}'][ti] = eff / inf if inf else np.nan
                r[f'p_vls_given_art_{sex}_{rng}'][ti] = eff / art if art else np.nan
        return
