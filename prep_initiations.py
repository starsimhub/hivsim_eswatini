"""PrEP driven by annual INITIATION counts, not a coverage (stock) target.

Added for exp 027. `sti.Prep` tops a pool up to a *number on PrEP* each step.
Eswatini's programme data (`data/prep_initiations.csv`: PrEPWatch for oral
TDF, MOH for lenacapavir) are the number of people *started* on PrEP each
year, which is a flow. Feeding a flow into a stock target would understate
PrEP use whenever courses last under a year and overstate it otherwise, so
this subclass starts a fixed number of people per step instead:

  - each year's count is spread evenly over that year's timesteps, with the
    fractional remainder carried forward so the annual total is exact;
  - the count is national, so it is divided by `pop_scale`, like every
    'n'-format coverage input in stisim;
  - people currently on PrEP (either product -- stisim keeps a single PrEP
    state on HIV, so the two products are mutually exclusive) and HIV-positive
    people are never initiated. Someone whose 1-year course has ended can be
    initiated again in a later year, which is how programme initiation counts
    behave (returning users are counted again);
  - within the eligible pool, people are picked by the inherited per-agent
    random `willingness`, exactly as `sti.Prep` does. So who is reached is
    random with respect to risk *within* the eligible group.

The number on PrEP at any time is an output, not an input, and can be above
or below the annual initiation count.

Only `step()` and `init_pre()` differ from `sti.Prep`; expiry, efficacy and
results (`new_prep`, `n_on_prep`) are inherited unchanged.
"""

import numpy as np
import pandas as pd
import starsim as ss
import stisim as sti
from stisim.interventions.utils import top_n_by


def prep_priority(sim):
    """Who PrEP is offered to: the groups named in DTC_check/context_dtc.md.

    Sexually debuted people aged 15+ who are any of:
      - female sex workers (current), or male clients of FSW (current);
      - people with more than one current partner;
      - "unmarried" people. The model has no marriage state, so this is
        approximated as having no current *stable* partnership in the
        structured sexual network (casual and one-time partners don't count).
    HIV-negative is applied by `PrepInitiations` itself, not here.
    """
    net = sim.networks.structuredsexual
    ppl = sim.people
    # Everyone in at least one current stable partnership.
    in_stable = np.zeros(len(ppl.uid.raw), dtype=bool)
    stable = net.edges.edge_type == net.edge_types['stable']
    in_stable[net.edges.p1[stable]] = True
    in_stable[net.edges.p2[stable]] = True
    adult_active = ppl.alive.raw & net.over_debut.raw & (ppl.age.raw >= 15)
    high_risk = net.fsw.raw | net.client.raw | (net.partners.raw > 1)
    return ss.uids(np.flatnonzero(adult_active & (high_risk | ~in_stable)))


class PrepInitiations(sti.Prep):
    """PrEP that starts a given number of people per year (see module docstring).

    Args:
        annual_n (dict): {year: number of national initiations in that year}.
            Years absent from the dict get zero initiations.
        prep_eff, prep_dur, prep_adh, eligibility, name: as `sti.Prep`.
    """

    def __init__(self, annual_n, **kwargs):
        # coverage=0 keeps sti.Prep from installing its built-in FSW ramp
        # (see the NO PrEP note in interventions.py); it is never used here.
        super().__init__(coverage=0, **kwargs)
        self.annual_n = {int(y): float(n) for y, n in annual_n.items()}
        self._carry = 0.0
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self._dt_years = float(self.t.yearvec[1] - self.t.yearvec[0])
        return

    def step(self):
        sim = self.sim
        hiv = sim.diseases.hiv

        expiring = hiv.on_prep & (hiv.prep_source == self._source_id) & (hiv.ti_prep_stop <= self.ti)
        if expiring.any():
            hiv.stop_prep(expiring.uids)

        year = int(np.floor(self.t.now('year') + 1e-9))
        n_float = self.annual_n.get(year, 0.0) / sim.pars.pop_scale * self._dt_years + self._carry
        n_want = int(n_float)
        self._carry = n_float - n_want

        n_new = 0
        if n_want > 0:
            elig = self.eligibility(sim)
            candidates = elig[~hiv.infected[elig] & ~hiv.on_prep[elig]]
            if len(candidates):
                new_uids = top_n_by(candidates, self.willingness[candidates], n_want)
                started = hiv.start_prep(new_uids, eff=self.pars.prep_eff, dur=self.pars.prep_dur,
                                         source_id=self._source_id, adh=self.pars.prep_adh)
                n_new = len(started)
            # Unmet demand (pool smaller than the count) is dropped, not carried:
            # carrying it would let a later year exceed its own programme total.

        self.results['new_prep'][self.ti] = n_new
        self.results['n_on_prep'][self.ti] = (hiv.on_prep & (hiv.prep_source == self._source_id)).count()
        return


def make_prep(path='data/prep_initiations.csv'):
    """The two Eswatini PrEP products, from programme initiation counts.

    Oral TDF: PrEPWatch annual initiations; 1-year course; 75% efficacy
      (Baeten et al., NEJM 2012, doi:10.1056/NEJMoa1108524 -- NB 75% is that
      trial's TDF-FTC arm; TDF alone was 67%).
    LA-PrEP (lenacapavir): MOH counts (personal communication); 1-year
      course; 92% efficacy.
    Both use `prep_priority` eligibility. Values after 2025 in the file are
    flat carry-forwards, relevant only to projections.
    """
    d = pd.read_csv(path, encoding='utf-8-sig')
    def counts(kind):
        s = d[d.PrEP_type == kind]
        return dict(zip(s.Year, s.Number))
    oral = PrepInitiations(counts('OralTDF'), prep_eff=0.75, prep_dur=ss.years(1),
                           prep_adh=1.0, eligibility=prep_priority, name='prep_oral')
    la = PrepInitiations(counts('LA-PREP_LEN'), prep_eff=0.92, prep_dur=ss.years(1),
                         prep_adh=1.0, eligibility=prep_priority, name='prep_la')
    return [oral, la]
