"""PMTCT prophylaxis for diagnosed pregnant women who are not on ART.

Added for exp 027 step 2c (DTC_check/context_dtc.md, 2026-10-28 revision:
"begin PMTCT in 2003 with max at 2010").

The gap this closes
-------------------
In stisim, mother-to-child transmission is only reduced when the mother is on
ART (`sti.ART.step`, `pmtct_efficacy`). In this model ART is a forced stock --
coverage among all PLHIV by age and sex, from `data/art_coverage.csv` -- and
pregnancy gives no priority for it, so a woman diagnosed at ANC gets nothing for
her baby unless she happens to hold one of the forced ART places.

Eswatini's PMTCT programme began in 2003 with short-course ARV prophylaxis
(not lifelong ART) and reached most diagnosed pregnant women by about 2010.
This module represents that prophylaxis directly:

  - once per pregnancy, a diagnosed HIV-positive pregnant woman who is not on
    ART is enrolled with probability `coverage(t)`, a linear ramp from 0 in
    2003 to `p_max` in 2010, flat after;
  - while she is pregnant and enrolled (and still not on ART), her unborn
    child's susceptibility via the maternal network is multiplied by
    1 - `efficacy`, the same mechanism `sti.ART` uses for mothers on ART;
  - it does not touch the mother's own state, so the forced ART stock and all
    adult transmission are unchanged.

Women diagnosed later in a pregnancy get their draw when they are diagnosed.
Mothers on ART are left to `sti.ART`, so protection is never applied twice.
Must run after `sti.ART` in the intervention list, like the ART PMTCT step,
because `HIV.update_transmission` resets `rel_sus` at the start of each step.

`efficacy` defaults to stisim's ART `pmtct_efficacy` (0.96). Pre-2010
short-course regimens were less effective than that; this only affects
paediatric infections, and children infected in 2003 or later turn 15 in 2018
or later, so no fitted target (ages 15+, 2007-2017) is affected.
"""

import numpy as np
import starsim as ss


class PMTCT(ss.Intervention):
    """PMTCT prophylaxis for diagnosed pregnant women not on ART (see module docstring).

    Args:
        years (list): knot years of the coverage ramp (default 2003, 2010)
        coverage (list): share of eligible pregnancies enrolled at each knot
            (default 0, 0.9); interpolated linearly, flat outside the knots
        efficacy (float): reduction in the unborn child's susceptibility
    """

    def __init__(self, years=(2003, 2010), coverage=(0.0, 0.9), efficacy=0.96, **kwargs):
        super().__init__(name=kwargs.pop('name', 'pmtct'), **kwargs)
        self.define_pars(
            efficacy=efficacy,
            p_enroll=ss.bernoulli(p=0),  # p is set each step from the ramp
        )
        self.years = np.asarray(years, dtype=float)
        self.coverage = np.asarray(coverage, dtype=float)
        self.define_states(
            ss.BoolState('enrolled', label='Enrolled in PMTCT this pregnancy'),
            ss.FloatArr('ti_drawn', label='ti_pregnant of the pregnancy last drawn for'),
        )
        return

    def step(self):
        sim = self.sim
        hiv = sim.diseases.hiv
        preg = sim.demographics.pregnancy

        # Enrolment ends with the pregnancy
        self.enrolled[~preg.pregnant] = False

        # One draw per pregnancy, once the woman is diagnosed
        p = float(np.interp(self.t.now('year'), self.years, self.coverage))
        eligible = preg.pregnant & hiv.infected & hiv.diagnosed & ~hiv.on_art
        undrawn = eligible.uids[self.ti_drawn[eligible.uids] != preg.ti_pregnant[eligible.uids]]
        if len(undrawn):
            self.ti_drawn[undrawn] = preg.ti_pregnant[undrawn]
            if p > 0:
                self.pars.p_enroll.set(p=p)
                self.enrolled[self.pars.p_enroll.filter(undrawn)] = True

        # Protect the unborn children of enrolled mothers not on ART
        protected = (self.enrolled & preg.pregnant & ~hiv.on_art).uids
        if len(protected) and hasattr(sim.networks, 'maternalnet'):
            unborn = sim.networks.maternalnet.find_contacts(protected)
            hiv.rel_sus[ss.uids(unborn)] *= 1 - self.pars.efficacy
        return
