"""Viral suppression as a stock target, so it can improve for people already on ART.

The gap this closes
-------------------
`sti.ART(vls_coverage=...)` accepts a time-varying series, but it is applied
**only at ART initiation**: `HIV.on_effective_art` is written in exactly three
places -- `start_art()`, and cleared on death and on ART discontinuation. An
agent who started treatment in 2010 therefore keeps their 2010 suppression
status for life.

That makes the real-world mechanism unrepresentable. Eswatini's rise in the
third 90 came substantially from **better regimens and better adherence
support** -- the dolutegravir (TLD) transition from around 2019, and adherence
programmes -- which improve suppression among people *already on treatment*,
not only among new starters. Under flow-only semantics the existing treated
stock never benefits, and by 2021 the stock is most of the treated population.

Exp 021 measured the symptom: realized suppression among the treated lagged its
own input (0.956 against 0.967 in 2021) because the stock carries earlier
initiation years' values.

This is structurally the same defect exp 015 fixed for VMMC -- coverage applied
as a per-step hazard on new entrants rather than as a prevalence (stock) target
-- and it takes the same shape of fix.

How it works
------------
Each step, for each sex, the agents currently on ART are ranked by a persistent
per-agent `suppression_propensity` (a random score drawn once, mirroring
`sti.VMMC`'s `willingness`), and the top `target x n` are marked suppressed.

Ranking by a *persistent* score rather than re-drawing each step matters twice
over. It keeps agents from churning randomly between suppressed and
unsuppressed, and it means a rising target promotes people in a stable order --
so improving suppression adds to the suppressed pool rather than reshuffling
it. It also reads as a plausible individual trait: adherence propensity is
persistent, not redrawn monthly.

Because `HIV.update_transmission()` recomputes `rel_trans` from
`on_effective_art` each step, and `get_art_mortality_hazard()` reads the same
state, flipping these states is sufficient -- both transmission and on-ART
mortality follow automatically. Interventions run before the disease's state
update within a step, so a change lands the same step it is made.

Usage
-----
    from vls_stock_target import VLSStockTarget
    from vls_construction import build, to_vls_coverage
    tbl = to_vls_coverage(build(), fill_back_to=1985)
    sim = make_sim(..., art_vls_coverage=tbl,
                   extra_interventions=[VLSStockTarget(vls_coverage=tbl)])

Passing the same table to both is deliberate: `sti.ART` sets the status at
initiation and this re-targets the stock, so the two agree at every year and
the difference is only whether existing patients are updated. Used by exp 022.
"""

import numpy as np
import pandas as pd
import starsim as ss


class VLSStockTarget(ss.Intervention):
    """Re-target viral suppression among the treated stock each step.

    Args:
        vls_coverage: DataFrame with Year / Gender / p_vls, as produced by
            `vls_construction.to_vls_coverage`. Gender is 'm'/'f'. Values
            outside the year range are held flat at the nearest endpoint.
    """

    def __init__(self, vls_coverage=None, name='vls_stock_target', **kwargs):
        super().__init__(**kwargs)
        self.name = name
        if vls_coverage is None:
            raise ValueError('VLSStockTarget requires vls_coverage; passing None '
                             'would silently leave suppression at whatever '
                             'sti.ART assigned at initiation')
        self._raw = pd.DataFrame(vls_coverage)
        # Persistent per-agent adherence propensity, mirroring sti.VMMC's
        # willingness. Drawn once per agent, never redrawn.
        self.suppression_propensity = ss.FloatArr('suppression_propensity',
                                                  default=ss.random())
        self.targets = {}

    def init_pre(self, sim):
        super().init_pre(sim)
        need = {'Year', 'Gender', 'p_vls'}
        missing = need - set(self._raw.columns)
        if missing:
            raise ValueError(f'vls_coverage missing columns: {sorted(missing)}')
        for sex, g in self._raw.groupby('Gender'):
            g = g.sort_values('Year')
            self.targets[str(sex)] = (g.Year.to_numpy(dtype=float),
                                      g.p_vls.to_numpy(dtype=float))
        if set(self.targets) != {'m', 'f'}:
            raise ValueError(f"vls_coverage must cover both sexes, got "
                             f"{sorted(self.targets)}; any sex left out would "
                             f"keep stisim's 100%-suppressed default")
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('n_promoted', dtype=int, label='Newly suppressed by retargeting'),
            ss.Result('n_demoted', dtype=int, label='Newly unsuppressed by retargeting'),
        )
        return

    def target_for(self, sex, year):
        """Suppression target for one sex at a given year, held flat outside range."""
        years, vals = self.targets[sex]
        return float(np.interp(year, years, vals))   # np.interp clamps at both ends

    def step(self):
        sim, ti = self.sim, self.ti
        hiv, ppl = sim.diseases.hiv, sim.people
        year = float(sim.t.yearvec[ti])
        promoted = demoted = 0

        for sex, is_sex in (('f', ppl.female), ('m', ppl.male)):
            uids = (hiv.on_art & is_sex & ppl.alive).uids
            if not len(uids):
                continue
            target = self.target_for(sex, year)
            n_supp = int(round(target * len(uids)))

            # Highest propensity suppressed first. Stable across steps, so a
            # rising target adds to the suppressed pool rather than reshuffling.
            order = np.argsort(-self.suppression_propensity[uids])
            supp, nonsupp = uids[order[:n_supp]], uids[order[n_supp:]]

            promoted += int((~hiv.on_effective_art[supp]).sum()) if len(supp) else 0
            demoted += int(hiv.on_effective_art[nonsupp].sum()) if len(nonsupp) else 0

            if len(supp):
                hiv.on_effective_art[supp] = True
                hiv.on_nonsuppressive_art[supp] = False
            if len(nonsupp):
                hiv.on_effective_art[nonsupp] = False
                hiv.on_nonsuppressive_art[nonsupp] = True

        self.results['n_promoted'][ti] = promoted
        self.results['n_demoted'][ti] = demoted
        return
