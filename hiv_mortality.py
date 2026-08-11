"""HIV mortality multiplier (calibration knob, not a model fix).

Upstream `sti.HIV.make_p_hiv_death()` hard-codes the CD4-stratified annual
death rates `[0.003, 0.003, 0.005, 0.01, 0.05, 0.30]` inside the method body,
so there is no parameter to calibrate. Experiment 009's coverage check failed
with too little mortality flow (AIDS deaths below UNAIDS at the 1995-2005 peak
and near-zero in the post-2010 plateau, with prevalence drifting up as PLHIV
accumulate), and named this rate block as one of three fixed-parameter suspects.
This subclass exposes a single scalar `mort_mult` over those rates so the
coverage check can ask whether the data is reachable.

Note this is a *calibration parameter*, not a bug fix — unlike `vmmc.py`, which
corrects wrong behaviour. `mort_mult=1.0` reproduces upstream exactly, so this
class is opt-in via `make_sim(hiv_class=...)` rather than the model default.

Kept in this repo (not as a patch to the editable stisim checkout) so a stisim
`git pull` cannot silently wipe it — which is what happened to the original
exp-005 VMMC patch.

Used by experiments/014_prior_expansion/.
"""

import numpy as np
import starsim as ss
import stisim as sti


class HIVMortalityMultiplier(sti.HIV):
    """HIV with a calibratable multiplier on the CD4-stratified death rates.

    The multiplier scales the **annual rate** before the per-timestep
    probability conversion, not the probability itself::

        p = peryear(base_rates * mort_mult).to_prob(dt)

    This matters: scaling the probability directly would let a multiplier above
    ~3 push the highest CD4 bin (0.30/yr) past 1.0 and silently saturate. Rate
    scaling stays a well-defined hazard at any positive multiplier.

    Everything else — CD4 dynamics, stage durations, transmission, ART — is
    inherited from upstream `sti.HIV` unchanged, so `mort_mult` is the only
    behavioural difference and `mort_mult=1.0` is an exact no-op.
    """

    # Upstream values from sti.HIV.make_p_hiv_death (stisim 1.5.8). Kept here so
    # a change upstream shows up as a diff rather than silently rescaling.
    CD4_BINS = np.array([1000, 500, 350, 200, 50, 0])
    BASE_RATES = np.array([0.003, 0.003, 0.005, 0.01, 0.05, 0.300])

    def __init__(self, *args, mort_mult=1.0, name='hiv', **kwargs):
        # mort_mult is consumed here rather than passed through: the parent's
        # update_pars() rejects keys not in HIVPars.
        super().__init__(*args, **kwargs)
        self.define_pars(mort_mult=float(mort_mult))
        # Keep the upstream name 'hiv' so this is a true drop-in — analyzers,
        # interventions and results all look the module up by that key.
        self.name = name

    def make_p_hiv_death(self, uids=None):
        """Per-timestep HIV death probability, with rates scaled by mort_mult."""
        rates = self.BASE_RATES * self.pars.mort_mult
        p_hiv_death = ss.peryear(rates).to_prob(self.dt)
        return p_hiv_death[np.digitize(self.cd4[uids], self.CD4_BINS)]
