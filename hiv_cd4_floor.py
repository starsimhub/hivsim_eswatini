"""CD4 floor during the acute phase (model fix, not a calibration knob).

Upstream `sti.HIV` has three CD4-decline functions. Two floor the result at the
stage's end value; `acute_decline` does not::

    acute_decline:     cd4 = self.cd4[uids] - per_timestep_decline
    falling_decline:   cd4 = np.maximum(cd4_end, self.cd4[uids] - per_timestep_decline)
    post_art_decline:  cd4[...] = np.maximum(cd4_end, ...)

`acute_decline` sizes a per-timestep decrement to carry CD4 from `cd4_start` to
`cd4_end` in `acute_dur` steps, where `acute_dur = ti_latent - ti_acute` is
"rounded to timesteps" (upstream's own comment). An agent that remains in the
acute state for one step more than that rounding assumed keeps subtracting the
same decrement, so CD4 falls past the latent set-point and can cross zero.

Negative CD4 is then fatal to the run rather than merely wrong, because
`make_p_hiv_death` looks the rate up positionally::

    rate = pars.cd4_death_rates[np.digitize(cd4, pars.cd4_death_bins)]

with **descending** bins `[1000, 500, 350, 200, 50, 0]` and a 6-element rate
array. For descending bins `np.digitize` returns `len(bins)` — here 6 — only
when the value is below the last edge, i.e. only for `cd4 < 0`. Every
non-negative CD4 maps to 0-5. So index 6 is unreachable unless CD4 goes
negative, and when it does the lookup raises

    IndexError: index 6 is out of bounds for axis 0 with size 6

`step_state` guards NaN CD4 ("Invalid entry for CD4") but not negative CD4.

**How it presented.** Experiment 025 crashed at wave 2 of the re-identification
pre-flight, 22 minutes in, after wave 1's 1000 points had completed cleanly.
None of the seven calibration parameters touch CD4 dynamics — `s_f_young` is
`rel_sus_age`, susceptibility, not survival — so this is not a bad region of
parameter space. It is exposure: wave 2 draws inside the NROY, which
concentrates on higher transmission, so more agents are infected and more pass
through the acute phase, and the run was at N = 20,000 against experiment 024's
10,000. 024 ran 1000 points without hitting it.

**It crashes rather than corrupting.** Any negative CD4 raises, on both the
off-ART path (`make_p_hiv_death`) and the on-ART path
(`get_art_mortality_hazard`, which does the same positional lookup). No result
produced before this fix can have been silently distorted by it.

Present in stisim 1.5.11 **and 1.6.1** — upgrading is not a way around it.

Kept in this repo rather than patched into the editable stisim checkout so a
`git pull` cannot silently wipe it, which is what happened to the original
exp-005 VMMC patch (see exp 015). Once the one-line fix is upstream, delete this
module and drop the `hiv_class` default in `run_sims.make_sim` — the same
retirement path `vmmc.VMMCPrevalenceTarget` took in exp 018.
"""

import numpy as np
import sciris as sc
import stisim as sti


class HIVCD4Floor(sti.HIV):
    """`sti.HIV` with the acute-phase CD4 decline floored at the latent set-point.

    The one behavioural difference from upstream is `np.maximum(cd4_end, ...)`
    in `acute_decline`, matching what `falling_decline` and `post_art_decline`
    already do. For every agent whose acute phase divides evenly into timesteps
    the floor never binds and this is an exact no-op; it changes behaviour only
    for the agents upstream would have driven below their latent set-point,
    which is the defect.
    """

    def acute_decline(self, uids):
        """Acute CD4 decline, floored at the latent set-point.

        Upstream body, with the floor added. `cd4_end` is the agent's latent
        set-point, so flooring here means "the acute decline stops where the
        latent phase begins" — it cannot overshoot into the falling phase's
        range, let alone below zero.
        """
        acute_start = self.ti_acute[uids]   # Time of infection
        acute_end = self.ti_latent[uids]    # Time to latent infection
        acute_dur = acute_end - acute_start  # Total time in acute phase, rounded to timesteps
        cd4_start = self.cd4_start[uids]
        cd4_end = self.cd4_latent[uids]
        per_timestep_decline = sc.safedivide(cd4_start - cd4_end, acute_dur)
        cd4 = np.maximum(cd4_end, self.cd4[uids] - per_timestep_decline)
        return cd4


# Note for anyone reusing hiv_mortality.HIVMortalityMultiplier (exp 014): it
# subclasses sti.HIV directly, so it does not inherit this floor, and its own
# make_p_hiv_death repeats the same 6-element positional lookup. If that class
# is ever brought back into use, give it this floor too — e.g. by having it
# inherit from HIVCD4Floor instead of sti.HIV.
