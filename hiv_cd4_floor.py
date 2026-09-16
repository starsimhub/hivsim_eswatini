"""Non-negative CD4 invariant (model fix, not a calibration knob).

## The crash

`sti.HIV` looks the per-timestep death rate up positionally::

    rate = pars.cd4_death_rates[np.digitize(cd4, pars.cd4_death_bins)]

`cd4_death_bins` is **descending** — [1000, 500, 350, 200, 50, 0] — against a
6-element `cd4_death_rates`. For descending bins `np.digitize` returns
`len(bins)` only when the value falls below the last edge, so **index 6 is
reachable only for `cd4 < 0`**; every non-negative CD4 maps to 0-5. Verified
across the range. When it happens::

    IndexError: index 6 is out of bounds for axis 0 with size 6

`step_state` guards NaN CD4 ("Invalid entry for CD4") but not negative CD4.
Both the off-ART path (`make_p_hiv_death`) and the on-ART path
(`get_art_mortality_hazard`) do the same positional lookup, so both raise.

This crashes rather than corrupting: any negative CD4 raises immediately, so no
completed run can have been silently distorted by it.

## Two defects, one of which is still unidentified

**Defect 1 — the acute decline had no floor (found, fixed, insufficient).**
`acute_decline` was the only one of the three CD4-decline functions without a
floor at its stage's end value::

    acute_decline     cd4 = self.cd4[uids] - per_timestep_decline          # no floor
    falling_decline   cd4 = np.maximum(cd4_end, self.cd4[uids] - ...)      # floored
    post_art_decline  cd4[...] = np.maximum(cd4_end, ...)                  # floored

It sizes a decrement to carry CD4 from `cd4_start` to `cd4_end` in `acute_dur`
steps, where `acute_dur` is rounded to timesteps. An agent staying acute one
step longer keeps subtracting the same decrement. Demonstrated on real acute
agents: one extra step overshoots the set-point by ~350 CD4, the second crosses
zero.

**Defect 2 — unidentified.** Exp 027 ran with defect 1 fixed (verified: the fix
commits are ancestors of the launched tree, and `make_sim` defaulted to this
class) and **crashed with the identical IndexError**. So the acute phase is not
the only route to a negative CD4. `falling_decline` and `post_art_decline` both
floor at 1, and the latent assignment is direct, so the remaining suspects are
`cd4_increase` (the on-ART logistic, which can go negative if `dur_art < 0` or
if `cd4_potential < 2 * cd4_preart`) and any path that leaves CD4 unwritten
while state flags move. **Not yet proven — that is what the logging below is
for.**

## What this class does

1. **Floors the acute decline** (defect 1), matching the other two functions.
2. **Enforces the invariant at both lookups** (defect 2 and anything else):
   clips CD4 to `CD4_FLOOR` immediately before the rate lookup, so no path can
   crash the run regardless of how the value was produced.
3. **Records every clip**, with the offending agents' state flags, so the real
   source is identified from data rather than guessed at. One run both completes
   and diagnoses.

`CD4_FLOOR = 1.0` is not invented: it is the `cd4_end = 1  # To avoid divide by
zero problems` that upstream's own `falling_decline` and `post_art_decline` use
as their floor.

**The clip is a repair of nonphysical state, and it is loud by design.** A
negative CD4 is not a small numerical wobble — the agent's death hazard is
undefined. If `guard_events` is non-empty in a production run, that is a defect
to fix upstream, not a result to accept quietly.

Kept in this repo rather than patched into the editable stisim checkout so a
`git pull` cannot silently wipe it, which is what happened to the original
exp-005 VMMC patch (see exp 015). Once upstream carries the invariant, delete
this module and the `hiv_class` default in `run_sims.make_sim` — the retirement
path `vmmc.VMMCPrevalenceTarget` took in exp 018.
"""

import os
import sys

import numpy as np
import sciris as sc
import stisim as sti

#: Floor applied to nonphysical CD4. Matches upstream's own `cd4_end = 1`.
CD4_FLOOR = 1.0

#: Cap on stderr diagnostic lines per process, so a systematic defect cannot
#: flood a 1000-point wave's log. Structured events are still recorded in full
#: on the module instance regardless of this cap.
MAX_REPORTS_PER_PROCESS = 10

_n_reported = 0


class HIVCD4Floor(sti.HIV):
    """`sti.HIV` with a non-negative CD4 invariant, plus clip diagnostics.

    Three differences from upstream, all narrow:

    - `acute_decline` floors at the latent set-point (defect 1).
    - `make_p_hiv_death` and `get_art_mortality_hazard` clip CD4 to
      `CD4_FLOOR` before their positional rate lookup (defect 2 and any other
      path).
    - Each clip is recorded in `self.guard_events`.

    Where no clip is needed and the acute floor does not bind, this is an exact
    no-op: an identical seeded run gives the same `cum_infections` as upstream.
    """

    def __init__(self, *args, name='hiv', **kwargs):
        # Keep the upstream name 'hiv' so this is a true drop-in. Without it
        # starsim names the module after the class and every intervention doing
        # `sim.diseases.hiv` raises AttributeError — the same reason
        # HIVMortalityMultiplier pins its name.
        super().__init__(*args, **kwargs)
        self.name = name
        #: One dict per clip event. Read by the experiment's run.py.
        self.guard_events = []

    # --- Defect 1: the missing floor -------------------------------------------

    def acute_decline(self, uids):
        """Acute CD4 decline, floored at the latent set-point.

        Upstream body with the floor added. `cd4_end` is the agent's latent
        set-point, so the acute decline stops where the latent phase begins and
        cannot overshoot into the falling phase's range, let alone below zero.
        """
        acute_start = self.ti_acute[uids]   # Time of infection
        acute_end = self.ti_latent[uids]    # Time to latent infection
        acute_dur = acute_end - acute_start  # Total time in acute phase, rounded to timesteps
        cd4_start = self.cd4_start[uids]
        cd4_end = self.cd4_latent[uids]
        per_timestep_decline = sc.safedivide(cd4_start - cd4_end, acute_dur)
        cd4 = np.maximum(cd4_end, self.cd4[uids] - per_timestep_decline)
        return cd4

    # --- Defect 2 and everything else: the invariant ---------------------------

    def _enforce_cd4_floor(self, uids, where):
        """Clip nonphysical CD4 to `CD4_FLOOR`, recording what was clipped.

        Runs immediately before a positional rate lookup. Returns the number of
        agents clipped (0 in the overwhelming majority of steps, so the cost is
        one comparison over the passed uids).

        The state flags are the point of this method. Knowing *that* CD4 went
        negative does not locate the defect; knowing the agents were, say,
        on-ART with a negative `dur_art` does.
        """
        global _n_reported

        if uids is None or len(uids) == 0:
            return 0

        vals = self.cd4[uids]
        bad = np.nonzero(vals < CD4_FLOOR)[0]
        if len(bad) == 0:
            return 0

        bad_uids = uids[bad]

        def _frac(state):
            """Fraction of the clipped agents carrying a given boolean state."""
            try:
                return float(np.count_nonzero(np.asarray(state[bad_uids]))) / len(bad_uids)
            except Exception:  # a state may not exist across stisim versions
                return float('nan')

        event = dict(
            where=where,
            ti=int(self.ti),
            n=int(len(bad_uids)),
            min_cd4=float(np.min(vals[bad])),
            median_cd4=float(np.median(vals[bad])),
            f_acute=_frac(self.acute),
            f_latent=_frac(self.latent),
            f_falling=_frac(self.falling),
            f_on_art=_frac(self.on_art),
            f_art_discontinued=_frac(self.art_discontinued),
        )
        # dur_art < 0 is the leading hypothesis for defect 2: the on-ART
        # logistic in cd4_increase evaluates to 2*cd4_preart - cd4_potential as
        # dur_art -> -inf, which is negative whenever cd4_potential exceeds
        # twice cd4_preart.
        try:
            dur_art = np.asarray(self.ti - self.ti_art[bad_uids], dtype=float)
            finite = dur_art[np.isfinite(dur_art)]
            event['min_dur_art'] = float(np.min(finite)) if finite.size else float('nan')
            event['f_dur_art_negative'] = (
                float(np.count_nonzero(finite < 0)) / finite.size if finite.size else float('nan')
            )
        except Exception:
            event['min_dur_art'] = float('nan')
            event['f_dur_art_negative'] = float('nan')

        self.guard_events.append(event)

        if _n_reported < MAX_REPORTS_PER_PROCESS:
            _n_reported += 1
            print(
                f"[cd4-guard] pid={os.getpid()} ti={event['ti']} where={where} "
                f"n={event['n']} min_cd4={event['min_cd4']:.3f} "
                f"acute={event['f_acute']:.2f} latent={event['f_latent']:.2f} "
                f"falling={event['f_falling']:.2f} on_art={event['f_on_art']:.2f} "
                f"art_disc={event['f_art_discontinued']:.2f} "
                f"min_dur_art={event['min_dur_art']:.2f} "
                f"f_dur_art_neg={event['f_dur_art_negative']:.2f}",
                file=sys.stderr, flush=True,
            )

        self.cd4[bad_uids] = CD4_FLOOR
        return len(bad_uids)

    def make_p_hiv_death(self, uids=None):
        """Off-ART death probability, with the CD4 invariant enforced first."""
        self._enforce_cd4_floor(uids, 'make_p_hiv_death')
        return super().make_p_hiv_death(uids)

    def get_art_mortality_hazard(self, uids):
        """On-ART death probability, with the CD4 invariant enforced first."""
        self._enforce_cd4_floor(uids, 'get_art_mortality_hazard')
        return super().get_art_mortality_hazard(uids)


# Note for anyone reusing hiv_mortality.HIVMortalityMultiplier (exp 014): it
# subclasses sti.HIV directly, so it inherits neither the floor nor the
# invariant, and its own make_p_hiv_death repeats the same 6-element positional
# lookup. If that class is brought back into use, have it inherit from
# HIVCD4Floor instead of sti.HIV.
