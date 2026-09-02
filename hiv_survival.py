"""Age-dependent untreated survival, and deaths split by route.

Two things upstream `sti.HIV` does not provide, both needed by experiment 019.

**1. An age gradient on untreated survival.** Untreated survival in stisim is
`dur_acute + dur_latent + dur_falling` = lognorm_ex(10 y, 3 y) latency plus
lognorm_ex(3 y, 1 y) late stage, mean ~13.1 y, applied identically to a
17-year-old and a 55-year-old. Real untreated survival falls steeply with age
at seroconversion (ALPHA-network pooled: ~12.5 y at 15-24 down to ~7.5 y at
45+). `latent_mult` rescales the latent interval as a function of age at
infection, leaving the drawn `dur_falling` untouched so the CD4 decline rate --
and hence `rel_trans_falling`, the 8x late-infection transmissibility -- is
unchanged. Only the latency stretches or compresses.

**2. Deaths split by route.** The HIV module kills agents two ways and
`hiv.new_deaths` pools them:

  - `ti_zero`      death date drawn at infection; tunable by nothing
  - `p_hiv_death`  per-timestep Bernoulli on CD4; tunable by rel_death,
                   rel_death_f, and the on-ART multipliers

017 observation 7 estimated the split at ~80/20 by integrating the hazard
table, and every mortality knob in stisim acts on the minority route -- which
is why 017's two candidate explanations for the AIDS-death deficit were both
inert. That estimate has never been measured. This class measures it.

Two properties worth preserving, both relied on by 019's design:

  - **`latent_mult=1.0` is an exact no-op.** The rescaling consumes no random
    numbers -- it rescales intervals already drawn by `super()` -- so arm A is
    bit-identical to upstream at the same seed, and the instrumentation is
    free. Contrast `hiv_mortality.py`, which also leaves dynamics unchanged at
    `mort_mult=1.0` but for the simpler reason that it multiplies by one.
  - **Under-15s are never rescaled, in any arm.** The ALPHA estimates are for
    adult seroconverters; vertically infected infants have a different natural
    history entirely (real untreated survival ~2 y, which stisim already
    overstates badly). Fixing that is out of scope for 019 and is not silently
    bundled in here. Because the exclusion is identical across arms, the
    arm-to-arm contrasts are unaffected.

Kept in this repo, not as a patch to the editable stisim checkout, so a stisim
`git pull` cannot silently wipe it -- which is what happened to the original
exp-005 VMMC patch.

Used by experiments/019_age_dependent_survival/.
"""

import numpy as np
import pandas as pd
import starsim as ss
import stisim as sti


class AgeDependentSurvival(sti.HIV):
    """HIV with an age gradient on untreated latency, and deaths split by route.

    Args:
        latent_mult: multiplier on the latent interval. Either a scalar applied
            to everyone aged >= 15, or a sequence of one value per band in
            `age_bands`. 1.0 reproduces upstream exactly.
        age_bands: (lo, hi) age-at-infection bands, left-closed. Defaults to
            ALPHA's 15-24 / 25-34 / 35-44 / 45+.
    """

    AGE_BANDS = ((15, 25), (25, 35), (35, 45), (45, 200))
    MIN_AGE = 15  # below this, never rescale -- see module docstring

    def __init__(self, *args, latent_mult=1.0, age_bands=None, name='hiv',
                 **kwargs):
        # latent_mult is consumed here rather than passed through: the parent's
        # update_pars() rejects keys not in HIVPars. Same reason as
        # hiv_mortality.py's mort_mult.
        super().__init__(*args, **kwargs)
        self.age_bands = tuple(age_bands) if age_bands is not None else self.AGE_BANDS

        mult = np.atleast_1d(np.asarray(latent_mult, dtype=float))
        if mult.size not in (1, len(self.age_bands)):
            errormsg = (f'latent_mult must be a scalar or one value per age band '
                        f'({len(self.age_bands)}), got {mult.size}')
            raise ValueError(errormsg)
        if (mult <= 0).any():
            raise ValueError(f'latent_mult must be positive, got {latent_mult}')
        self.latent_mult = mult

        # Age at infection is needed both to pick the multiplier and to report
        # realized survival by band. Upstream clears ti_infected in step_die(),
        # so this cannot be reconstructed after the fact.
        self.define_states(ss.FloatArr('age_at_infection'))

        # Keep the upstream name 'hiv' so this is a true drop-in -- analyzers,
        # interventions and results all look the module up by that key.
        self.name = name

        self._death_log = []

    # --- Age gradient --------------------------------------------------------

    def multiplier_for_age(self, age):
        """Latent-interval multiplier per agent, given age at infection."""
        mult = np.ones(len(age), dtype=float)
        if self.latent_mult.size == 1:
            mult[age >= self.MIN_AGE] = self.latent_mult[0]
        else:
            for (lo, hi), value in zip(self.age_bands, self.latent_mult):
                mult[(age >= max(lo, self.MIN_AGE)) & (age < hi)] = value
        return mult

    def set_prognoses(self, uids, sources=None, ti=None):
        """Upstream prognoses, then rescale latency by age at infection.

        Consumes no random numbers of its own: `super()` draws dur_acute,
        dur_latent and dur_falling, and this only rescales the interval between
        ti_latent and ti_falling, sliding ti_zero along with it so the drawn
        dur_falling is preserved exactly.
        """
        super().set_prognoses(uids, sources=sources, ti=ti)

        age = self.sim.people.age[uids]
        self.age_at_infection[uids] = age

        if self.latent_mult.size == 1 and self.latent_mult[0] == 1.0:
            return  # exact no-op, and skips the arithmetic entirely

        mult = self.multiplier_for_age(np.asarray(age))
        ti_latent = self.ti_latent[uids]
        dur_latent = self.ti_falling[uids] - ti_latent
        dur_falling = self.ti_zero[uids] - self.ti_falling[uids]

        self.ti_falling[uids] = ti_latent + np.rint(dur_latent * mult)
        self.ti_zero[uids] = self.ti_falling[uids] + dur_falling
        return

    # --- Instrumentation -----------------------------------------------------

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_deaths_progression', dtype=int,
                      label='HIV deaths via ti_zero'),
            ss.Result('new_deaths_hazard', dtype=int,
                      label='HIV deaths via p_hiv_death'),
        )
        return

    def step_die(self, uids):
        """Record route and realized survival, then let upstream clear states.

        This has to happen before `super().step_die()`, which sets ti_infected
        and ti_zero to NaN -- after that the route is unrecoverable.

        Upstream sets ti_dead from the Bernoulli hazard first and then
        overwrites it for anyone whose ti_zero has arrived, so an agent can
        satisfy both routes on the same step. Ties go to `progression`: those
        agents were scheduled to die this step regardless of the coin flip.
        """
        ti = self.ti
        if len(uids):
            hiv_dead = uids[self.ti_dead[uids] == ti]
            if len(hiv_dead):
                ti_zero = self.ti_zero[hiv_dead]
                via_prog = np.isfinite(ti_zero) & (ti_zero <= ti)

                self.results['new_deaths_progression'][ti] += int(via_prog.sum())
                self.results['new_deaths_hazard'][ti] += int((~via_prog).sum())

                survival_y = (ti - self.ti_infected[hiv_dead]) * self.sim.t.dt_year
                self._death_log.append(np.column_stack([
                    np.full(len(hiv_dead), ti, dtype=float),
                    self.age_at_infection[hiv_dead],
                    np.asarray(survival_y, dtype=float),
                    via_prog.astype(float),
                    self.on_art[hiv_dead].astype(float),
                    self.sim.people.female[hiv_dead].astype(float),
                ]))

        super().step_die(uids)
        return

    def death_records(self):
        """Agent-level death records as a DataFrame; empty if nobody died.

        Columns: ti, age_at_infection, survival_years, via_progression, on_art,
        female. Age at death is age_at_infection + survival_years.

        Agent-level, so *unscaled* -- unlike the results arrays, which starsim
        multiplies by pop_scale. Use for survival distributions (019 metric 3);
        scale by pop_scale before comparing counts to UNAIDS or to 016's
        implied AIDS deaths by age.
        """
        cols = ['ti', 'age_at_infection', 'survival_years', 'via_progression',
                'on_art', 'female']
        if not self._death_log:
            return pd.DataFrame(columns=cols)
        df = pd.DataFrame(np.vstack(self._death_log), columns=cols)
        for col in ('via_progression', 'on_art', 'female'):
            df[col] = df[col].astype(bool)
        return df


# The four arms of experiment 019. Bands are 15-24 / 25-34 / 35-44 / 45+.
#
# B, C and D all sit at ~11.5 y mean survival by construction -- 0.84 is
# approximately the population mean of the ALPHA gradient once weighted by the
# model's age-at-infection distribution -- so they differ in *shape*, not
# level. That makes B->D the gradient effect at constant level, and A->B the
# level effect at zero gradient. Realized means are measured per arm rather
# than assumed, since the age distribution of incident infections is endogenous
# and shifts between arms.
ARMS = {
    'A_flat_13':    1.0,
    'B_flat_11.5':  0.84,
    'C_grad_mild':  (0.89, 0.84, 0.73, 0.61),
    'D_grad_alpha': (0.94, 0.84, 0.64, 0.44),
    # Derived from the EMOD model of Akullian et al. 2020 (Lancet HIV), whose
    # appendix Table A2 gives age-dependent untreated survival as a Weibull with
    # shape 2 and scale lambda = 21.182 - 0.2717 * age_at_infection. Mean
    # survival = lambda * Gamma(1.5), so 13.96 y at age 20, 11.55 at 30, 9.14 at
    # 40 and 6.73 at 50. Converting to a latency multiplier at the band
    # midpoints, holding stisim's non-latent portion (dur_acute 1.7 months +
    # dur_falling 3 y = 3.14 y) fixed:  mult = (mean - 3.14) / 10.
    #
    # The result PIVOTS rather than compresses: >1 at the young end, far below
    # 1 at the old end. Every arm in exp 019 had multipliers <= 1.0, so 019
    # tested only shortening gradients and its conclusion that a gradient is not
    # worth a parameter was drawn from too narrow a family. A pivot pushes both
    # of the model's prevalence defects the right way at once -- longer young
    # survival adds prevalent infections where the model is 49-65% low, shorter
    # old survival removes the excess where women 35-44 sit 11% high.
    'E_grad_emod':  (1.086, 0.845, 0.604, 0.363),
}
