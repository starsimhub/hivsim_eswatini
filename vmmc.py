"""Prevalence-target VMMC (in-repo model fix).

Upstream `sti.VMMC` in stisim 1.5.6–1.5.8 applies coverage as a per-step *hazard*
on the uncircumcised pool and ignores age stratification, so circumcision coverage
overshoots to ~100%. This subclass treats coverage as an age-stratified circumcision
*prevalence* (stock) target, matching SHIMS3/PHIA cross-sectional data. Validated in
experiments/015_vmmc_prevalence_target/ (coverage → SHIMS3 targets; male HIV
prevalence roughly doubles vs the broken upstream).

Kept in this repo (not as a patch to the editable stisim checkout) so a stisim
`git pull` cannot silently wipe it — which is what happened to the original
exp-005 patch. The fix has also been submitted upstream
(starsimhub/stisim `fix/vmmc-prevalence-target`); once stisim is bumped to a
release that includes it, this module can be removed and `make_interventions()`
can revert to plain `sti.VMMC`.
"""

import numpy as np
import starsim as ss
import stisim as sti
from stisim.interventions.utils import age_sex_mask
from stisim.utils import count


class VMMCPrevalenceTarget(sti.VMMC):
    """VMMC that hits an age-stratified circumcision *prevalence* target.

    Semantics (per timestep, per age/sex stratum):
      target_count   = coverage_proportion × (all alive males in the stratum)
      current_count  = males already circumcised in the stratum
      if current < target: circumcise the (target − current) highest-willingness
                           *uncircumcised* men in the stratum.
    Never removes (circumcision is irreversible). Everything else — the states
    (`circumcised`, `willingness`, `ti_circumcised`), coverage parsing, and the
    rel_sus efficacy reduction — is inherited from upstream `sti.VMMC` unchanged,
    so the *only* behavioural difference from upstream is which men are chosen.

    Assumes proportion ('p') coverage format (our `vmmc_coverage.csv` uses
    `p_vmmc`). Absolute-count ('n') format is not handled here.
    """

    def __init__(self, *args, name='vmmc', **kwargs):
        # Keep the upstream name 'vmmc' so this is a true drop-in — the
        # dashboard's VMMCPrevByAge analyzer looks the intervention up by that
        # key (sim.interventions['vmmc']).
        super().__init__(*args, **kwargs)
        self.name = name

    def _topup(self, stratum_mask, p, ti):
        """Top the stratum up to a p-fraction circumcised. Returns n newly done."""
        n_target = int(p * stratum_mask.count())              # stock target
        n_current = (stratum_mask & self.circumcised).count()  # already circ'd
        n_add = n_target - n_current
        if n_add <= 0:
            return 0
        candidates = (stratum_mask & ~self.circumcised).uids   # uncircumcised only
        if len(candidates) == 0:
            return 0
        n_add = min(n_add, len(candidates))
        # highest willingness first (same per-agent score upstream uses)
        pick = candidates[np.argsort(-self.willingness[candidates])[:n_add]]
        self.circumcised[pick] = True
        self.ti_circumcised[pick] = ti
        return len(pick)

    def step(self):
        sim = self.sim
        ppl = sim.people
        hiv = sim.diseases.hiv
        ti = self.ti

        n_new = 0
        if self.coverage is not None:
            # alive males ∧ user eligibility (if any). BoolArr must lead so the
            # result stays a BoolArr (uids & BoolArr would drop .count()).
            base = (ppl.male & ppl.alive) & self.check_eligibility()

            def cov_at(cov):
                return cov[ti] if len(cov) > ti else cov[-1]

            if isinstance(self.coverage, dict):
                # Stratified: hit each (age-bin[, sex]) target independently,
                # mirroring the key logic in compute_coverage_target.
                for ab in self.age_bins:
                    for sex in (self.sex_keys or [None]):
                        key = (ab, sex) if sex is not None else ab
                        p = cov_at(self.coverage.get(key, np.zeros(1)))
                        stratum = base & age_sex_mask(ab, 1, ppl)  # males in bin
                        n_new += self._topup(stratum, p, ti)
            else:
                # Aggregate single prevalence target across all eligible males
                n_new += self._topup(base, cov_at(self.coverage), ti)

        self.results['new_circumcisions'][ti] = n_new
        self.results['n_circumcised'][ti] = count(self.circumcised)

        # Efficacy: identical to upstream — reduction re-applied each step to all
        # circumcised (upstream relies on rel_sus being reset each step).
        hiv.rel_sus[self.circumcised] *= 1 - self.pars.eff_circ
        return
