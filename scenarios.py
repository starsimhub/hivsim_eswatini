"""Scenario definitions for the decision analysis.

Sits beside `interventions.py`, which builds the *calibration* baseline. This
module builds the things layered on top of it, and is passed to
`run_sims.make_sim(extra_interventions=...)` -- the hook at run_sims.py:88 that
appends after the defaults, so anything re-targeting a stock runs last.

The research question (CLAUDE.md)
--------------------------------
Relative value of scaling up long-acting PrEP versus improving the treatment
cascade to raise population viral suppression -- and **where prevention still
matters when suppression is high.** That last clause is the novel part, and it
is a subtraction: (PrEP + high suppression) minus (high suppression alone).

Why PrEP is not in the calibration baseline
-------------------------------------------
Exp 018 removed it. `sti.Prep()` with `coverage=None` does not mean "off" -- it
falls back to a built-in ramp reaching 80% of FSW by 2025, **starting in 2004**,
a decade before PrEP had efficacy evidence. Every experiment from 001 to 017 ran
with that undeclared programme. `interventions.py` therefore ships no PrEP, and
every scenario here specifies `coverage` explicitly. In a counterfactual
comparison a phantom baseline intervention is worse than a phantom fitted one:
it biases the contrast, not just the fit.

Why the baseline suppression is trustworthy, which the comparison depends on
---------------------------------------------------------------------------
stisim's `vls_coverage=None` means 100% of ART initiators are virally
suppressed, which exp 021 measured as overstating population suppression by
8.8 percentage points in 2016. Measuring a suppression-raising scenario against
a baseline that already overstates suppression **understates what cascade
improvement buys, and so systematically favours PrEP.** The PHIA-derived
`vls_coverage` plus the exp 022 stock target is what makes this comparison
valid; it is not a detail, and it belongs in the abstract's methods.
"""

import numpy as np
import pandas as pd
import starsim as ss
import stisim as sti

# --- Lenacapavir ---------------------------------------------------------------
# From Akullian et al. LenOptim manuscript (Lancet HIV, R1), Methods and
# Table 1: efficacy base case 0.95, DSA range 0.90-1.00 (refs 29, 30);
# 6-monthly dosing, with the pppy cost covering the two doses needed for
# continuous cover over one year.
LEN_EFF = 0.95
LEN_EFF_RANGE = (0.90, 1.00)
LEN_DUR_MONTHS = 6

# --- Eligibility functions -----------------------------------------------------
# `sti.Prep` applies HIV-negative and not-already-on-PrEP filters on top, so
# these only express *who to target*.
#
# IMPORTANT, and it changes the answer: within an eligible group `sti.Prep`
# fills its coverage target by ranking on `willingness`, a RANDOM per-agent
# score (hiv_interventions.py:830). So by default PrEP reaches a random subset
# of the eligible -- relative incidence in users ~1x the group average.
#
# The LenOptim manuscript's uptake-heterogeneity scenario analysis found the
# >2x quartile "most representative of self-selecting PrEP adopters consistent
# with the ECHO PrEP implementation study". Random allocation within a broad
# group therefore UNDERSTATES impact relative to a real programme. The way to
# express self-selection with the existing API is to narrow eligibility onto
# the higher-risk strata rather than to widen coverage -- hence the
# `_risk` variants below. Reporting both is the honest treatment, because the
# difference between them is a programme-delivery assumption, not biology.


def agyw(sim):
    """Adolescent girls and young women, 15-24."""
    ppl = sim.people
    return ppl.female & (ppl.age >= 15) & (ppl.age < 25)


def agyw_risk(sim):
    """AGYW in the higher-activity risk groups -- a self-selecting programme."""
    net = sim.networks.structuredsexual
    return agyw(sim) & (net.risk_group >= 1)


def fsw(sim):
    """Female sex workers."""
    return sim.networks.structuredsexual.fsw


def women_15_34(sim):
    """Broader female target: 15-34, where incidence peaks in this model."""
    ppl = sim.people
    return ppl.female & (ppl.age >= 15) & (ppl.age < 35)


ELIGIBILITY = {"agyw": agyw, "agyw_risk": agyw_risk,
               "fsw": fsw, "women_15_34": women_15_34}


# --- PrEP ----------------------------------------------------------------------

def lenacapavir(coverage, eligibility="agyw", start=2026, scale_to=2030,
                eff=LEN_EFF, name="prep_len"):
    """6-monthly lenacapavir as a coverage (stock) target.

    Args:
        coverage:    final coverage, as a PROPORTION of the eligible group.
                     stisim treats this as a prevalence target, so it is the
                     fraction covered at a given time, not a fraction initiated.
        eligibility: key into ELIGIBILITY, or a callable (sim) -> BoolArr
        start:       first year of any coverage. Default 2026 -- the scenario
                     window opens after the last calibration target (2016) and
                     after the 2021 validation hold-out, so no scenario can
                     contaminate either.
        scale_to:    year the final coverage is reached (linear from `start`)
        eff:         efficacy; vary over LEN_EFF_RANGE for sensitivity
    """
    elig = ELIGIBILITY[eligibility] if isinstance(eligibility, str) else eligibility
    return sti.Prep(
        name=name,
        prep_eff=eff,
        prep_dur=ss.months(LEN_DUR_MONTHS),
        prep_adh=1.0,          # adherence is not the constraint for an injectable;
                               # discontinuation is, and that is `coverage` falling
        # Anchors start AT `start`, not at start-1. With a zero anchor at
        # start-1 and a non-zero one at `start`, stisim interpolates across the
        # intervening year's twelve monthly steps, so coverage is already ~90%
        # of the `start` value by that December -- exp 026's verification caught
        # 1,601 people on PrEP in 2025 for a scenario documented as starting in
        # 2026. Small (0.3% of PrEP person-years, and outside the outcome
        # window) but the design should be literally true.
        coverage={"year": [start, start + 1, scale_to],
                  "value": [0.0, coverage * 0.1, coverage]},
        eligibility=elig,
    )


# --- Treatment cascade ---------------------------------------------------------

def _ramp_long(df, value_col, target, start, reach, stop=2040):
    """Ramp `value_col` linearly toward `target` between `start` and `reach`.

    Both coverage tables are LONG format -- (Year, AgeBin, Gender, value) with
    one row per stratum-year, not year columns. So a scenario cannot just scale
    existing columns: the future years have to be ADDED as rows, because the
    calibration tables stop at the last observation (2021-ish) and stisim
    forward-fills from the final year. Without new rows a scenario would simply
    inherit the last observed value forever, and every arm would be identical.

    That failure would have been silent -- the sim runs, the arms differ only by
    PrEP, and the cascade axis quietly does nothing.
    """
    df = df.copy()
    strata = [c for c in df.columns if c not in ("Year", value_col)]
    last_year = int(df.Year.max())
    base = (df[df.Year == last_year]
            .drop(columns="Year")
            .set_index(strata)[value_col])

    # A scale-up RAISES strata that are below the target and LEAVES ALONE those
    # already above it. Ramping every stratum toward a flat target would
    # un-treat people in the strata that have already exceeded it -- Eswatini's
    # ART coverage is 0.650-0.970 across strata, so a flat ramp to 0.95 would
    # pull three of eight strata DOWN while raising five. Hence the per-stratum
    # max().
    if np.isfinite(target):
        eff_target = np.maximum(base.to_numpy(dtype=float), target)
    else:
        eff_target = base.to_numpy(dtype=float)

    # Guard: if NO stratum is below the target, this is not a scale-up at all.
    # Not hypothetical -- Eswatini has already achieved the third 95 (SHIMS3
    # 2021: suppression given ART is 0.959 in women, 0.967 in men), so
    # cascade_target(target=0.95) is a no-op at best and a regression if the
    # max() above were removed. Failing loudly beats an arm that silently does
    # nothing and gets reported as "cascade improvement".
    if np.isfinite(target) and not (base.to_numpy(dtype=float) < target).any():
        raise ValueError(
            f"target {target} is at or below the {last_year} baseline for ALL "
            f"{len(base)} strata (min {base.min():.3f}, max {base.max():.3f}). "
            f"This arm would do nothing. For Eswatini the third 95 is already "
            f"met -- the cascade headroom is in ART coverage, not in "
            f"suppression among the treated. Pick a higher target or ramp the "
            f"other axis.")

    rows = []
    for y in range(max(int(start), last_year + 1), int(stop) + 1):
        w = min((y - start) / max(reach - start, 1e-9), 1.0)
        w = max(w, 0.0)
        for (key, v0), tgt in zip(base.items(), eff_target):
            key = key if isinstance(key, tuple) else (key,)
            rows.append({**dict(zip(strata, key)), "Year": y,
                         value_col: float(v0) * (1 - w) + float(tgt) * w})
    return pd.concat([df, pd.DataFrame(rows)], ignore_index=True)


def cascade_target(base_vls, target=0.95, start=2026, reach=2030, stop=2040):
    """Raise the third 95 -- suppression among those on ART -- toward `target`.

    Returns a table for `make_sim(art_vls_coverage=...)`. Because exp 022 made
    suppression a STOCK target, this reaches patients already on treatment --
    the mechanism that actually drove Eswatini's third-95 gains (the TLD
    transition from ~2019). Without that fix this arm would only affect new
    initiators and would badly understate cascade improvement.

    Moves the third 95 ONLY. The second 95 (diagnosis -> ART) is a separate
    axis; see `art_coverage_target`.
    """
    return _ramp_long(base_vls, "p_vls", target, start, reach, stop)


def art_coverage_target(base_art, target=0.95, start=2026, reach=2030,
                        stop=2040):
    """Raise ART coverage among PLHIV -- the second 95."""
    return _ramp_long(base_art, "p_art", target, start, reach, stop)


def baseline_tables():
    """The calibration baseline's two coverage tables, for scenarios to modify."""
    import pandas as pd_
    from vls_construction import build, to_vls_coverage
    vls = to_vls_coverage(build(), fill_back_to=1985)
    art = pd_.read_csv("data/art_coverage.csv")
    return art, vls


def hold_flat(df, value_col, stop=2040):
    """Extend a table to `stop` at its last observed value -- the status quo arm.

    Needed for symmetry: the cascade arms add future rows, so the baseline must
    too, or the arms would differ in how far their tables extend as well as in
    their values. stisim forward-fills, so this is a no-op numerically -- but it
    makes the comparison explicit rather than relying on that behaviour.
    """
    return _ramp_long(df, value_col, target=np.nan, start=stop + 1,
                      reach=stop + 1, stop=stop).assign(
        **{value_col: lambda d: d[value_col].ffill()})
