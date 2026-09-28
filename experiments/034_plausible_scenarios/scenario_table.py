"""Exp 034 -- the full scenario x stratum table, for writing from.

One row per scenario x sex x broad age group, carrying the cascade, the
implied time spent at each cascade step, the infection burden, and incidence
against the status-quo baseline.

Cascade residence times -- mean years undiagnosed, diagnosed-but-untreated, and
on-ART-but-unsuppressed, via Little's Law (W = stock / throughput, with
incidence as the throughput) -- were built and then REMOVED at the researcher's
request.

They are straightforward to reinstate, but do not do so without handling three
biases, all of which inflate them and none of which is correctable from these
outputs alone:

  1. The epidemic is NOT in steady state -- incidence is falling throughout, so
     the denominator understates the historical flow that built each stock.
     Exp 030 hit the same bias and measured it as substantial.
  2. People LEAVE these compartments by dying as well as by progressing, so
     throughput exceeds new infections and residence time is overstated.
  3. Stocks are counted in the band a person is in NOW; infections are counted
     in the band they were in when infected. Ageing across band boundaries
     therefore breaks the per-stratum accounting, worst in the oldest group,
     which accumulates people who were infected decades earlier and younger.

Usage (repo root):  python experiments/034_plausible_scenarios/scenario_table.py
"""

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))

OUT = HERE / "outputs"
SIMS = OUT / "sims"
WINDOW = (2026, 2040)
CASC_YEAR = 2030          # cascade and residence times are read here

AGE_GROUPS = [(15, 25, "15-24"), (25, 50, "25-49"), (50, 200, "50+"),
              (15, 200, "15+ (all)")]

CASCADES = {
    "S0_status_quo": "Status quo. HIV testing, ART coverage and viral "
        "suppression held at their 2026 values throughout. Represents current "
        "policy continuing unchanged, not a frozen epidemic -- coverage among "
        "PLHIV still drifts as the PLHIV population ages.",
    "S1_testing_only": "Testing only. The routine HIV testing rate is tripled "
        "from 2026, reaching 0.75/yr among undiagnosed adults by 2030, with NO "
        "expansion of ART coverage. Isolates diagnosis from treatment.",
    "S2_unaids_95": "95-95-95 in every group. Testing tripled, and the ART "
        "coverage target raised to 90.25% of PLHIV (0.95 aware x 0.95 linked) "
        "in every age-sex stratum. Targets are per stratum, so laggards are "
        "raised and leaders untouched. Viral suppression is left at baseline: "
        "Eswatini's measured 96.2% already exceeds the third 95.",
    "S3_99_96_98": "Maximal plausible cascade. Testing quadrupled, ART "
        "coverage target 96% of PLHIV, viral suppression target 98% -- each at "
        "or below the best value any single age-sex group in SHIMS3 has "
        "achieved, but required here of every group simultaneously.",
}

PREPS = {
    "P0_none": "No PrEP.",
    "P1_fsw": "LA-PrEP, 60% coverage of female sex workers.",
    "P2_agyw_risk": "LA-PrEP: FSW at 60%, plus 30% of higher-activity AGYW "
        "aged 15-24.",
    "P3_agyw_all": "LA-PrEP: FSW at 60%, plus 30% of all AGYW aged 15-24.",
    "P4_women_25_34": "LA-PrEP: FSW at 60%, plus 30% of all women aged 15-34. "
        "The broad-population arm.",
}


def cols_for(df, stem, sex, lo, hi):
    """Columns of `stem` whose own band lies inside [lo, hi).

    Read from the data, never constructed: popagesex runs 5-year bands to
    95-100 while cascadeage lumps everything above 80 into one 80-200 band, so
    a name built from one assumed grid silently misses on the other.
    """
    out, pre = [], f"{stem}_{sex}_"
    for c in df.columns:
        if not c.startswith(pre):
            continue
        try:
            a, b = int(c.rsplit("_", 2)[-2]), int(c.rsplit("_", 2)[-1])
        except ValueError:
            continue
        if a >= lo and (b <= hi or hi >= 200):
            out.append(c)
    return out


def load(cascade, prep):
    f = glob.glob(str(SIMS / f"{cascade}__{prep}__*.parquet"))
    if not f:
        return None
    return pd.concat([pd.read_parquet(x) for x in f], ignore_index=True)


def stratum_metrics(d, sex, lo, hi):
    """Everything for one scenario x sex x age group."""
    S = lambda stem: d[cols_for(d, stem, sex, lo, hi)].sum(axis=1)  # noqa: E731
    plhiv = S("cascadeage.n_infected")
    dx = S("cascadeage.n_diagnosed")
    art = S("cascadeage.n_on_art")
    vls = S("cascadeage.n_effective_art")
    newi = S("popagesex.new_infections")
    alive = S("popagesex.n_alive")
    susc = alive - S("popagesex.n_infected")

    t = pd.DataFrame({"year": d.timevec.values, "seed": d.seed.values,
                      "plhiv": plhiv.values, "dx": dx.values, "art": art.values,
                      "vls": vls.values, "newi": newi.values,
                      "susc": susc.values})
    # Pool numerator and denominator across seeds before dividing: the mean of
    # per-seed ratios is biased in small strata.
    g = t.groupby("year").sum(numeric_only=True)
    nseed = t.seed.nunique()

    def at(y, col):
        return g.loc[y, col] if y in g.index else np.nan

    def rate(y, num, den):
        n, dd = at(y, num), at(y, den)
        return n / dd if dd and np.isfinite(dd) and dd > 0 else np.nan

    inc = {y: 100 * at(y, "newi") / at(y, "susc") for y in (2030, 2040)}
    lo_y, hi_y = WINDOW
    cum = g.loc[(g.index >= lo_y) & (g.index <= hi_y), "newi"].sum() / nseed

    y = CASC_YEAR
    return dict(
        aware=rate(y, "dx", "plhiv"),
        art_given_aware=rate(y, "art", "dx"),
        vls_given_art=rate(y, "vls", "art"),
        vls_of_plhiv=rate(y, "vls", "plhiv"),
        plhiv_2030=at(y, "plhiv") / nseed,
        cum_infections=cum,
        incidence_2030=inc[2030], incidence_2040=inc[2040])


def main():
    base = {}
    rows = []
    for casc in CASCADES:
        for prep in PREPS:
            d = load(casc, prep)
            if d is None:
                continue
            for lo, hi, lab in AGE_GROUPS:
                for sex in ("f", "m"):
                    m = stratum_metrics(d, sex, lo, hi)
                    key = (sex, lab)
                    if casc == "S0_status_quo" and prep == "P0_none":
                        base[key] = m
                    rows.append(dict(
                        cascade=casc, prep=prep,
                        cascade_name=casc.split("_", 1)[1],
                        prep_name=prep.split("_", 1)[1],
                        description=f"{CASCADES[casc]} {PREPS[prep]}",
                        sex={"f": "Women", "m": "Men"}[sex], age_group=lab,
                        **m))
    t = pd.DataFrame(rows)

    # Everything relative to the status-quo, no-PrEP cell of the SAME stratum.
    for col, newcol in [("cum_infections", "baseline_infections"),
                        ("incidence_2030", "baseline_incidence_2030"),
                        ("incidence_2040", "baseline_incidence_2040")]:
        t[newcol] = [base[(("f" if r.sex == "Women" else "m"), r.age_group)][col]
                     for r in t.itertuples()]
    t["infections_averted"] = t.baseline_infections - t.cum_infections
    t["pct_infections_averted"] = 100 * t.infections_averted / t.baseline_infections
    t["pct_diff_incidence_2030"] = (100 * (t.baseline_incidence_2030
                                           - t.incidence_2030)
                                    / t.baseline_incidence_2030)
    t["pct_diff_incidence_2040"] = (100 * (t.baseline_incidence_2040
                                           - t.incidence_2040)
                                    / t.baseline_incidence_2040)

    order = ["cascade_name", "prep_name", "sex", "age_group", "description",
             "aware", "art_given_aware", "vls_given_art", "vls_of_plhiv",
             "plhiv_2030", "baseline_infections", "cum_infections",
             "infections_averted", "pct_infections_averted",
             "baseline_incidence_2030", "incidence_2030",
             "pct_diff_incidence_2030",
             "baseline_incidence_2040", "incidence_2040",
             "pct_diff_incidence_2040"]
    t = t[order + [c for c in t.columns if c not in order]]
    t.to_csv(OUT / "scenario_table.csv", index=False)
    print(f"wrote outputs/scenario_table.csv: {len(t)} rows")
    print(f"  {t.cascade_name.nunique()} cascades x {t.prep_name.nunique()} "
          f"PrEP x {t.sex.nunique()} sexes x {t.age_group.nunique()} age groups")

    show = t[(t.prep_name == "none") & (t.age_group != "15+ (all)")]
    print(f"\n=== Conditional cascade at {CASC_YEAR}, no PrEP ===")
    print(show[["cascade_name", "sex", "age_group", "aware",
                "art_given_aware", "vls_given_art",
                "vls_of_plhiv"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
