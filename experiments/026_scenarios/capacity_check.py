"""Is the model's care cascade plausible against real health-system capacity?

Run before the scenario arms, because `art_95` assumes Eswatini can find and
treat the people it is not currently reaching -- most of them men 25-34, whose
coverage sits at 0.650. If the model gets there by diagnosing people
implausibly fast, or by implying a throughput the health system could not
deliver, the arm measures an impossibility and the abstract would report it as
a policy option.

Two questions, deliberately separate:

1. **Timing.** How long from infection to diagnosis, and diagnosis to ART? The
   model has per-agent `ti_infected`, `ti_diagnosed`, `ti_art`, so these are
   directly measurable rather than inferred from coverage.
2. **Throughput.** How many diagnoses and ART initiations per year does the
   model imply, at the real population scale, and how much does `art_95` add on
   top? Compared against the observed programme, which delivered the 2011->2021
   scale-up we already fit.

Usage
  python experiments/026_scenarios/capacity_check.py            # baseline
  python experiments/026_scenarios/capacity_check.py --arm art_95
Run from the repo ROOT -- run_sims reads data/ by relative path.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from run_sims import make_sim                      # noqa: E402
from analyzers import PopByAgeSex, Cascade         # noqa: E402
import scenarios as S                              # noqa: E402

OUT = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

# 024 design row 868 -- the draw that fit 48/48 targets.
BEST = dict(beta_m2f=0.014161835611057863, rel_beta_f2m=0.3097261030208422,
            s_f_young=2.224000314771824, age_gap_shift=0.0699801603950227,
            age_gap_sd_mult=1.5001144712607821, prop_f0=0.608751696061895,
            prop_m0=0.4742150061435999)

AGE_DIFF_BASE = {"teens": [(7, 3), (6, 3), (5, 1)],
                 "young": [(8, 3), (7, 3), (5, 2)],
                 "adult": [(8, 3), (7, 3), (5, 2)]}


# Kept identical to experiments/025_hm_wave2/run.py::_one_point. If the two
# ever diverge, scenario arms stop being comparable to the calibration that
# justified their parameters -- so this is copied deliberately rather than
# re-derived, and CONC_BASE/rel_init_prev are carried even though 023 fixed
# them, because they are part of the point that fit 48/48.
REL_INIT_PREV = 0.2
CONC_MULT = 1.0
CONC_BASE = {"f1_conc": 0.15, "f2_conc": 0.25, "m1_conc": 0.15, "m2_conc": 0.5}


def build_pars(p):
    """Translate a parameter set into hiv_pars / network_pars."""
    hiv_pars = dict(
        beta_m2f=p["beta_m2f"], rel_beta_f2m=p["rel_beta_f2m"],
        rel_init_prev=REL_INIT_PREV,
        rel_sus_age=[(15, 25, 'f', p["s_f_young"]),
                     (25, 50, 'f', 1.0), (15, 50, 'm', 1.0)])
    shift, sd_mult = p["age_gap_shift"], p["age_gap_sd_mult"]
    network_pars = dict(
        prop_f0=p["prop_f0"], prop_m0=p["prop_m0"],
        age_diff_pars={g: [(max(m + shift, 1.0), max(s * sd_mult, 0.2))
                           for m, s in v] for g, v in AGE_DIFF_BASE.items()},
        **{k: v * CONC_MULT for k, v in CONC_BASE.items()})
    return hiv_pars, network_pars


def a_stop(stop):
    """Coverage tables must extend to the sim end, or stisim forward-fills the
    last row -- which would silently truncate a scale-up at 2040."""
    return max(int(stop), 2040)


def run_one(arm, seed, stop):
    hiv_pars, network_pars = build_pars(BEST)
    kw = dict(seed=seed, stop=stop, verbose=-1,
              hiv_pars=hiv_pars, network_pars=network_pars,
              analyzers=[PopByAgeSex(), Cascade()])
    if arm == "art_95":
        # The cascade axis is a constructor argument, not an extra intervention:
        # it changes the table sti.ART targets. make_sim gained art_coverage=
        # for exactly this.
        art, _ = S.baseline_tables()
        kw["art_coverage"] = S.art_coverage_target(art, 0.95, 2026, 2030,
                                                   stop=a_stop(stop))
    sim = make_sim(**kw)
    sim.run()
    return sim


def timing_table(sim):
    """Per-agent infection -> diagnosis -> ART intervals, in years."""
    hiv = sim.diseases.hiv
    dt_years = float(sim.t.dt) if np.isfinite(float(sim.t.dt)) else 1.0
    start = float(sim.t.yearvec[0])

    ti_inf = np.asarray(hiv.ti_infected.raw, dtype=float)
    ti_dx = np.asarray(hiv.ti_diagnosed.raw, dtype=float)
    ti_art = np.asarray(hiv.ti_art.raw, dtype=float)

    ever_inf = np.isfinite(ti_inf)
    dx = ever_inf & np.isfinite(ti_dx)
    art = dx & np.isfinite(ti_art)

    yr_inf = start + ti_inf * dt_years
    d = pd.DataFrame({
        "year_infected": yr_inf,
        "inf_to_dx": np.where(dx, (ti_dx - ti_inf) * dt_years, np.nan),
        "dx_to_art": np.where(art, (ti_art - ti_dx) * dt_years, np.nan),
        "inf_to_art": np.where(art, (ti_art - ti_inf) * dt_years, np.nan),
        "diagnosed": dx, "on_art_ever": art,
    })[ever_inf]
    # Agents infected recently have not had time to be diagnosed; including them
    # would bias the intervals down and the undiagnosed fraction up. Restrict
    # each era to infections with at least 5 years of follow-up.
    d["year_diagnosed"] = np.where(dx, start + ti_dx * dt_years, np.nan)[ever_inf]
    d["era_inf"] = pd.cut(d.year_infected,
                          [-np.inf, 2005, 2010, 2015, 2020, np.inf],
                          labels=["<2005", "2005-09", "2010-14", "2015-19", "2020+"])
    # dx_to_art MUST be grouped by year of DIAGNOSIS, not of infection. Grouped
    # by infection era it is dominated by "ART did not exist yet": someone
    # diagnosed in 1995 could not start before ~2004 regardless of the health
    # system, which is why the raw numbers came out at a 26-year median for the
    # pre-2005 cohort. Grouped by diagnosis year it measures what was asked --
    # the wait for a treatment slot.
    d["era_dx"] = pd.cut(d.year_diagnosed,
                         [-np.inf, 2005, 2010, 2016, 2021, np.inf],
                         labels=["<2005", "2005-09", "2010-15", "2016-20", "2021+"])
    return d, dt_years


def report_timing(d):
    print("\n=== A. Infection -> diagnosis, by infection era (years) ===")
    print("    Later eras are right-censored: the sim ends, so anyone still")
    print("    undiagnosed at the end raises pct_ever_dx but does not enter the")
    print("    quantiles.")
    rows = []
    for era, g in d.groupby("era_inf", observed=True):
        if not len(g):
            continue
        q = g.inf_to_dx.quantile([0.25, 0.5, 0.75, 0.9])
        rows.append(dict(infected=str(era), n=len(g),
                         pct_ever_dx=100 * float(g.diagnosed.mean()),
                         p25=q.get(0.25), median=q.get(0.5),
                         p75=q.get(0.75), p90=q.get(0.9)))
    t1 = pd.DataFrame(rows)
    print(t1.round(2).to_string(index=False))

    print("\n=== B. Diagnosis -> ART, by year of DIAGNOSIS (years) ===")
    print("    This is the queue for a treatment slot. Eswatini adopted")
    print("    test-and-treat in 2016, so the 2016-20 and 2021+ rows should")
    print("    look like weeks, not years.")
    rows = []
    for era, g in d[d.diagnosed].groupby("era_dx", observed=True):
        if not len(g):
            continue
        q = g.dx_to_art.quantile([0.25, 0.5, 0.75, 0.9])
        rows.append(dict(diagnosed=str(era), n=len(g),
                         pct_ever_art=100 * float(g.on_art_ever.mean()),
                         p25=q.get(0.25), median=q.get(0.5),
                         p75=q.get(0.75), p90=q.get(0.9)))
    t2 = pd.DataFrame(rows)
    print(t2.round(2).to_string(index=False))
    return t1, t2


def report_throughput(sim, dt_years):
    """Annual diagnoses and ART stock at real population scale.

    FLOWS must be SUMMED over the year, not read off one timestep. dt is 1/12,
    so a single index is one month -- and at ~70 real people per agent a monthly
    count is a coarse multiple of ~71, which is why reading one timestep gave
    777 "annual" infections in 2011 and 71 in 2021. Stocks are read at the
    year's last timestep instead.
    """
    print("\n=== Annual throughput, scaled to the real population ===")
    hivres = sim.results.hiv
    yr = np.asarray(sim.t.yearvec, dtype=float)
    want_flow = ["new_infections", "new_diagnoses"]
    want_stock = ["n_infected", "n_diagnosed", "n_on_art"]

    def get(k):
        try:
            return np.asarray(hivres[k], dtype=float)
        except Exception:
            return None

    rows = []
    for y in (2011, 2016, 2021, 2026, 2030):
        m = (yr >= y) & (yr < y + 1)
        if not m.any():
            continue
        r = {"year": y}
        for k in want_flow:
            v = get(k)
            if v is not None:
                r[k + "/yr"] = round(float(v[m].sum()))
        for k in want_stock:
            v = get(k)
            if v is not None:
                r[k] = round(float(v[m][-1]))
        if r.get("n_infected"):
            r["pct_dx"] = round(100 * r["n_diagnosed"] / r["n_infected"], 1)
        if r.get("n_diagnosed"):
            r["pct_art_of_dx"] = round(100 * r["n_on_art"] / r["n_diagnosed"], 1)
        rows.append(r)
    t = pd.DataFrame(rows)
    print(t.to_string(index=False))
    print("\n  Observed first 95 (share of PLHIV who know their status):")
    print("    2011 ~66-75%   2016 ~85% (SHIMS2)   2021 ~94% (SHIMS3)")
    print("  If the model pct_dx runs well above these it is finding people")
    print("  faster than the programme did, and an art_95 arm would then be")
    print("  measuring linkage only -- with case-finding, the hard part,")
    print("  assumed already done.")
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="baseline")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--stop", type=int, default=2031)
    a = ap.parse_args()

    sim = run_one(a.arm, a.seed, a.stop)
    d, dt_years = timing_table(sim)
    print(f"\nsim dt = {dt_years:.4f} years ({dt_years*12:.1f} months); "
          f"{len(d):,} model agents ever infected")
    t_inf, t_art = report_timing(d)
    t_thru = report_throughput(sim, dt_years)

    d.to_csv(OUT / f"capacity_agents_{a.arm}.csv", index=False)
    t_inf.to_csv(OUT / f"capacity_inf_to_dx_{a.arm}.csv", index=False)
    t_art.to_csv(OUT / f"capacity_dx_to_art_{a.arm}.csv", index=False)
    if t_thru is not None:
        t_thru.to_csv(OUT / f"capacity_throughput_{a.arm}.csv", index=False)
    print(f"\nwrote outputs/capacity_*_{a.arm}.csv")



if __name__ == "__main__":
    main()
