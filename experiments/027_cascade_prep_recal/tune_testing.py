"""Exp 027 step 2: tune men's and women's testing so % diagnosed matches SHIMS2.

Row 868's 7 transmission parameters are held fixed. The new ANC ramp,
PrEP, and SHIMS2 inputs pinned to 2017 (shims2_year) are switched on. The general-population testing ramp is scaled by
`rel_test_m` for men and `rel_test_f` for women (interventions.py).

Method
------
1. Grid: one common multiplier for both sexes over GRID, N_SEEDS seeds each.
   Men's and women's testing are separate streams, so each sex's % diagnosed
   depends (almost entirely) on its own multiplier, and one grid serves both.
2. For each sex, interpolate the grid to the multiplier that minimises the
   squared error against SHIMS2 at 2017.0, over ages 15-49 and 15+ together.
3. Confirm: rerun at the chosen (rel_test_m, rel_test_f), N_CONFIRM seeds.

Run from anywhere:  python experiments/027_cascade_prep_recal/tune_testing.py
Outputs go to experiments/027_cascade_prep_recal/outputs/.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
os.chdir(REPO)  # run_sims reads data/ relative to cwd

from run_sims import make_sim  # noqa: E402

OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

# --- Row 868 (copied from DTC_check/build_and_check_calibration.py::BEST) ----
BEST = dict(beta_m2f=0.014161835611057863, rel_beta_f2m=0.3097261030208422,
            s_f_young=2.224000314771824, age_gap_shift=0.0699801603950227,
            age_gap_sd_mult=1.5001144712607821, prop_f0=0.608751696061895,
            prop_m0=0.4742150061435999)
AGE_DIFF_BASE = {"teens": [(7, 3), (6, 3), (5, 1)],
                 "young": [(8, 3), (7, 3), (5, 2)],
                 "adult": [(8, 3), (7, 3), (5, 2)]}
CONC_BASE = {"f1_conc": 0.15, "f2_conc": 0.25, "m1_conc": 0.15, "m2_conc": 0.5}


def build_pars(p):
    """Same mapping as 024/026 and the DTC check script."""
    hiv_pars = dict(beta_m2f=p["beta_m2f"], rel_beta_f2m=p["rel_beta_f2m"],
                    rel_init_prev=0.2,
                    rel_sus_age=[(15, 25, "f", p["s_f_young"]),
                                 (25, 50, "f", 1.0), (15, 50, "m", 1.0)])
    shift, sd = p["age_gap_shift"], p["age_gap_sd_mult"]
    network_pars = dict(prop_f0=p["prop_f0"], prop_m0=p["prop_m0"],
                        age_diff_pars={g: [(max(m + shift, 1.0), max(s * sd, 0.2))
                                           for m, s in v] for g, v in AGE_DIFF_BASE.items()},
                        **CONC_BASE)
    return hiv_pars, network_pars


# --- SHIMS2 first-90 targets (DTC_check/SHIMS_Datasets, labelled 2017) -------
SHIMS = pd.read_csv(REPO / "DTC_check" / "SHIMS_Datasets" / "SHIMS_909090_coverage.csv",
                    encoding="utf-8-sig")
SHIMS["sex"] = SHIMS.Gender.str.lower().map({"female": "f", "male": "m"})
TARGET_AGES = {"15-49": (15, 50), "15+": (15, 100)}
TARGETS = (SHIMS[(SHIMS.Year == 2017) & (SHIMS.Metric == "Diagnosed")
                 & SHIMS["Age group"].isin(TARGET_AGES)]
           .assign(p=lambda d: d.Percentage / 100)
           .set_index(["sex", "Age group"]).p)

SNAP_YEARS = (2011.0, 2017.0, 2021.0)   # 2021 recorded for the held-out check only
ART_BINS = [(15, 25), (25, 35), (35, 45), (45, 100)]  # data/art_coverage.csv bins
GRID = [0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
N_SEEDS, N_CONFIRM = 5, 10


class CascadeSnapshot(ss.Analyzer):
    """Counts at exact time points (not annual means): PLHIV, diagnosed, on ART.

    Rows: (year, sex, age group). Age groups are the SHIMS aggregates plus the
    ART-input bins, so an ART shortfall (forced target > diagnosed pool) shows
    up as on_art/diagnosed ~= 1 in a bin.
    """

    def init_pre(self, sim):
        super().init_pre(sim)
        self.rows = []

    def step(self):
        y = self.t.now("year")
        hit = [s for s in SNAP_YEARS if abs(y - s) < 1e-6]
        if not hit:
            return
        ppl, hiv = self.sim.people, self.sim.diseases.hiv
        groups = dict(TARGET_AGES)
        groups.update({f"{lo}-{hi - 1}" if hi < 100 else f"{lo}+": (lo, hi)
                       for lo, hi in ART_BINS})
        for sex, sb in (("f", ppl.female), ("m", ppl.male)):
            for name, (lo, hi) in groups.items():
                inb = ppl.alive & sb & (ppl.age >= lo) & (ppl.age < hi)
                self.rows.append(dict(
                    year=hit[0], sex=sex, age=name,
                    alive=inb.count(), plhiv=(inb & hiv.infected).count(),
                    dx=(inb & hiv.diagnosed).count(),
                    art=(inb & hiv.on_art).count()))


def run_one(rel_m, rel_f, seed, new_inputs=True):
    hiv_pars, network_pars = build_pars(BEST)
    kw = dict(testing_pars=dict(rel_test_m=rel_m, rel_test_f=rel_f, anc=True),
              prep=True, shims2_year=2017) if new_inputs else {}
    sim = make_sim(seed=seed, start=1985, stop=2026, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars,
                   analyzers=[CascadeSnapshot(name="snap")], **kw)
    sim.run()
    d = pd.DataFrame(sim.analyzers.snap.rows)
    d["rel_m"], d["rel_f"], d["seed"] = rel_m, rel_f, seed
    return d


def summarise(d):
    """Pool counts over seeds -> fractions by (rel_m, rel_f, year, sex, age)."""
    g = d.groupby(["rel_m", "rel_f", "year", "sex", "age"])[["plhiv", "dx", "art"]].sum()
    g["p_dx"] = g.dx / g.plhiv
    g["p_art_dx"] = g.art / g.dx
    return g.reset_index()


def pick(grid_sum):
    """Per sex: multiplier minimising squared error vs SHIMS2 (15-49 and 15+)."""
    fine = np.linspace(min(GRID), max(GRID), 401)
    chosen = {}
    for sex in ("f", "m"):
        sse = np.zeros_like(fine)
        for grp in TARGET_AGES:
            s = grid_sum[(grid_sum.year == 2017.0) & (grid_sum.sex == sex)
                         & (grid_sum.age == grp)].sort_values("rel_m")
            sse += (np.interp(fine, s.rel_m, s.p_dx) - TARGETS[(sex, grp)]) ** 2
        chosen[sex] = float(fine[np.argmin(sse)])
    return chosen


if __name__ == "__main__":
    # 1. Grid
    args = [(r, r, s) for r in GRID for s in range(1, N_SEEDS + 1)]
    grid = pd.concat(sc.parallelize(run_one, iterarg=args, ncpus=8))
    grid.to_csv(OUT / "tune_grid_raw.csv", index=False)
    gs = summarise(grid)
    gs.to_csv(OUT / "tune_grid.csv", index=False)

    # 2. Pick
    chosen = pick(gs)
    print("Chosen multipliers:", {k: round(v, 3) for k, v in chosen.items()})

    # 3. Confirm, alongside the old inputs for reference
    conf = pd.concat(sc.parallelize(
        run_one, iterarg=[(chosen["m"], chosen["f"], s) for s in range(1, N_CONFIRM + 1)],
        ncpus=8))
    old = pd.concat(sc.parallelize(
        run_one, iterarg=[(np.nan, np.nan, s, False) for s in range(1, N_CONFIRM + 1)],
        ncpus=8))
    cs = summarise(conf).assign(inputs="027")
    os_ = summarise(old.fillna({"rel_m": -1, "rel_f": -1})).assign(inputs="row868 as-is")
    both = pd.concat([os_, cs])
    both.to_csv(OUT / "tune_confirm.csv", index=False)
    sc.savejson(OUT / "chosen_testing.json", chosen)

    show = both[both.age.isin(TARGET_AGES)].copy()
    show["shims2"] = [TARGETS.get((s, a)) if y == 2017.0 else np.nan
                      for s, a, y in zip(show.sex, show.age, show.year)]
    with pd.option_context("display.width", 160, "display.float_format", "{:.3f}".format):
        print(show[["inputs", "year", "sex", "age", "p_dx", "p_art_dx", "shims2"]]
              .to_string(index=False))
        print("\nART-input bins at 2017.0 (on ART / diagnosed; ~1.0 flags a shortfall):")
        b = cs[(cs.year == 2017.0) & ~cs.age.isin(TARGET_AGES)]
        print(b[["sex", "age", "p_dx", "p_art_dx"]].to_string(index=False))
