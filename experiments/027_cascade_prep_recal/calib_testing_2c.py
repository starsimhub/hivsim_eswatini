"""Exp 027 step 2c: recalibrate general-population testing under the revised inputs.

Revised instructions (DTC_check/context_dtc.md, 2026-10-28 revision). Fixed
inputs, no longer calibrated:

    ANC testing     0 in 1992 -> 45% of pregnancies 2003 -> 90% 2010, flat after
                    (interventions.ANC_KNOTS_2C)
    PMTCT           prophylaxis for diagnosed pregnant women not on ART,
                    0 in 2003 -> 90% 2010 (pmtct.py)
    VLS among ART   logistic per sex through SHIMS2 (2017.0) and SHIMS3 (2021.0),
                    80% floor, 98% ceiling (art_vls_coverage='phia_logistic')

plus PrEP and shims2_year=2017 as in 2b; row 868's transmission parameters
fixed. Free (fixed testing rates, per goal 5):

    gp_start_m/f    year general-population testing begins, men / women (0 before)
    rel_test_m/f    general-population testing level, men / women
                    (x 0.5 annual probability in 2020)

Loss = SHIMS2 + ART-delivery penalty + Kasturi guide:

  - SHIMS2 first and second 90 at 2017.0, by sex, 15-49 and 15+, sigma from
    the Wilson CI -- the targets, exactly as in 2b (calib_testing.py).
  - The 2b ART-delivery penalty: forced ART (data/art_coverage.csv) can only go
    to diagnosed people, so the share not delivered, pooled over the ART-input
    bins at the worst of 2007.5/2011.5/2014.5, is penalised above 5%.
  - Kasturi guide, 15-49, at 2005.5/2008.5/2011.5/2014.5. Researcher,
    2026-10-28: Kasturi is a guide on the men-vs-women balance of diagnosis
    (no "massive spike" in men from lopsided diagnosis), not a target. So:
      * the men:women ratio of % diagnosed vs Kasturi's, sigma KAS_RATIO_SIGMA
        (0.10) -- the guide proper;
      * each sex's % diagnosed vs Kasturi's level, sigma KAS_LEVEL_SIGMA (0.20)
        -- a loose tie-breaker toward later diagnosis (goal 4.1).
    Why the ratio leads: under the fixed ANC schedule women 15-49 are already
    34-37% diagnosed in 2005 with NO general-population testing for women
    (Kasturi: 7%), so women's level is out of reach, and pulling men down to
    Kasturi's level would widen exactly the gap the researcher objected to.
    Kasturi's 5-year bands are weighted by the run's own PLHIV. Kasturi's ART
    exceeds its own % aware at young ages, so the ART-delivery penalty, not
    Kasturi, sets the floor on early diagnosis.

SHIMS3 (2021) is held out and only reported. Search as 2b: Latin hypercube,
a second hypercube around the best points, best points re-run on fresh seeds.
Outputs: outputs/calib2c_*.csv and calib2c_chosen.json.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
os.chdir(REPO)

from run_sims import make_sim  # noqa: E402
from interventions import ANC_KNOTS_2C  # noqa: E402
from tune_testing import BEST, build_pars, SHIMS  # noqa: E402

OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

BOUNDS = dict(gp_start_m=(1990, 2012), gp_start_f=(1990, 2012),
              rel_test_m=(0.2, 2.5), rel_test_f=(0.05, 2.5))
NAMES = list(BOUNDS)
SHORT_TIMES = (2007.5, 2011.5, 2014.5)
KAS_TIMES = (2005.5, 2008.5, 2011.5, 2014.5)
SNAPS = sorted(set(SHORT_TIMES + KAS_TIMES + (2017.0, 2021.0)))
GROUPS = {"15-49": (15, 50), "15+": (15, 100)}
ART_BINS = [(15, 25), (25, 35), (35, 45), (45, 100)]
KAS_BINS = [(a, a + 5) for a in range(15, 50, 5)]
SHORT_TOL, SHORT_SCALE = 0.05, 0.02
KAS_RATIO_SIGMA, KAS_LEVEL_SIGMA = 0.10, 0.20
N_P1, N_P2, SEEDS_P12, TOP_P3, SEEDS_P3 = 280, 140, 2, 10, 8
NCPUS = 9

# --- SHIMS2 targets: first and second 90 (as calib_testing.py) ----------------
STAGE = {"Diagnosed": "dx", "Diagnosed_OnART": "art_dx", "OnART_VLS": "vls"}
_s = SHIMS[(SHIMS.Year == 2017) & SHIMS["Age group"].isin(GROUPS)].copy()
_s["stage"] = _s.Metric.map(STAGE)
_p, _n, _z = _s.Percentage / 100, _s.Number, 1.96
_s["p"] = _p
_s["sigma"] = _z * np.sqrt(_p * (1 - _p) / _n + _z**2 / (4 * _n**2)) / (1 + _z**2 / _n) / _z
TARGETS = {(r.stage, r.sex, r["Age group"]): (r.p, r.sigma) for _, r in _s.iterrows()}
SCORED = ("dx", "art_dx")

# --- Kasturi % aware by 5-year band (guide) -----------------------------------
_k = pd.read_csv(REPO / "DTC_check" / "SHIMS_Datasets" /
                 "aware_and_art_no_interactions_modeled_output_1995_2021-From_Kasturi.csv",
                 index_col=0)
_k["sex"] = _k.gender.str.lower().map({"female": "f", "male": "m"})
_k["lo"] = _k.agecat.str.split("-").str[0].astype(int)
KAS = {(int(r.year), r.sex, r.lo): r.aware_extrapolated for _, r in _k.iterrows()}


def art_targets(t):
    """Forced ART coverage among PLHIV at time t, per (sex, bin), shims2_year=2017."""
    a = pd.read_csv(REPO / "data" / "art_coverage.csv")
    a.loc[a.Year == 2016, "Year"] = 2017
    a["sex"] = a.Gender.map({0: "f", 1: "m"})
    out = {}
    for (sex, ab), s in a.groupby(["sex", "AgeBin"]):
        lo, hi = (int(x) for x in ab.strip("[)").split(","))
        s = s.sort_values("Year")
        out[(sex, (lo, hi))] = float(np.interp(t, s.Year, s.p_art))
    return out


ART_TARGETS = {t: art_targets(t) for t in SHORT_TIMES}


class Snap(ss.Analyzer):
    """PLHIV / diagnosed / on-ART / suppressed counts at exact times."""

    def init_pre(self, sim):
        super().init_pre(sim)
        self.rows = []

    def step(self):
        y = self.t.now("year")
        hit = [s for s in SNAPS if abs(y - s) < 1e-6]
        if not hit:
            return
        ppl, hiv = self.sim.people, self.sim.diseases.hiv
        groups = {**GROUPS,
                  **{f"art_{lo}_{hi}": (lo, hi) for lo, hi in ART_BINS},
                  **{f"k_{lo}_{hi}": (lo, hi) for lo, hi in KAS_BINS}}
        for sex, sb in (("f", ppl.female), ("m", ppl.male)):
            for name, (lo, hi) in groups.items():
                inb = ppl.alive & sb & (ppl.age >= lo) & (ppl.age < hi)
                self.rows.append(dict(year=hit[0], sex=sex, age=name,
                                      plhiv=(inb & hiv.infected).count(),
                                      dx=(inb & hiv.diagnosed).count(),
                                      art=(inb & hiv.on_art).count(),
                                      vls=(inb & hiv.on_effective_art).count()))


def testing_pars(x):
    """Parameter dict (names as BOUNDS) -> make_sim testing_pars."""
    return dict(gp_start_m=x["gp_start_m"], gp_start_f=x["gp_start_f"],
                rel_test_m=x["rel_test_m"], rel_test_f=x["rel_test_f"],
                anc=True, anc_knots=ANC_KNOTS_2C)


MODEL_KW = dict(prep=True, shims2_year=2017, pmtct=True, art_vls_coverage="phia_logistic")


def run_point(pid, x, seed):
    hiv_pars, network_pars = build_pars(BEST)
    sim = make_sim(seed=seed, start=1985, stop=2022, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars,
                   testing_pars=testing_pars(x), analyzers=[Snap(name="snap")],
                   **MODEL_KW)
    sim.run()
    d = pd.DataFrame(sim.analyzers.snap.rows)
    d["pid"], d["seed"] = pid, seed
    return d


def score(raw):
    """Pool seeds per point -> loss, its parts, ART shortfall, Kasturi comparison."""
    g = raw.groupby(["pid", "year", "sex", "age"])[["plhiv", "dx", "art", "vls"]].sum().reset_index()
    g["p_dx"], g["p_art_dx"], g["p_vls"] = g.dx / g.plhiv, g.art / g.dx, g.vls / g.art
    rows = []
    for pid, gp in g.groupby("pid"):
        r = dict(pid=pid, loss_shims=0.0, loss_kas=0.0)
        at = lambda y, sex, age: gp[(gp.year == y) & (gp.sex == sex) & (gp.age == age)].iloc[0]  # noqa: E731
        for sex in ("f", "m"):
            for grp in GROUPS:
                a17, a21 = at(2017.0, sex, grp), at(2021.0, sex, grp)
                for stage, col in (("dx", "p_dx"), ("art_dx", "p_art_dx"), ("vls", "p_vls")):
                    r[f"{stage}17_{sex}_{grp}"] = a17[col]
                    r[f"{stage}21_{sex}_{grp}"] = a21[col]
                    if stage in SCORED:
                        p, sig = TARGETS[(stage, sex, grp)]
                        r["loss_shims"] += ((a17[col] - p) / sig) ** 2
            # Kasturi level guide: 15-49, PLHIV-weighted by this point's own model
            for t in KAS_TIMES:
                bands = [at(t, sex, f"k_{lo}_{hi}") for lo, hi in KAS_BINS]
                w = np.array([b.plhiv for b in bands], dtype=float)
                kas = float(sum(KAS[(int(t), sex, lo)] * wi
                                for (lo, _), wi in zip(KAS_BINS, w)) / w.sum())
                a = at(t, sex, "15-49")
                r[f"dx{int(t)}_{sex}"], r[f"kas{int(t)}_{sex}"] = a.p_dx, kas
                r[f"artdx{int(t)}_{sex}"] = a.p_art_dx
                r["loss_kas"] += ((a.p_dx - kas) / KAS_LEVEL_SIGMA) ** 2
        # Kasturi balance guide: men:women ratio of % diagnosed
        for t in KAS_TIMES:
            y = int(t)
            ratio = r[f"dx{y}_m"] / r[f"dx{y}_f"]
            r[f"ratio{y}"], r[f"kasratio{y}"] = ratio, r[f"kas{y}_m"] / r[f"kas{y}_f"]
            r["loss_kas"] += ((ratio - r[f"kasratio{y}"]) / KAS_RATIO_SIGMA) ** 2
        short = []
        for t, tg in ART_TARGETS.items():
            want = miss = 0.0
            for (sex, (lo, hi)), target in tg.items():
                row = at(t, sex, f"art_{lo}_{hi}")
                want += target * row.plhiv
                miss += max(0.0, target * row.plhiv - row.art)
            short.append(miss / want)
        r["art_short"] = max(short)
        r["loss_short"] = (max(0.0, r["art_short"] - SHORT_TOL) / SHORT_SCALE) ** 2
        r["loss"] = r["loss_shims"] + r["loss_short"] + r["loss_kas"]
        rows.append(r)
    return pd.DataFrame(rows)


def lhs(n, bounds, seed):
    u = qmc.LatinHypercube(d=len(NAMES), seed=seed).random(n)
    lo = np.array([bounds[k][0] for k in NAMES])
    hi = np.array([bounds[k][1] for k in NAMES])
    return pd.DataFrame(lo + u * (hi - lo), columns=NAMES)


def evaluate(design, seeds, offset):
    design = design.reset_index(drop=True)
    design["pid"] = design.index + offset
    args = [(int(r.pid), {k: r[k] for k in NAMES}, s)
            for _, r in design.iterrows() for s in seeds]
    raw = pd.concat(sc.parallelize(run_point, iterarg=args, ncpus=NCPUS))
    return design, raw, score(raw)


if __name__ == "__main__":
    # Phase 1
    d1, raw1, s1 = evaluate(lhs(N_P1, BOUNDS, 1), range(1, SEEDS_P12 + 1), 0)
    p1 = d1.merge(s1, on="pid"); p1["phase"] = 1

    # Phase 2: hypercube in the box spanned by the best 15 phase-1 points
    top = p1.nsmallest(15, "loss")
    box = {k: (max(BOUNDS[k][0], top[k].min() - 0.1 * np.ptp(top[k])),
               min(BOUNDS[k][1], top[k].max() + 0.1 * np.ptp(top[k]))) for k in NAMES}
    d2, raw2, s2 = evaluate(lhs(N_P2, box, 2), range(1, SEEDS_P12 + 1), 10_000)
    p2 = d2.merge(s2, on="pid"); p2["phase"] = 2
    allp = pd.concat([p1, p2], ignore_index=True)
    allp.to_csv(OUT / "calib2c_search.csv", index=False)

    # Phase 3: best points again on fresh seeds
    best = allp.nsmallest(TOP_P3, "loss")[NAMES]
    d3, raw3, s3 = evaluate(best, range(101, 101 + SEEDS_P3), 20_000)
    p3 = d3.merge(s3, on="pid"); p3["phase"] = 3
    p3.to_csv(OUT / "calib2c_final.csv", index=False)
    chosen = p3.nsmallest(1, "loss").iloc[0]
    sc.savejson(OUT / "calib2c_chosen.json",
                {**{k: float(chosen[k]) for k in NAMES},
                 "loss": float(chosen.loss), "loss_shims": float(chosen.loss_shims),
                 "loss_kas": float(chosen.loss_kas), "art_short": float(chosen.art_short)})

    cols = NAMES + ["loss", "loss_shims", "loss_kas", "loss_short", "art_short"]
    with pd.option_context("display.width", 220, "display.float_format", "{:.3f}".format):
        print("\nPhase 3 (8 fresh seeds), best first:")
        print(p3.sort_values("loss")[cols].to_string(index=False))
        print("\nChosen point vs SHIMS2 (2017.0); 2021 held out:")
        for sex in ("f", "m"):
            for grp in GROUPS:
                line = f"  {sex} {grp:5s}"
                for stage in ("dx", "art_dx", "vls"):
                    line += (f"  {stage} {chosen[f'{stage}17_{sex}_{grp}']:.3f}"
                             f" ({TARGETS[(stage, sex, grp)][0]:.3f})")
                line += f"   | 2021 dx {chosen[f'dx21_{sex}_{grp}']:.3f}"
                print(line)
        print("\nChosen point, 15-49 % diagnosed vs Kasturi (guide), and on ART of diagnosed:")
        for sex in ("f", "m"):
            print(f"  {sex}  " + "  ".join(
                f"{int(t)}: {chosen[f'dx{int(t)}_{sex}']:.2f} (K {chosen[f'kas{int(t)}_{sex}']:.2f},"
                f" art|dx {chosen[f'artdx{int(t)}_{sex}']:.2f})" for t in KAS_TIMES))
        print("  m:f ratio  " + "  ".join(
            f"{int(t)}: {chosen[f'ratio{int(t)}']:.2f} (K {chosen[f'kasratio{int(t)}']:.2f})"
            for t in KAS_TIMES))
