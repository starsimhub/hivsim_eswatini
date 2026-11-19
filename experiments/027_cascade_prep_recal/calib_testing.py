"""Exp 027 step 2b: calibrate testing start years and ramps to SHIMS2.

Step 2 (tune_testing.py) only scaled the testing ramp. Here six testing
inputs are calibrated together, with row 868's transmission parameters fixed
and PrEP + shims2_year=2017 on:

    gp_start        year general-population testing begins (0 before)
    rel_test_m/f    general-population testing level, men / women
                    (x 0.5 annual probability in 2020)
    anc_start       year ANC testing begins
    anc_lo          share of pregnant women tested at anc_start
    anc_dur         years for ANC to ramp from anc_lo to 90%

Targets: SHIMS2 (SHIMS_909090_coverage.csv, compared at 2017.0), by sex, ages
15-49 and 15+:
    first 90   % of PLHIV diagnosed
    second 90  % of diagnosed on ART
sigma = Wilson 95% CI half-width / 1.96 on the file's unweighted n. The third
90 is a model input (vls_construction.py), so testing cannot move it; it is
reported, not scored. SHIMS3 (2021) is held out.

Constraint: ART must be able to reach its forced input (data/art_coverage.csv,
itself SHIMS1/SHIMS2 data). ART only goes to diagnosed people, so testing that
starts too late leaves it short. `art_short` is the share of the forced ART
(summed over the ART-input age/sex bins, in agents) that is not delivered, at
the worst of 2007.5, 2011.5 and 2014.5; loss adds ((art_short - 0.05) /
0.02)^2 above a 5% tolerance. A pooled measure, not a worst-bin one: with
the original testing (diagnosis never limiting) men 15-24 still run up to
~20% short in 2007 because that bin holds a handful of agents and targets
are rounded to whole agents, so a per-bin measure scores noise.

Kasturi's modelled series is NOT a target (researcher, 2026-10-28): it is
reported per point (dx11_*) as a sanity check only.

History: a first version of this script also scored Kasturi's 2011 values;
its outputs are kept as outputs/calib_v1_kasturi_*.

Search: Latin hypercube (phase 1), a second hypercube in the box spanned by
the best phase-1 points (phase 2), then the best points re-run on fresh seeds
(phase 3). Results in outputs/calib_*.csv and calib_chosen.json.
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
from tune_testing import BEST, build_pars, SHIMS  # noqa: E402

OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

BOUNDS = dict(gp_start=(1990, 2012), rel_test_m=(0.2, 2.0), rel_test_f=(0.2, 2.0),
              anc_start=(1990, 2012), anc_lo=(0.1, 0.8), anc_dur=(2, 20))
NAMES = list(BOUNDS)
SHORT_TIMES = (2007.5, 2011.5, 2014.5)
SNAPS = SHORT_TIMES + (2017.0, 2021.0)
GROUPS = {"15-49": (15, 50), "15+": (15, 100)}
ART_BINS = [(15, 25), (25, 35), (35, 45), (45, 100)]
SHORT_TOL, SHORT_SCALE = 0.05, 0.02
N_P1, N_P2, SEEDS_P12, TOP_P3, SEEDS_P3 = 320, 160, 2, 10, 8
NCPUS = 9

# --- SHIMS2 targets: first and second 90 --------------------------------------
STAGE = {"Diagnosed": "dx", "Diagnosed_OnART": "art_dx", "OnART_VLS": "vls"}
_s = SHIMS[(SHIMS.Year == 2017) & SHIMS["Age group"].isin(GROUPS)].copy()
_s["stage"] = _s.Metric.map(STAGE)
_p, _n, _z = _s.Percentage / 100, _s.Number, 1.96
_s["p"] = _p
_s["sigma"] = _z * np.sqrt(_p * (1 - _p) / _n + _z**2 / (4 * _n**2)) / (1 + _z**2 / _n) / _z
TARGETS = {(r.stage, r.sex, r["Age group"]): (r.p, r.sigma) for _, r in _s.iterrows()}
SCORED = ("dx", "art_dx")


def kasturi_2011():
    """Kasturi 2011 % aware, 15-49 and 15+, weighted by model PLHIV (sanity check only)."""
    k = pd.read_csv(REPO / "DTC_check" / "SHIMS_Datasets" /
                    "aware_and_art_no_interactions_modeled_output_1995_2021-From_Kasturi.csv",
                    index_col=0)
    k = k[k.year == 2011].copy()
    k["sex"] = k.gender.str.lower().map({"female": "f", "male": "m"})
    k["lo"] = k.agecat.str.split("-").str[0].astype(int)
    df = sc.loadobj(REPO / "DTC_check" / "outputs" / "027" / "dtc_check_sim.df")
    g = df[df.timevec == 2011]
    out = {}
    for sex in ("f", "m"):
        for grp, (lo, hi) in GROUPS.items():
            kk = k[(k.sex == sex) & (k.lo >= lo) & (k.lo < hi)]
            w = np.array([g[f"popagesex.n_infected_{sex}_{a}_{a + 5}"].mean() for a in kk.lo])
            out[(sex, grp)] = float((kk.aware_extrapolated.values * w).sum() / w.sum())
    return out


KASTURI = kasturi_2011()


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
        groups = {**GROUPS, **{f"art_{lo}_{hi}": (lo, hi) for lo, hi in ART_BINS}}
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
    return dict(gp_start=x["gp_start"], rel_test_m=x["rel_test_m"], rel_test_f=x["rel_test_f"],
                anc=True, anc_start=x["anc_start"], anc_lo=x["anc_lo"],
                anc_end=min(x["anc_start"] + x["anc_dur"], 2030))


def run_point(pid, x, seed):
    hiv_pars, network_pars = build_pars(BEST)
    sim = make_sim(seed=seed, start=1985, stop=2022, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars,
                   testing_pars=testing_pars(x), prep=True, shims2_year=2017,
                   analyzers=[Snap(name="snap")])
    sim.run()
    d = pd.DataFrame(sim.analyzers.snap.rows)
    d["pid"], d["seed"] = pid, seed
    return d


def score(raw):
    """Pool seeds per point -> loss, its parts, ART shortfall, Kasturi check."""
    g = raw.groupby(["pid", "year", "sex", "age"])[["plhiv", "dx", "art", "vls"]].sum().reset_index()
    g["p_dx"], g["p_art_dx"], g["p_vls"] = g.dx / g.plhiv, g.art / g.dx, g.vls / g.art
    g["p_art"] = g.art / g.plhiv
    rows = []
    for pid, gp in g.groupby("pid"):
        r = dict(pid=pid, loss_shims=0.0)
        at = lambda y, sex, age: gp[(gp.year == y) & (gp.sex == sex) & (gp.age == age)].iloc[0]  # noqa: E731
        for sex in ("f", "m"):
            for grp in GROUPS:
                a17, a11, a21 = at(2017.0, sex, grp), at(2011.5, sex, grp), at(2021.0, sex, grp)
                for stage, col in (("dx", "p_dx"), ("art_dx", "p_art_dx"), ("vls", "p_vls")):
                    r[f"{stage}17_{sex}_{grp}"] = a17[col]
                    r[f"{stage}21_{sex}_{grp}"] = a21[col]
                    if stage in SCORED:
                        p, sig = TARGETS[(stage, sex, grp)]
                        r["loss_shims"] += ((a17[col] - p) / sig) ** 2
                r[f"dx11_{sex}_{grp}"] = a11["p_dx"]
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
        r["loss"] = r["loss_shims"] + r["loss_short"]
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
    print("SHIMS2 2017 targets:", {k: round(v[0], 3) for k, v in TARGETS.items()})

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
    allp.to_csv(OUT / "calib_search.csv", index=False)

    # Phase 3: best points again on fresh seeds
    best = allp.nsmallest(TOP_P3, "loss")[NAMES]
    d3, raw3, s3 = evaluate(best, range(101, 101 + SEEDS_P3), 20_000)
    p3 = d3.merge(s3, on="pid"); p3["phase"] = 3
    p3.to_csv(OUT / "calib_final.csv", index=False)
    chosen = p3.nsmallest(1, "loss").iloc[0]
    sc.savejson(OUT / "calib_chosen.json",
                {**{k: float(chosen[k]) for k in NAMES},
                 "loss": float(chosen.loss), "art_short": float(chosen.art_short)})

    cols = NAMES + ["loss", "loss_shims", "loss_short", "art_short"]
    with pd.option_context("display.width", 220, "display.float_format", "{:.3f}".format):
        print("\nPhase 3 (8 fresh seeds), best first:")
        print(p3.sort_values("loss")[cols].to_string(index=False))
        print("\nChosen point vs SHIMS2 (2017.0); 2021 held out; Kasturi 2011 sanity check:")
        for sex in ("f", "m"):
            for grp in GROUPS:
                line = f"  {sex} {grp:5s}"
                for stage in ("dx", "art_dx", "vls"):
                    line += (f"  {stage} {chosen[f'{stage}17_{sex}_{grp}']:.3f}"
                             f" ({TARGETS[(stage, sex, grp)][0]:.3f})")
                line += (f"   | dx 2011.5 {chosen[f'dx11_{sex}_{grp}']:.3f}"
                         f" (Kasturi {KASTURI[(sex, grp)]:.3f})")
                print(line)
        good = allp[allp.loss < allp.loss.min() + 10]
        print(f"\nSpread among search points within 10 loss units of the best (n={len(good)}):")
        print(good[NAMES + [f"dx11_{s}_15-49" for s in 'fm']].describe().loc[["min", "50%", "max"]]
              .to_string())
