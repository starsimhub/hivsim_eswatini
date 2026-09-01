"""Exp 021 — Viral suppression as an input: removing the vls_coverage = 1.0 default.

Two arms, fixed parameters and seeds:

  A_vls_1.0    vls_coverage=None  -> stisim's default: every ART initiator is
                                    virally suppressed. Current behaviour.
  B_vls_phia   vls_coverage from  -> PHIA sex x year suppression among the
               data/eswatini_vls.csv  treated: 0.913/0.922 (2016), 0.967/0.959
                                    (2021), held flat before 2016.

Outputs
  outputs/sims/{arm}__{seed}.parquet   per-run, written as each finishes
  outputs/results.parquet              concatenated
  outputs/cascade_vs_phia.csv          metric 1: population VLS vs SHIMS2/3
  outputs/ab.csv                       metrics 2-4: arm B vs arm A with z
  outputs/prevalence.csv               metric 2: by age/sex vs PHIA
"""

import os
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                  NUMEXPR_NUM_THREADS="1", MKL_NUM_THREADS="1")

import argparse
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sciris as sc

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
OUT_DIR, FIG_DIR = EXP_DIR / "outputs", EXP_DIR / "figures"
SIM_DIR = OUT_DIR / "sims"
for d in (OUT_DIR, FIG_DIR, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
from run_sims import make_sim  # noqa: E402
from analyzers import PopByAgeSex, Cascade  # noqa: E402
import vls_construction as vc  # noqa: E402

N_AGENTS = 10_000
STOP = 2026
N_SEEDS = 10
PARS = {"beta_m2f": 0.0139, "rel_init_prev": 0.49}

# Held flat back to this year -- an assumption, not a measurement. The surveys
# start in 2016; ART starts around 2004. Early-ART-era suppression was
# plausibly worse than 0.92, so holding flat understates the effect of removing
# the 1.0 default rather than overstating it. See README.
VLS_FILL_BACK_TO = 1985

PHIA_SEX = {0: "m", 1: "f"}   # prevalence_by_age_sex.csv -- NB opposite of stisim's

# PHIA population VLS among all PLHIV, ages 15-49, for metric 1. From
# data/eswatini_vls.csv (SHIMS2 Table 9.4.A / SHIMS3 Table 8.2 aggregate rows).
PHIA_POP_VLS = {(2016, "m"): 0.623, (2016, "f"): 0.748,
                (2021, "m"): 0.824, (2021, "f"): 0.886}

KEEP = ("timevec", "hiv.new_deaths", "hiv.prevalence_15_49", "hiv.n_infected",
        "hiv.new_infections", "hiv.p_on_art", "n_alive", "popagesex.",
        "cascade.")


def vls_table():
    """The PHIA input table, in the shape sti.ART(vls_coverage=...) expects."""
    df = vc.build()
    tbl = vc.to_vls_coverage(df, fill_back_to=VLS_FILL_BACK_TO)
    # vls_coverage silently defaults any stratum it is NOT given to 1.0, which
    # would reintroduce the default this experiment removes. Assert coverage of
    # both sexes over the full adult range before anything runs.
    assert set(tbl.Gender) == {"m", "f"}, f"both sexes required, got {set(tbl.Gender)}"
    assert (tbl.AgeBin == "[15,100)").all(), "input must span the adult range"
    assert tbl.p_vls.between(0, 1).all(), "p_vls must be a proportion"
    assert tbl.Year.min() <= VLS_FILL_BACK_TO, "series must reach the model start"
    return tbl


ARMS = {"A_vls_1.0": None, "B_vls_phia": "phia"}


def _run(arm, seed):
    path = SIM_DIR / f"{arm}__{seed:03d}.parquet"
    if path.exists():
        return

    vls = vls_table() if ARMS[arm] == "phia" else None
    t0 = time.perf_counter()
    sim = make_sim(seed=seed, stop=STOP, verbose=-1, hiv_pars=dict(PARS),
                   art_vls_coverage=vls,
                   analyzers=[PopByAgeSex(), Cascade()])
    sim.pars.n_agents = N_AGENTS
    sim.run()
    elapsed = time.perf_counter() - t0

    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns if any(c == k or c.startswith(k) for k in KEEP)]
    out = df[keep].copy()
    out["arm"], out["seed"], out["runtime_s"] = arm, seed, elapsed
    out["pop_scale"] = float(out["popagesex.n_alive_total"].iloc[0]) / N_AGENTS
    out.to_parquet(path, index=False)


def load():
    files = sorted(SIM_DIR.glob("*.parquet"))
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


# --- Metric 1: population VLS vs PHIA ---------------------------------------

def cascade_vs_phia(df):
    """Model VLS among all PLHIV vs the PHIA survey figures.

    Not circular in the fitting sense: ART coverage and suppression-given-ART
    are both inputs, so suppression among *all* PLHIV is an emergent output of
    the cascade -- diagnosis, initiation, retention -- and a genuine check.
    """
    rows = []
    for (arm, year), g in df[df.timevec.isin([2016, 2021])].groupby(["arm", "timevec"]):
        for sex in ("m", "f"):
            for stem, label in (("p_vls", "vls_among_plhiv"),
                                ("p_vls_given_art", "vls_given_art"),
                                ("p_on_art", "art_coverage")):
                col = f"cascade.{stem}_{sex}_15_49"
                if col not in df.columns:
                    continue
                rows.append(dict(
                    arm=arm, year=int(year), sex=sex, quantity=label,
                    model=g[col].mean(), model_sd=g[col].std(ddof=1),
                    phia=PHIA_POP_VLS.get((int(year), sex)) if label == "vls_among_plhiv" else np.nan))
    out = pd.DataFrame(rows)
    out["diff"] = out.model - out.phia
    return out


# --- Metrics 2-4: the A/B ----------------------------------------------------

def ab_table(df):
    a, b = df[df.arm == "A_vls_1.0"], df[df.arm == "B_vls_phia"]
    metrics = ["hiv.prevalence_15_49", "hiv.new_infections", "hiv.new_deaths"]
    rows = []
    for yr in (2005, 2011, 2016, 2021, 2025):
        for m in metrics:
            if m not in df.columns:
                continue
            x = a[a.timevec == yr][m]
            y = b[b.timevec == yr][m]
            if not len(x) or not len(y):
                continue
            se = np.sqrt(x.std(ddof=1)**2 / len(x) + y.std(ddof=1)**2 / len(y))
            rows.append(dict(year=yr, metric=m, arm_A=x.mean(), arm_B=y.mean(),
                             rel=(y.mean() - x.mean()) / x.mean() if x.mean() else np.nan,
                             z=(y.mean() - x.mean()) / se if se else np.nan))
    return pd.DataFrame(rows)


def prevalence_table(df):
    tg = pd.read_csv(REPO_ROOT / "calibration_data" / "prevalence_by_age_sex.csv")
    rows = []
    for arm, g_all in df.groupby("arm"):
        for _, t in tg.iterrows():
            yr, lo = int(t.Year), int(t["start age"])
            sex = PHIA_SEX[int(t.Gender)]
            acol = f"popagesex.n_alive_{sex}_{lo}_{lo + 5}"
            icol = f"popagesex.n_infected_{sex}_{lo}_{lo + 5}"
            if acol not in df.columns:
                continue
            g = g_all[g_all.timevec == yr]
            if not len(g):
                continue
            rows.append(dict(arm=arm, year=yr, sex=sex, age_low=lo,
                             model_prev=(g[icol] / g[acol].replace(0, np.nan)).mean(),
                             phia_prev=t.NationalPrevalence,
                             phia_lb=t.lb, phia_ub=t.ub))
    return pd.DataFrame(rows)


# --- Figures -----------------------------------------------------------------

COLORS = {"A_vls_1.0": "#4c78a8", "B_vls_phia": "#e45756"}


def plot_main(df, cvp, path):
    fig, axes = plt.subplots(1, 4, figsize=(21, 4.4))
    panels = [("hiv.prevalence_15_49", "HIV prevalence 15-49"),
              ("hiv.new_infections", "New infections"),
              ("hiv.new_deaths", "AIDS deaths")]
    for ax, (col, title) in zip(axes, panels):
        for arm, g in df.groupby("arm"):
            m = g.groupby("timevec")[col].agg(["mean", "std"])
            ax.plot(m.index, m["mean"], label=arm, color=COLORS.get(arm))
            ax.fill_between(m.index, m["mean"] - m["std"], m["mean"] + m["std"],
                            alpha=0.15, color=COLORS.get(arm))
        ax.set(title=title, xlabel="year")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    ax = axes[3]
    for arm, g in df.groupby("arm"):
        for sex, ls in (("f", "-"), ("m", "--")):
            col = f"cascade.p_vls_{sex}_15_49"
            if col not in df.columns:
                continue
            m = g.groupby("timevec")[col].mean()
            ax.plot(m.index, m, ls, color=COLORS.get(arm),
                    label=f"{arm} {sex}", alpha=0.85)
    for (yr, sex), v in PHIA_POP_VLS.items():
        ax.plot(yr, v, "ko", ms=6)
        ax.annotate(sex, (yr, v), fontsize=7, xytext=(3, 3),
                    textcoords="offset points")
    ax.set(title="Population VLS among PLHIV 15-49\n(black = PHIA)",
           xlabel="year", xlim=(2004, 2026))
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.suptitle("Exp 021 — viral suppression as an input", y=1.03)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_prevalence(pt, path):
    years = sorted(pt.year.unique())
    fig, axes = plt.subplots(2, len(years), figsize=(4.6 * len(years), 7.2),
                             squeeze=False)
    for j, yr in enumerate(years):
        for i, sex in enumerate(("f", "m")):
            ax = axes[i][j]
            sub = pt[(pt.year == yr) & (pt.sex == sex)]
            for arm, g in sub.groupby("arm"):
                g = g.sort_values("age_low")
                ax.plot(g.age_low, g.model_prev, label=arm, color=COLORS.get(arm))
            t = sub.drop_duplicates(["age_low"]).sort_values("age_low")
            ax.errorbar(t.age_low, t.phia_prev,
                        yerr=[t.phia_prev - t.phia_lb, t.phia_ub - t.phia_prev],
                        fmt="ko", ms=4, capsize=2, label="PHIA")
            ax.set(title=f"{'Female' if sex == 'f' else 'Male'} {yr}", xlabel="age")
            if j == 0:
                ax.set_ylabel("HIV prevalence")
            ax.grid(alpha=0.3)
            if i == 0 and j == 0:
                ax.legend(fontsize=7)
    fig.suptitle("Exp 021 — prevalence by age and sex vs PHIA", y=1.0)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_seeds", type=int, default=N_SEEDS)
    p.add_argument("--n_workers", type=int, default=None)
    p.add_argument("--plot_only", action="store_true")
    args = p.parse_args()

    print("VLS input table (arm B):")
    print(vls_table().to_string(index=False))

    if not args.plot_only:
        work = [(arm, s) for arm in ARMS for s in range(args.n_seeds)]
        todo = [w for w in work
                if not (SIM_DIR / f"{w[0]}__{w[1]:03d}.parquet").exists()]
        print(f"\n{len(work) - len(todo)}/{len(work)} on disk; running {len(todo)}")
        if todo:
            t0 = sc.tic()
            sc.parallelize(_run, iterarg=todo, ncpus=args.n_workers)
            sc.toc(t0, label="021")

    df = load()
    if not len(df):
        print("no results")
        return
    df.to_parquet(OUT_DIR / "results.parquet", index=False)

    print("\n=== Metric 1: cascade vs PHIA (ages 15-49) ===")
    cvp = cascade_vs_phia(df)
    cvp.to_csv(OUT_DIR / "cascade_vs_phia.csv", index=False)
    print(cvp.round(4).to_string(index=False))

    print("\n=== Metrics 2-4: arm B (PHIA VLS) vs arm A (default 1.0) ===")
    ab = ab_table(df)
    ab.to_csv(OUT_DIR / "ab.csv", index=False)
    print(ab.round(4).to_string(index=False))

    pt = prevalence_table(df)
    if len(pt):
        pt.to_csv(OUT_DIR / "prevalence.csv", index=False)
        print("\n=== Metric 2: prevalence fit vs PHIA ===")
        print(pt.assign(err=pt.model_prev - pt.phia_prev).groupby("arm").agg(
            mean_abs_err=("err", lambda s: s.abs().mean()),
            mean_bias=("err", "mean")).round(4).to_string())
        plot_prevalence(pt, FIG_DIR / "prevalence_vs_phia.png")

    plot_main(df, cvp, FIG_DIR / "vls_effect.png")
    print(f"\nfigures -> {FIG_DIR}")


if __name__ == "__main__":
    main()
