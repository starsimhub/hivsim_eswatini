"""Exp 019 — Age-dependent untreated survival: is that the missing AIDS mortality?

Four arms at fixed parameters and seeds, varying only a multiplier on
`dur_latent` as a function of age at infection (see ../../hiv_survival.py):

  A_flat_13     1.0 everywhere              13.1 y   baseline (= upstream exactly)
  B_flat_11.5   0.84 everywhere            ~11.5 y   LEVEL effect (A->B)
  C_grad_mild   0.89/0.84/0.73/0.61        ~11.5 y   gradient, half strength
  D_grad_alpha  0.94/0.84/0.64/0.44        ~11.5 y   GRADIENT effect (B->D), full ALPHA

B, C and D sit at the same level and differ only in shape, so B->D isolates the
gradient at constant level and A->B isolates level at zero gradient.

Metrics, in the order the README sets out (metric 1 gates the rest):
  1 deaths by route -- ti_zero vs p_hiv_death. 017 estimated ~80/20 by
    integrating the hazard table; this measures it.
  2 AIDS deaths vs UNAIDS -- trajectory and peak-to-peak. The target.
  3 realized survival per arm -- validates the subclass and quantifies the
    level/gradient separation.
  4 age distribution of AIDS deaths vs 016's implied deaths by age/sex.
  5 prevalence 15-49 and by age/sex vs PHIA -- the cost side of the trade-off.
  6 population, as a demographic sanity check.

Outputs
  outputs/sims/{arm}__{seed}.parquet          per-run results, written as each finishes
  outputs/sims/{arm}__{seed}__deaths.parquet  per-run agent-level death records
  outputs/results.parquet, deaths.parquet     concatenated
  outputs/routes.csv          metric 1: deaths by route, per arm
  outputs/deaths_vs_unaids.csv metric 2: peak and trajectory vs UNAIDS
  outputs/survival.csv        metric 3: realized survival by arm and age band
  outputs/deaths_by_age.csv   metric 4: age distribution vs 016
  outputs/prevalence.csv      metric 5: prevalence vs PHIA
"""

import os
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                  NUMEXPR_NUM_THREADS="1", MKL_NUM_THREADS="1")

import argparse
import sys
import time
from functools import partial
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sciris as sc

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
EXP_016 = REPO_ROOT / "experiments" / "016_double_counted_mortality"
OUT_DIR, FIG_DIR = EXP_DIR / "outputs", EXP_DIR / "figures"
SIM_DIR = OUT_DIR / "sims"
for d in (OUT_DIR, FIG_DIR, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
from run_sims import make_sim  # noqa: E402
from analyzers import PopByAgeSex  # noqa: E402
from hiv_survival import AgeDependentSurvival, ARMS  # noqa: E402

N_AGENTS = 10_000
STOP = 2026

# High-transmission set only. 017 established that the default-parameter arm
# cannot carry quantitative claims -- though see 018 obs 4, which found that was
# a property of the 1.5.8 stack (CV 45.9%) and not of model-v1.1 (7.8%). This
# experiment is entirely about effect sizes, so the lower-variance point is
# still the right one.
PARS = {"beta_m2f": 0.0139, "rel_init_prev": 0.49}

# Seeds: 10, from 018 obs 5 -- CV is 4.4% at this parameter point on
# model-v1.1, which sits at the 3-5 boundary of the standard bands. 10 is
# comfortable.
N_SEEDS = 10

AGE_BANDS = [(15, 25), (25, 35), (35, 45), (45, 200)]
# Restrict survival statistics to infections early enough to have died within
# the window. Without this, mean survival is biased downward by right-censoring
# -- the smoke test showed 10.0 y for an arm whose true mean is 13.1 y, purely
# because only the fast progressors had died by the cutoff.
SURV_INFECTION_CUTOFF = 2000

KEEP = ("timevec", "hiv.new_deaths", "hiv.new_deaths_progression",
        "hiv.new_deaths_hazard", "hiv.prevalence_15_49", "hiv.n_infected",
        "hiv.p_on_art", "n_alive", "popagesex.")

# PHIA Gender coding: 0 = Male, 1 = Female. Per experiments/008 (the
# authoritative ingestion) and confirmed against the data itself -- Gender=1
# gives 0.101 at 15-19 in 2007, which is SDHS 2006-07's figure for women.
# NB experiments/018/run.py:144 has this inverted; see this experiment's SUMMARY.
PHIA_SEX = {0: "m", 1: "f"}


# --- Runs --------------------------------------------------------------------

def _run(arm, seed):
    """One sim; writes before returning so an interrupted batch is recoverable."""
    path = SIM_DIR / f"{arm}__{seed:03d}.parquet"
    dpath = SIM_DIR / f"{arm}__{seed:03d}__deaths.parquet"
    if path.exists() and dpath.exists():
        return

    t0 = time.perf_counter()
    hiv_class = partial(AgeDependentSurvival, latent_mult=ARMS[arm])
    sim = make_sim(seed=seed, stop=STOP, verbose=-1, hiv_pars=dict(PARS),
                   hiv_class=hiv_class, analyzers=[PopByAgeSex()])
    sim.pars.n_agents = N_AGENTS
    sim.run()
    elapsed = time.perf_counter() - t0

    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns if any(c == k or c.startswith(k) for k in KEEP)]
    out = df[keep].copy()

    # pop_scale, needed to put agent-level death records on the same footing as
    # the scale=True results arrays. Derived rather than read off pars so it
    # cannot silently disagree with what the results were actually scaled by.
    scale = float(out["popagesex.n_alive_total"].iloc[0]) / N_AGENTS
    out["arm"], out["seed"] = arm, seed
    out["runtime_s"], out["pop_scale"] = elapsed, scale
    out.to_parquet(path, index=False)

    deaths = sim.diseases.hiv.death_records()
    deaths["arm"], deaths["seed"], deaths["pop_scale"] = arm, seed, scale
    deaths["year"] = sim.t.yearvec[deaths["ti"].astype(int)] if len(deaths) else []
    deaths["age_at_death"] = deaths["age_at_infection"] + deaths["survival_years"]
    deaths["infection_year"] = deaths["year"] - deaths["survival_years"]
    deaths.to_parquet(dpath, index=False)


def load(pattern="*.parquet", exclude_deaths=True):
    files = sorted(SIM_DIR.glob(pattern))
    if exclude_deaths:
        files = [f for f in files if not f.name.endswith("__deaths.parquet")]
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


def load_deaths():
    files = sorted(SIM_DIR.glob("*__deaths.parquet"))
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


# --- Metric 1: deaths by route ----------------------------------------------

def routes_table(df):
    rows = []
    for arm, g in df.groupby("arm"):
        per_seed = g.groupby("seed")[["hiv.new_deaths_progression",
                                     "hiv.new_deaths_hazard"]].sum()
        prog, haz = per_seed.mean()
        tot = prog + haz
        rows.append(dict(arm=arm, progression=prog, hazard=haz, total=tot,
                         progression_share=prog / tot if tot else np.nan,
                         check_vs_new_deaths=g.groupby("seed")["hiv.new_deaths"]
                                              .sum().mean()))
    return pd.DataFrame(rows).sort_values("arm")


# --- Metric 2: AIDS deaths vs UNAIDS ----------------------------------------

def unaids_deaths():
    d = pd.read_csv(REPO_ROOT / "data" / "eswatini_hiv_calib.csv")
    return d[["time", "hiv.new_deaths"]].dropna().rename(
        columns={"time": "year", "hiv.new_deaths": "unaids_deaths"})


def deaths_vs_unaids(df, un):
    rows = []
    un_peak = un.loc[un.unaids_deaths.idxmax()]
    for arm, g in df.groupby("arm"):
        m = g.groupby("timevec")["hiv.new_deaths"].mean()
        pk_year = m.idxmax()
        rows.append(dict(
            arm=arm, model_peak=m.max(), model_peak_year=pk_year,
            unaids_peak=un_peak.unaids_deaths, unaids_peak_year=un_peak.year,
            pct_of_unaids_peak=100 * m.max() / un_peak.unaids_deaths,
            deficit_at_peak=un_peak.unaids_deaths - m.max(),
        ))
    return pd.DataFrame(rows).sort_values("arm")


# --- Metric 3: realized survival --------------------------------------------

def survival_table(dd):
    """Untreated survival by arm and age band, on an uncensored subset."""
    u = dd[(~dd.on_art) & (dd.age_at_infection >= 15)
           & (dd.infection_year <= SURV_INFECTION_CUTOFF)]
    rows = []
    for arm, g in u.groupby("arm"):
        rows.append(dict(arm=arm, band="all", n=len(g),
                         mean_years=g.survival_years.mean(),
                         median_years=g.survival_years.median()))
        for lo, hi in AGE_BANDS:
            b = g[(g.age_at_infection >= lo) & (g.age_at_infection < hi)]
            rows.append(dict(arm=arm, band=f"{lo}-{hi if hi < 200 else '+'}",
                             n=len(b),
                             mean_years=b.survival_years.mean() if len(b) else np.nan,
                             median_years=b.survival_years.median() if len(b) else np.nan))
    return pd.DataFrame(rows)


# --- Metric 4: age distribution of AIDS deaths ------------------------------

def deaths_by_age(dd):
    """Model AIDS deaths by 5-year age band and sex, scaled to population.

    Compared against 016's implied AIDS deaths, which were reconstructed from
    all-cause mortality with no HIV information as input -- an independent
    age-stratified target the model has never been checked against.
    """
    tgt_path = EXP_016 / "outputs" / "excess_deaths_by_age_sex.csv"
    years = [1995, 2000, 2005, 2010, 2015, 2020]
    dd = dd.copy()
    dd["age_bin"] = (dd.age_at_death // 5 * 5).clip(upper=80)
    # floor, not round: rounding straddles the calendar-year boundary, pulling in
    # the back half of the previous year and dropping the back half of this one.
    # It inflated late-epidemic counts by ~35% and made the model look like it
    # over-killed in 2020. With floor, these counts reconcile with the
    # hiv.new_deaths trajectory to within a few percent.
    dd["year_r"] = np.floor(dd.year).astype(int)

    rows = []
    for (arm, seed), g in dd.groupby(["arm", "seed"]):
        scale = g.pop_scale.iloc[0]
        for yr in years:
            sub = g[g.year_r == yr]
            for female, sg in sub.groupby("female"):
                counts = sg.groupby("age_bin").size() * scale
                for age, n in counts.items():
                    rows.append(dict(arm=arm, seed=seed, year=yr,
                                     sex="Female" if female else "Male",
                                     age=int(age), model_deaths=n))
    model = pd.DataFrame(rows)
    if not len(model):
        return model

    # Fill absent (seed, year, sex, age) cells with zero before averaging over
    # seeds. Without this, a cell with a death in 3 of 10 seeds averages over
    # those 3 only, which inflates every sparse bin -- badly enough in the
    # late epidemic to invent an apparent excess of deaths in 2020.
    keys = ["arm", "seed", "year", "sex", "age"]
    grid = pd.MultiIndex.from_product(
        [sorted(model.arm.unique()), sorted(model.seed.unique()), years,
         ["Female", "Male"], sorted(model.age.unique())], names=keys)
    model = (model.set_index(keys)["model_deaths"]
                  .reindex(grid, fill_value=0.0).reset_index())
    model = (model.groupby(["arm", "year", "sex", "age"])["model_deaths"]
                  .mean().reset_index())

    if tgt_path.exists():
        tgt = pd.read_csv(tgt_path)
        tgt = tgt.melt(id_vars=["sex", "age"], var_name="year",
                       value_name="implied_deaths_016")
        tgt["year"] = tgt.year.astype(int)
        model = model.merge(tgt, on=["sex", "age", "year"], how="left")
    return model


# --- Metric 5: prevalence vs PHIA -------------------------------------------

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
            per_seed = g[icol] / g[acol].replace(0, np.nan)
            rows.append(dict(
                arm=arm, year=yr, sex=sex, age_low=lo,
                model_prev=per_seed.mean(), model_sd=per_seed.std(ddof=1),
                phia_prev=t.NationalPrevalence, phia_lb=t.lb, phia_ub=t.ub,
                # expected infected AGENTS, for the rare-event floor
                expected_agents=(g[icol] / g["pop_scale"]).mean(),
            ))
    return pd.DataFrame(rows)


# --- Figures -----------------------------------------------------------------

ARM_COLORS = {"A_flat_13": "#4c78a8", "B_flat_11.5": "#f58518",
              "C_grad_mild": "#54a24b", "D_grad_alpha": "#e45756"}


def plot_deaths(df, un, path):
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
    for arm, g in df.groupby("arm"):
        m = g.groupby("timevec")["hiv.new_deaths"].agg(["mean", "std"])
        c = ARM_COLORS.get(arm)
        axes[0].plot(m.index, m["mean"], label=arm, color=c)
        axes[0].fill_between(m.index, m["mean"] - m["std"], m["mean"] + m["std"],
                             alpha=0.15, color=c)
    axes[0].plot(un.year, un.unaids_deaths, "k--", lw=2, label="UNAIDS")
    axes[0].set(title="AIDS deaths vs UNAIDS", xlabel="year", ylabel="deaths/year")
    axes[0].legend(fontsize=8)

    r = routes_table(df).set_index("arm")
    bottom = np.zeros(len(r))
    for col, lbl, c in (("progression", "ti_zero (untunable)", "#e45756"),
                        ("hazard", "p_hiv_death (tunable)", "#4c78a8")):
        axes[1].bar(r.index, r[col], bottom=bottom, label=lbl, color=c)
        bottom += r[col].values
    axes[1].set(title="Cumulative deaths by route", ylabel="deaths")
    axes[1].tick_params(axis="x", rotation=20)
    axes[1].legend(fontsize=8)

    for arm, g in df.groupby("arm"):
        m = g.groupby("timevec")["hiv.prevalence_15_49"].mean()
        axes[2].plot(m.index, m, label=arm, color=ARM_COLORS.get(arm))
    axes[2].set(title="HIV prevalence 15-49 (the trade-off)", xlabel="year")
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle("Exp 019 — age-dependent untreated survival", y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_survival(st, path):
    bands = [b for b in st.band.unique() if b != "all"]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    width = 0.2
    x = np.arange(len(bands))
    for i, (arm, g) in enumerate(st.groupby("arm")):
        vals = [g[g.band == b].mean_years.iloc[0] if len(g[g.band == b]) else np.nan
                for b in bands]
        ax.bar(x + i * width, vals, width, label=arm, color=ARM_COLORS.get(arm))
    ax.set(xticks=x + 1.5 * width, xlabel="age at infection",
           ylabel="mean untreated survival (years)",
           title=f"Realized survival, infections up to {SURV_INFECTION_CUTOFF} "
                 f"(uncensored subset)")
    ax.set_xticklabels(bands)
    ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_deaths_by_age(ma, path):
    if not len(ma) or "implied_deaths_016" not in ma.columns:
        return
    years = [y for y in (2000, 2005, 2010) if y in set(ma.year)]
    if not years:
        return
    fig, axes = plt.subplots(2, len(years), figsize=(5 * len(years), 7.5),
                             squeeze=False)
    for j, yr in enumerate(years):
        for i, sex in enumerate(("Female", "Male")):
            ax = axes[i][j]
            sub = ma[(ma.year == yr) & (ma.sex == sex)]
            for arm, g in sub.groupby("arm"):
                ax.plot(g.age, g.model_deaths, label=arm, color=ARM_COLORS.get(arm))
            t = sub.drop_duplicates(["age"])
            ax.plot(t.age, t.implied_deaths_016, "k--", lw=2, label="016 implied")
            ax.set(title=f"{sex} {yr}", xlabel="age")
            if j == 0:
                ax.set_ylabel("AIDS deaths")
            ax.grid(alpha=0.3)
            if i == 0 and j == 0:
                ax.legend(fontsize=7)
    fig.suptitle("Exp 019 metric 4 — age distribution of AIDS deaths vs 016", y=1.0)
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
                ax.plot(g.age_low, g.model_prev, label=arm,
                        color=ARM_COLORS.get(arm))
            t = sub.drop_duplicates(["age_low"]).sort_values("age_low")
            ax.errorbar(t.age_low, t.phia_prev,
                        yerr=[t.phia_prev - t.phia_lb, t.phia_ub - t.phia_prev],
                        fmt="ko", ms=4, capsize=2, label="PHIA")
            ax.set(title=f"{'Female' if sex == 'f' else 'Male'} {yr}",
                   xlabel="age")
            if j == 0:
                ax.set_ylabel("HIV prevalence")
            ax.grid(alpha=0.3)
            if i == 0 and j == 0:
                ax.legend(fontsize=7)
    fig.suptitle("Exp 019 metric 5 — prevalence by age and sex vs PHIA", y=1.0)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


# --- Main --------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_seeds", type=int, default=N_SEEDS)
    p.add_argument("--n_workers", type=int, default=None)
    p.add_argument("--arms", nargs="*", default=list(ARMS))
    p.add_argument("--plot_only", action="store_true")
    args = p.parse_args()

    if not args.plot_only:
        work = [(arm, s) for arm in args.arms for s in range(args.n_seeds)]
        todo = [w for w in work
                if not (SIM_DIR / f"{w[0]}__{w[1]:03d}__deaths.parquet").exists()]
        print(f"{len(work) - len(todo)}/{len(work)} on disk; running {len(todo)} "
              f"(n_workers={args.n_workers or 'all'})")
        if todo:
            t0 = sc.tic()
            sc.parallelize(_run, iterarg=todo, ncpus=args.n_workers)
            sc.toc(t0, label="019")

    df, dd = load(), load_deaths()
    if not len(df):
        print("no results on disk")
        return
    df.to_parquet(OUT_DIR / "results.parquet", index=False)
    if len(dd):
        dd.to_parquet(OUT_DIR / "deaths.parquet", index=False)

    un = unaids_deaths()

    print("\n=== Metric 1: deaths by route (mean per seed, whole run) ===")
    rt = routes_table(df)
    rt.to_csv(OUT_DIR / "routes.csv", index=False)
    print(rt.to_string(index=False))

    print("\n=== Metric 2: AIDS deaths vs UNAIDS ===")
    dv = deaths_vs_unaids(df, un)
    dv.to_csv(OUT_DIR / "deaths_vs_unaids.csv", index=False)
    print(dv.to_string(index=False))

    print(f"\n=== Metric 3: realized untreated survival "
          f"(infections <= {SURV_INFECTION_CUTOFF}) ===")
    st = survival_table(dd) if len(dd) else pd.DataFrame()
    if len(st):
        st.to_csv(OUT_DIR / "survival.csv", index=False)
        print(st.to_string(index=False))

    print("\n=== Metric 4: age distribution of AIDS deaths vs 016 ===")
    ma = deaths_by_age(dd) if len(dd) else pd.DataFrame()
    if len(ma):
        ma.to_csv(OUT_DIR / "deaths_by_age.csv", index=False)
        if "implied_deaths_016" in ma:
            agg = ma.dropna(subset=["implied_deaths_016"]).groupby(["arm", "year"])[
                ["model_deaths", "implied_deaths_016"]].sum()
            agg["ratio"] = agg.model_deaths / agg.implied_deaths_016
            print(agg.round(1).to_string())

    print("\n=== Metric 5: prevalence vs PHIA ===")
    pt = prevalence_table(df)
    if len(pt):
        pt.to_csv(OUT_DIR / "prevalence.csv", index=False)
        summ = pt.assign(err=pt.model_prev - pt.phia_prev).groupby("arm").agg(
            mean_abs_err=("err", lambda s: s.abs().mean()),
            mean_bias=("err", "mean"),
            min_expected_agents=("expected_agents", "min"),
            n_strata_under_10=("expected_agents", lambda s: int((s < 10).sum())),
        )
        print(summ.round(4).to_string())

    print("\n=== Metric 6: population at 2026 ===")
    print(df[df.timevec == df.timevec.max()].groupby("arm")["n_alive"]
          .mean().round(0).to_string())

    plot_deaths(df, un, FIG_DIR / "deaths_and_routes.png")
    if len(st):
        plot_survival(st, FIG_DIR / "survival.png")
    if len(ma):
        plot_deaths_by_age(ma, FIG_DIR / "deaths_by_age.png")
    if len(pt):
        plot_prevalence(pt, FIG_DIR / "prevalence_vs_phia.png")
    print(f"\nfigures -> {FIG_DIR}")


if __name__ == "__main__":
    main()
