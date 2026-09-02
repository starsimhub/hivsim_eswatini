"""Exp 022 — Pivoting survival gradient x suppression as a stock target.

2 x 2 factorial, fixed parameters and seeds:

  A_base    flat survival (upstream)  + VLS flow only   -> model-v1.2, control
  B_pivot   EMOD gradient             + VLS flow only    -> the survival pivot
  C_stock   flat survival             + VLS stock target -> the capability
  D_both    EMOD gradient             + VLS stock target -> interaction / v1.3

Outputs
  outputs/sims/{arm}__{seed}.parquet        per-run, written as each finishes
  outputs/results.parquet, deaths.parquet   concatenated
  outputs/fit_by_band.csv    prevalence bias by sex x age band, per arm
  outputs/scorecard.csv      MAE / bias / within-CI, both resolutions, per arm
  outputs/deaths.csv         peak AIDS deaths vs UNAIDS, per arm
  outputs/survival.csv       realized untreated survival by age at infection
  outputs/cascade.csv        population VLS and suppression-given-ART vs PHIA
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
OUT_DIR, FIG_DIR = EXP_DIR / "outputs", EXP_DIR / "figures"
SIM_DIR = OUT_DIR / "sims"
for d in (OUT_DIR, FIG_DIR, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
from run_sims import make_sim  # noqa: E402
from analyzers import PopByAgeSex, Cascade  # noqa: E402
from hiv_survival import AgeDependentSurvival, ARMS as SURV_ARMS  # noqa: E402
from vls_construction import build as vls_build, to_vls_coverage  # noqa: E402
from standard_figures import (load_targets, plot_prevalence_fit, scorecard,
                              fit_by_stratum, fit_15_49)  # noqa: E402

N_AGENTS = 10_000
STOP = 2026
N_SEEDS = 10
PARS = {"beta_m2f": 0.0139, "rel_init_prev": 0.49}
VLS_FILL_BACK_TO = 1985

AGE_BANDS = [(15, 25), (25, 35), (35, 45), (45, 65)]
SURV_CUTOFF = 2000   # infections up to here, so survival is not right-censored
UNAIDS_PEAK, UNAIDS_PEAK_YEAR = 11_000, 2004

# (survival arm in hiv_survival.ARMS, VLS applied as a stock target)
ARMS = {
    "A_base":  ("A_flat_13",   False),
    "B_pivot": ("E_grad_emod", False),
    "C_stock": ("A_flat_13",   True),
    "D_both":  ("E_grad_emod", True),
}

KEEP = ("timevec", "hiv.new_deaths", "hiv.prevalence_15_49", "hiv.n_infected",
        "hiv.new_infections", "n_alive", "popagesex.", "cascade.",
        "vls_stock_target.")


def vls_table():
    return to_vls_coverage(vls_build(), fill_back_to=VLS_FILL_BACK_TO)


def _run(arm, seed):
    path = SIM_DIR / f"{arm}__{seed:03d}.parquet"
    dpath = SIM_DIR / f"{arm}__{seed:03d}__deaths.parquet"
    if path.exists() and dpath.exists():
        return

    surv_arm, use_stock = ARMS[arm]
    tbl = vls_table()
    # Explicit since exp 022 closed: VLSStockTarget became the make_sim default
    # (model-v1.3), so leaving this implicit would give arms A and B a stock
    # target they were defined without, and arms C and D two of them.

    t0 = time.perf_counter()
    sim = make_sim(seed=seed, stop=STOP, verbose=-1, hiv_pars=dict(PARS),
                   hiv_class=partial(AgeDependentSurvival,
                                     latent_mult=SURV_ARMS[surv_arm]),
                   art_vls_coverage=tbl, vls_stock_target=use_stock,
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

    deaths = sim.diseases.hiv.death_records()
    deaths["arm"], deaths["seed"] = arm, seed
    deaths["year"] = (sim.t.yearvec[deaths["ti"].astype(int)]
                      if len(deaths) else [])
    deaths["infection_year"] = deaths["year"] - deaths["survival_years"]
    deaths.to_parquet(dpath, index=False)


def load(deaths=False):
    pat = "*__deaths.parquet" if deaths else "*.parquet"
    files = sorted(SIM_DIR.glob(pat))
    if not deaths:
        files = [f for f in files if not f.name.endswith("__deaths.parquet")]
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


# --- Metric 1: the bias decomposition that MAE hides -------------------------

def fit_by_band(df, tg):
    rows = []
    for arm, g in df.groupby("arm"):
        st = fit_by_stratum(g, tg=tg)
        if not len(st):
            continue
        st["err"] = st.model - st.phia
        st["band"] = pd.cut(st.age_low, [14, 24, 34, 44, 64],
                            labels=[f"{lo}-{hi - 1}" for lo, hi in AGE_BANDS])
        piv = st.pivot_table(index="band", columns="sex", values="err",
                             observed=True)
        for band in piv.index:
            rows.append(dict(arm=arm, band=str(band),
                             bias_f=piv.loc[band].get("f", np.nan),
                             bias_m=piv.loc[band].get("m", np.nan)))
    return pd.DataFrame(rows)


# --- Metric 2: deaths --------------------------------------------------------

def deaths_table(df):
    rows = []
    for arm, g in df.groupby("arm"):
        m = g.groupby("timevec")["hiv.new_deaths"].mean()
        rows.append(dict(arm=arm, peak=m.max(), peak_year=int(m.idxmax()),
                         pct_of_unaids=100 * m.max() / UNAIDS_PEAK,
                         deficit=UNAIDS_PEAK - m.max()))
    return pd.DataFrame(rows).sort_values("arm")


# --- Metric 3: realized survival ---------------------------------------------

def survival_table(dd):
    u = dd[(~dd.on_art) & (dd.age_at_infection >= 15)
           & (dd.infection_year <= SURV_CUTOFF)]
    rows = []
    for arm, g in u.groupby("arm"):
        for lo, hi in AGE_BANDS + [(15, 200)]:
            b = g[(g.age_at_infection >= lo) & (g.age_at_infection < hi)]
            rows.append(dict(arm=arm,
                             band="all" if lo == 15 and hi == 200 else f"{lo}-{hi - 1}",
                             n=len(b),
                             mean_years=b.survival_years.mean() if len(b) else np.nan))
    return pd.DataFrame(rows)


# --- Metric 4: cascade -------------------------------------------------------

PHIA_POP_VLS = {(2016, "m"): 0.623, (2016, "f"): 0.748,
                (2021, "m"): 0.824, (2021, "f"): 0.886}


def cascade_table(df):
    tbl = vls_table()
    rows = []
    for (arm, yr), g in df[df.timevec.isin([2016, 2021])].groupby(["arm", "timevec"]):
        for sex in ("m", "f"):
            t = tbl[(tbl.Gender == sex) & (tbl.Year <= yr)].sort_values("Year")
            rows.append(dict(
                arm=arm, year=int(yr), sex=sex,
                vls_plhiv=g[f"cascade.p_vls_{sex}_15_49"].mean(),
                phia_vls_plhiv=PHIA_POP_VLS[(int(yr), sex)],
                vls_given_art=g[f"cascade.p_vls_given_art_{sex}_15plus"].mean(),
                input_vls_given_art=float(t.p_vls.iloc[-1]) if len(t) else np.nan,
            ))
    out = pd.DataFrame(rows)
    out["plhiv_err"] = out.vls_plhiv - out.phia_vls_plhiv
    out["given_art_err"] = out.vls_given_art - out.input_vls_given_art
    return out


# --- Figures -----------------------------------------------------------------

COLORS = {"A_base": "#4c78a8", "B_pivot": "#e45756",
          "C_stock": "#54a24b", "D_both": "#b279a2"}


def plot_overview(df, dt, path):
    fig, axes = plt.subplots(1, 4, figsize=(21, 4.5))
    for col, title, ax in (("hiv.prevalence_15_49", "HIV prevalence 15-49", axes[0]),
                           ("hiv.new_deaths", "AIDS deaths", axes[1]),
                           ("hiv.new_infections", "New infections", axes[2])):
        for arm, g in df.groupby("arm"):
            m = g.groupby("timevec")[col].agg(["mean", "std"])
            ax.plot(m.index, m["mean"], label=arm, color=COLORS.get(arm))
            ax.fill_between(m.index, m["mean"] - m["std"], m["mean"] + m["std"],
                            alpha=0.12, color=COLORS.get(arm))
        if col == "hiv.new_deaths":
            ax.plot([UNAIDS_PEAK_YEAR], [UNAIDS_PEAK], "k*", ms=14,
                    label="UNAIDS peak")
        ax.set(title=title, xlabel="year")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)

    ax = axes[3]
    for arm, g in df.groupby("arm"):
        for sex, ls in (("f", "-"), ("m", "--")):
            col = f"cascade.p_vls_given_art_{sex}_15plus"
            if col not in df.columns:
                continue
            m = g.groupby("timevec")[col].mean()
            ax.plot(m.index, m, ls, color=COLORS.get(arm), alpha=0.85,
                    label=f"{arm} {sex}")
    tbl = vls_table()
    for sex, mk in (("f", "ks"), ("m", "k^")):
        s = tbl[tbl.Gender == sex].sort_values("Year")
        ax.plot(s.Year, s.p_vls, mk, ms=8, mfc="none", label=f"input {sex}")
    ax.set(title="Suppression among the treated\n(black = input series)",
           xlabel="year", xlim=(2004, 2026), ylim=(0.85, 1.005))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, ncol=2)

    fig.suptitle("Exp 022 — pivoting survival x suppression stock target", y=1.03)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_bias(fb, path):
    bands = [f"{lo}-{hi - 1}" for lo, hi in AGE_BANDS]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    width = 0.2
    x = np.arange(len(bands))
    for ax, sex, lbl in ((axes[0], "bias_f", "Women"), (axes[1], "bias_m", "Men")):
        for i, (arm, g) in enumerate(fb.groupby("arm")):
            vals = [g[g.band == b][sex].iloc[0] if len(g[g.band == b]) else np.nan
                    for b in bands]
            ax.bar(x + i * width, vals, width, label=arm, color=COLORS.get(arm))
        ax.axhline(0, color="k", lw=0.8)
        ax.set(title=f"{lbl}: model - PHIA by age band", xticks=x + 1.5 * width,
               xlabel="age band")
        ax.set_xticklabels(bands)
        ax.grid(alpha=0.3, axis="y")
    axes[0].set_ylabel("bias in prevalence")
    axes[0].legend(fontsize=8)
    fig.suptitle("Exp 022 — the defect MAE hides: bias by sex and age band", y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_seeds", type=int, default=N_SEEDS)
    p.add_argument("--n_workers", type=int, default=None)
    p.add_argument("--arms", nargs="*", default=list(ARMS))
    p.add_argument("--plot_only", action="store_true")
    args = p.parse_args()

    if not args.plot_only:
        work = [(a, s) for a in args.arms for s in range(args.n_seeds)]
        todo = [w for w in work
                if not (SIM_DIR / f"{w[0]}__{w[1]:03d}__deaths.parquet").exists()]
        print(f"{len(work) - len(todo)}/{len(work)} on disk; running {len(todo)}")
        if todo:
            t0 = sc.tic()
            sc.parallelize(_run, iterarg=todo, ncpus=args.n_workers)
            sc.toc(t0, label="022")

    df, dd = load(), load(deaths=True)
    if not len(df):
        print("no results")
        return
    df.to_parquet(OUT_DIR / "results.parquet", index=False)
    if len(dd):
        dd.to_parquet(OUT_DIR / "deaths.parquet", index=False)
    tg = load_targets()

    print("\n=== Metric 1: prevalence bias by sex and age band (model - PHIA) ===")
    fb = fit_by_band(df, tg)
    fb.to_csv(OUT_DIR / "fit_by_band.csv", index=False)
    print(fb.pivot_table(index="band", columns="arm",
                         values=["bias_f", "bias_m"]).round(3).to_string())

    print("\n=== Scorecard, both resolutions ===")
    cards = []
    for arm, g in df.groupby("arm"):
        cards.append(dict(arm=arm, **scorecard(g, tg=tg)))
    sc_df = pd.DataFrame(cards)
    sc_df.to_csv(OUT_DIR / "scorecard.csv", index=False)
    print(sc_df.round(4).to_string(index=False))

    print("\n=== Metric 2: AIDS deaths vs UNAIDS ===")
    dt = deaths_table(df)
    dt.to_csv(OUT_DIR / "deaths.csv", index=False)
    print(dt.round(1).to_string(index=False))

    if len(dd):
        print(f"\n=== Metric 3: realized untreated survival "
              f"(infections <= {SURV_CUTOFF}) ===")
        st = survival_table(dd)
        st.to_csv(OUT_DIR / "survival.csv", index=False)
        print(st.pivot_table(index="band", columns="arm", values="mean_years",
                             observed=True).round(2).to_string())

    print("\n=== Metric 4: cascade ===")
    ct = cascade_table(df)
    ct.to_csv(OUT_DIR / "cascade.csv", index=False)
    print(ct.round(4).to_string(index=False))

    for arm, g in df.groupby("arm"):
        plot_prevalence_fit(g, f"022 {arm}", FIG_DIR / f"prevalence_fit_{arm}.png",
                            tg=tg)
    plot_overview(df, dt, FIG_DIR / "overview.png")
    plot_bias(fb, FIG_DIR / "bias_by_band.png")
    print(f"\nfigures -> {FIG_DIR}")


if __name__ == "__main__":
    main()
