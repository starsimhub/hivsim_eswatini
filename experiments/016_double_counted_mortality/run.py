"""Exp 016 — A/B: all-cause vs HIV-deleted background mortality.

data/eswatini_deaths.csv feeds ss.Deaths, which kills agents regardless of HIV
status, and its adult rates carry a 5.5x hump peaking in 2005 — the AIDS
epidemic. The HIV module then kills agents again via p_hiv_death and ti_zero.
This runs the model with and without that hump, same seeds and parameters, and
measures the difference.

See README.md for the question and the construction's assumptions.

Outputs:
  outputs/data_hiv_deleted/        alternate datafolder (HIV-deleted mortality)
  outputs/mortality_construction.csv  per age/sex/year: rate deleted, AIDS share
  outputs/sims/{arm}_{pset}_{seed}.parquet  per-run results, as each finishes
  outputs/results.parquet          consolidated
  outputs/scorecard.csv            headline A/B numbers
  figures/mortality_ab_{pset}.png
  figures/aids_share_by_age.png
  figures/mortality_curves_by_age.png
  figures/implied_vs_unaids.png

Usage (from repo root):
    python experiments/016_double_counted_mortality/run.py
    python experiments/016_double_counted_mortality/run.py --n_seeds 2
    python experiments/016_double_counted_mortality/run.py --plot_only
"""

from __future__ import annotations

import os
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                  NUMEXPR_NUM_THREADS="1", MKL_NUM_THREADS="1")

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import yaml

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
OUT_DIR, FIG_DIR = EXP_DIR / "outputs", EXP_DIR / "figures"
SIM_DIR = OUT_DIR / "sims"
for d in (OUT_DIR, FIG_DIR, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(EXP_DIR))
from run_sims import make_sim  # noqa: E402
from mortality_construction import (  # noqa: E402
    build_hiv_deleted, deleted_fraction, make_datafolder,
)

CFG = yaml.safe_load((EXP_DIR / "config.yaml").read_text())
SIM_STOP = CFG["model"]["stop"]
DEATHS_FILE = "eswatini_deaths.csv"
ALT_DATAFOLDER = OUT_DIR / "data_hiv_deleted"
ARMS = ["all_cause", "hiv_deleted"]

# Two parameter sets, because the fix's net effect depends on epidemic size.
# The HIV-deleted arm removes AIDS mortality from the whole population and
# expects the HIV module to supply it back — which it can only do if the
# epidemic is roughly the right size. At defaults the model produces ~1.4%
# prevalence in 2005 against a real ~26%, so the fix would be judged against a
# model that cannot compensate. `high_transmission` takes the two parameters
# 014 found actually drive the epidemic, at its best-fitting draw.
PARAM_SETS = {
    "default": {},
    "high_transmission": {"beta_m2f": 0.0139, "rel_init_prev": 0.49},
}
# Metric 5 validates the mortality *data*, not the model, so it needs a
# population denominator that tracks reality. Use the arm/pset that does.
PSET_REF = "default"


class PopByAgeSex(ss.Analyzer):
    """Population and infection counts needed for the mortality attribution.

    Everything here could in principle be pieced together from existing results,
    but the pieces are scattered across bins that don't line up with the
    mortality data's 5-year age/sex structure. Recording them directly keeps the
    metric-5 denominator exact.
    """
    AGE_BINS = [(a, a + 5) for a in range(0, 100, 5)]

    def __init__(self, *args, name="popagesex", **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name

    def init_results(self):
        super().init_results()
        res = []
        for sex in ("f", "m"):
            for lo, hi in self.AGE_BINS:
                res.append(ss.Result(f"n_alive_{sex}_{lo}_{hi}", dtype=int, scale=True))
            res.append(ss.Result(f"prevalence_{sex}_15_49", dtype=float, scale=False))
        res.append(ss.Result("n_alive_total", dtype=int, scale=True))
        res.append(ss.Result("n_infected_total", dtype=int, scale=True))
        res.append(ss.Result("new_infections_total", dtype=int, scale=True))
        self.define_results(*res)

    def step(self):
        sim, ti = self.sim, self.ti
        ppl, hiv = sim.people, sim.diseases.hiv
        alive = ppl.alive

        for sex, sex_bool in (("f", ppl.female), ("m", ppl.male)):
            for lo, hi in self.AGE_BINS:
                in_bin = alive & sex_bool & (ppl.age >= lo) & (ppl.age < hi)
                self.results[f"n_alive_{sex}_{lo}_{hi}"][ti] = in_bin.count()
            adults = alive & sex_bool & (ppl.age >= 15) & (ppl.age < 50)
            if adults.count() > 0:
                self.results[f"prevalence_{sex}_15_49"][ti] = float(np.mean(hiv.infected[adults]))

        self.results["n_alive_total"][ti] = alive.count()
        self.results["n_infected_total"][ti] = (alive & hiv.infected).count()
        self.results["new_infections_total"][ti] = (alive & (hiv.ti_infected == ti)).count()


# --- Runs -------------------------------------------------------------------

KEEP = ("timevec", "hiv.new_deaths", "hiv.prevalence_15_49", "popagesex.")


def run_one(arm: str, pset: str, seed: int) -> pd.DataFrame:
    """One run; writes before returning so a spot reclamation is recoverable."""
    path = SIM_DIR / f"{arm}_{pset}_{seed:03d}.parquet"
    if path.exists():
        return pd.read_parquet(path)

    datafolder = str(ALT_DATAFOLDER) if arm == "hiv_deleted" else None
    sim = make_sim(seed=seed, stop=SIM_STOP, verbose=-1, datafolder=datafolder,
                   hiv_pars=dict(PARAM_SETS[pset]) or None,
                   analyzers=[PopByAgeSex()])
    sim.run()
    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns if any(c == k or c.startswith(k) for k in KEEP)]
    out = df[keep].copy()
    out["arm"], out["pset"], out["seed"] = arm, pset, seed
    out.to_parquet(path, index=False)
    return out


# --- Metrics ----------------------------------------------------------------

def unaids_deaths() -> pd.DataFrame:
    """UNAIDS Spectrum annual AIDS deaths — the target, from the correct file."""
    d = pd.read_csv(REPO_ROOT / "data" / "eswatini_hiv_calib.csv",
                    usecols=["time", "hiv.new_deaths"]).dropna()
    return d.rename(columns={"time": "year", "hiv.new_deaths": "unaids_deaths"})


def attribute_deaths(df: pd.DataFrame) -> pd.DataFrame:
    """Split deaths among PLHIV into HIV-module vs background-module.

    Uses the stock identity, since the HIV module is the only route out of the
    infected state apart from death:
        deaths_among_plhiv[t] = n_infected[t-1] + new_infections[t] - n_infected[t]
    then subtracts the HIV module's own reported deaths. Caveat for the SUMMARY:
    migration is on, so agents entering or leaving while infected also move this
    quantity and are not separated out here.
    """
    rows = []
    for (arm, pset, seed), g in df.groupby(["arm", "pset", "seed"]):
        g = g.sort_values("timevec")
        n_inf = g["popagesex.n_infected_total"].values
        new_inf = g["popagesex.new_infections_total"].values
        hiv_deaths = g["hiv.new_deaths"].values
        total = np.full(len(g), np.nan)
        total[1:] = n_inf[:-1] + new_inf[1:] - n_inf[1:]
        rows.append(pd.DataFrame({
            "arm": arm, "pset": pset, "seed": seed, "year": g["timevec"].values,
            "plhiv_deaths_total": total,
            "hiv_module_deaths": hiv_deaths,
            "background_deaths_plhiv": total - hiv_deaths,
        }))
    return pd.concat(rows, ignore_index=True)


def implied_aids_deaths(df: pd.DataFrame, construction: pd.DataFrame) -> pd.DataFrame:
    """Metric 5: deleted mortality rate x population, summed — vs UNAIDS.

    Validates the construction. If the rates we deleted imply roughly the AIDS
    deaths UNAIDS reports, the interpolation is picking up the right magnitude;
    substantially more means it is over-deleting.
    """
    base = df[(df.arm == "all_cause") & (df.pset == PSET_REF)]
    pop_cols = {c: c.replace("popagesex.n_alive_", "")
                for c in df.columns if c.startswith("popagesex.n_alive_")
                and c != "popagesex.n_alive_total"}
    pop = (base.groupby("timevec")[list(pop_cols)].mean()
           .rename(columns=pop_cols).stack().rename("pop").reset_index())
    pop = pop.rename(columns={"level_1": "stratum", "timevec": "year"})
    pop[["sex", "age_low", "age_high"]] = pop["stratum"].str.extract(r"([fm])_(\d+)_(\d+)")
    pop["Sex"] = pop["sex"].map({"f": "Female", "m": "Male"})
    pop["AgeStart"] = pop["age_low"].astype(int)

    # The mortality data is by SINGLE year of age; the population is in 5-year
    # bins. Merging them directly applies age 0's rate to the whole 0-4 bin,
    # which massively over-counts infants (age 0 mortality is ~5x ages 1-4) and
    # inflates the apparent child share of AIDS deaths. Average the single-year
    # rates within each bin first.
    c = construction[["Time", "Sex", "AgeStart", "deleted_rate"]].copy()
    c["AgeStart"] = (c["AgeStart"] // 5) * 5
    c = c.groupby(["Time", "Sex", "AgeStart"], as_index=False)["deleted_rate"].mean()

    # The mortality data is decadal, so the deleted rate has to be spread across
    # intervening years. Linear interpolation, matching how starsim interpolates
    # the rate data itself — a step/carry-forward would misstate the years
    # between anchors, which is most of the epidemic.
    years = np.array(sorted(pop.year.unique()), dtype=float)
    filled = []
    for (sex, age), g in c.groupby(["Sex", "AgeStart"]):
        g = g.sort_values("Time")
        filled.append(pd.DataFrame({
            "Sex": sex, "AgeStart": age, "year": years,
            "deleted_rate": np.interp(years, g["Time"].values, g["deleted_rate"].values),
        }))
    c = pd.concat(filled, ignore_index=True)

    pop["year"] = pop["year"].astype(float)
    m = pop.merge(c, on=["year", "Sex", "AgeStart"], how="left").fillna({"deleted_rate": 0.0})
    m["implied_aids_deaths"] = m["deleted_rate"] * m["pop"]
    return m.groupby("year")["implied_aids_deaths"].sum().reset_index()


# --- Figures ----------------------------------------------------------------

def plot_ab(df, deaths_attr, targets, path, pset):
    df = df[df.pset == pset]
    deaths_attr = deaths_attr[deaths_attr.pset == pset]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    colors = {"all_cause": "C3", "hiv_deleted": "C0"}
    labels = {"all_cause": "All-cause (as shipped)", "hiv_deleted": "HIV-deleted"}

    ax = axes[0, 0]
    for arm in ARMS:
        g = df[df.arm == arm].groupby("timevec")["hiv.prevalence_15_49"]
        ax.plot(g.mean().index, g.mean().values, color=colors[arm], lw=2, label=labels[arm])
        ax.fill_between(g.mean().index, g.min(), g.max(), color=colors[arm], alpha=0.15)
    ax.set_title("A. HIV prevalence 15–49", fontsize=10)
    ax.set_ylabel("Prevalence"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    ax = axes[0, 1]
    for arm in ARMS:
        g = df[df.arm == arm].groupby("timevec")["hiv.new_deaths"]
        ax.plot(g.mean().index, g.mean().values, color=colors[arm], lw=2, label=labels[arm])
    ax.plot(targets["year"], targets["unaids_deaths"], "k--", lw=1.5, label="UNAIDS target")
    ax.set_title("B. AIDS deaths reported by the HIV module", fontsize=10)
    ax.set_ylabel("Annual deaths"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    ax = axes[1, 0]
    for arm in ARMS:
        g = df[df.arm == arm].groupby("timevec")["popagesex.n_alive_total"]
        ax.plot(g.mean().index, g.mean().values, color=colors[arm], lw=2, label=labels[arm])
    ax.set_title("C. Total population (refutation check)", fontsize=10)
    ax.set_ylabel("Alive"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    ax = axes[1, 1]
    for arm in ARMS:
        g = deaths_attr[deaths_attr.arm == arm].groupby("year")
        ax.plot(g["background_deaths_plhiv"].mean().index,
                g["background_deaths_plhiv"].mean().values,
                color=colors[arm], lw=2, label=f"{labels[arm]} — background")
        ax.plot(g["hiv_module_deaths"].mean().index,
                g["hiv_module_deaths"].mean().values,
                color=colors[arm], lw=1.2, ls=":", label=f"{labels[arm]} — HIV module")
    ax.set_title("D. Deaths among PLHIV, by which module killed them", fontsize=10)
    ax.set_ylabel("Annual deaths"); ax.grid(alpha=0.3); ax.legend(fontsize=7)

    for ax in axes.flat:
        ax.set_xlabel("Year")
    fig.suptitle("Exp 016 — all-cause vs HIV-deleted background mortality  |  "
                 f"parameters: {pset}  ({df.seed.nunique()} seeds)", fontsize=12)
    fig.tight_layout(); fig.savefig(path, dpi=120, bbox_inches="tight"); plt.close(fig)


def plot_aids_share(construction, path):
    fig, ax = plt.subplots(figsize=(8, 5))
    for year, ls in [(1995, ":"), (2005, "-"), (2015, "--")]:
        for sex, col in [("Female", "C3"), ("Male", "C0")]:
            s = construction[(construction.Time == year) & (construction.Sex == sex)]
            s = s.sort_values("AgeStart")
            ax.plot(s.AgeStart, 100 * s.aids_share, ls=ls, color=col,
                    label=f"{sex} {year}", lw=1.6)
    ax.set_xlabel("Age"); ax.set_ylabel("Implied AIDS share of all-cause mortality (%)")
    ax.set_title("Implied AIDS share by age — shape is the validation\n"
                 "(peaks mid-30s, vanishes at 80+; nothing imposed this)", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=3)
    fig.tight_layout(); fig.savefig(path, dpi=120, bbox_inches="tight"); plt.close(fig)


def plot_mortality_curves(construction, path, xmax=2035):
    """Observed all-cause vs the non-AIDS counterfactual, per age and sex.

    The audit figure for the construction: the gap between the solid and dashed
    line in each panel is what gets attributed to AIDS. Reading it by age is how
    the misallocation shows up — the gap should be widest in the 25-45 panels and
    near-absent in the very young and very old, and it is not.
    """
    ages = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
    ncol, nrow = 6, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(19, 9), sharex=True)
    c = construction[construction.Time <= xmax]

    for ax, age in zip(axes.flat, ages):
        for sex, col in (("Female", "C3"), ("Male", "C0")):
            g = c[(c.AgeStart == age) & (c.Sex == sex)].sort_values("Time")
            if not len(g):
                continue
            ax.plot(g.Time, g.Value_all_cause, color=col, lw=1.8, marker="o", ms=3,
                    label=f"{sex} all-cause")
            ax.plot(g.Time, g.Value_non_aids, color=col, lw=1.4, ls="--", marker="s",
                    ms=2.5, alpha=0.85, label=f"{sex} non-AIDS")
            ax.fill_between(g.Time, g.Value_non_aids, g.Value_all_cause,
                            color=col, alpha=0.12)
        ax.set_yscale("log")
        ax.set_title(f"age {age}", fontsize=9)
        ax.grid(alpha=0.3, which="both")
        ax.tick_params(labelsize=7)

    for ax in axes.flat[len(ages):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Year", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("Mortality rate (log)", fontsize=8)
    axes.flat[0].legend(fontsize=6, loc="upper left")

    fig.suptitle("Exp 016 — observed all-cause mortality (solid) vs the non-AIDS "
                 "counterfactual (dashed). Shaded gap is what the construction "
                 "attributes to AIDS.", fontsize=12)
    fig.tight_layout(); fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)


def plot_implied(implied, targets, path):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(implied["year"], implied["implied_aids_deaths"], color="C0", lw=2,
            label="Implied by the deleted rates")
    ax.plot(targets["year"], targets["unaids_deaths"], "k--", lw=1.5, label="UNAIDS Spectrum")
    ax.set_xlabel("Year"); ax.set_ylabel("Annual AIDS deaths")
    ax.set_title("Metric 5 — does the deleted mortality imply the AIDS deaths\n"
                 "UNAIDS actually reports?", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=120, bbox_inches="tight"); plt.close(fig)


# --- Main -------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_seeds", type=int, default=CFG["model"]["n_seeds"])
    p.add_argument("--n_workers", type=int, default=None)
    p.add_argument("--plot_only", action="store_true")
    args = p.parse_args()

    print("Building HIV-deleted mortality...")
    all_cause = pd.read_csv(REPO_ROOT / "data" / DEATHS_FILE)
    hiv_deleted = build_hiv_deleted(all_cause)
    construction = deleted_fraction(all_cause, hiv_deleted)
    construction.to_csv(OUT_DIR / "mortality_construction.csv", index=False)
    make_datafolder(REPO_ROOT / "data", ALT_DATAFOLDER, hiv_deleted, DEATHS_FILE)
    n_changed = (construction.deleted_rate > 1e-12).sum()
    print(f"  {n_changed} of {len(construction)} rows lowered; "
          f"years {sorted(construction.loc[construction.deleted_rate > 1e-12, 'Time'].unique())}")

    if args.plot_only:
        df = pd.read_parquet(OUT_DIR / "results.parquet")
    else:
        work = [(arm, pset, seed) for arm in ARMS for pset in PARAM_SETS
                for seed in range(args.n_seeds)]
        done = len(list(SIM_DIR.glob("*.parquet")))
        if done:
            print(f"Resuming: {done} runs already on disk")
        print(f"Running {len(work)} sims (n_workers={args.n_workers or 'all'})...")
        t0 = sc.tic()
        dfs = sc.parallelize(run_one, iterarg=work, ncpus=args.n_workers)
        sc.toc(t0, label="A/B run")
        df = pd.concat(dfs, ignore_index=True)
        df.to_parquet(OUT_DIR / "results.parquet", index=False)

    targets = unaids_deaths()
    attr = attribute_deaths(df)
    attr.to_csv(OUT_DIR / "death_attribution.csv", index=False)
    implied = implied_aids_deaths(df, construction)
    implied.to_csv(OUT_DIR / "implied_aids_deaths.csv", index=False)

    # Scorecard at the epidemic peak and at a recent year
    rows = []
    for pset in PARAM_SETS:
      for year in (2005, 2021):
        r = {"pset": pset, "year": year}
        for arm in ARMS:
            g = df[(df.arm == arm) & (df.pset == pset) & (df.timevec == year)]
            r[f"prev_15_49_{arm}"] = g["hiv.prevalence_15_49"].mean()
            r[f"prev_f_{arm}"] = g["popagesex.prevalence_f_15_49"].mean()
            r[f"prev_m_{arm}"] = g["popagesex.prevalence_m_15_49"].mean()
            r[f"hiv_deaths_{arm}"] = g["hiv.new_deaths"].mean()
            r[f"pop_{arm}"] = g["popagesex.n_alive_total"].mean()
        t = targets[targets.year == year]
        r["unaids_deaths"] = float(t.unaids_deaths.iloc[0]) if len(t) else np.nan
        i = implied[implied.year == year]
        r["implied_aids_deaths"] = float(i.implied_aids_deaths.iloc[0]) if len(i) else np.nan
        rows.append(r)
    score = pd.DataFrame(rows)
    score.to_csv(OUT_DIR / "scorecard.csv", index=False)

    print("\n=== Scorecard ===")
    for _, r in score.iterrows():
        print(f"\n{int(r.year)}:")
        print(f"  prevalence 15-49 : {r.prev_15_49_all_cause:.4f} -> {r.prev_15_49_hiv_deleted:.4f}")
        print(f"    female         : {r.prev_f_all_cause:.4f} -> {r.prev_f_hiv_deleted:.4f}")
        print(f"    male           : {r.prev_m_all_cause:.4f} -> {r.prev_m_hiv_deleted:.4f}")
        print(f"  HIV-module deaths: {r.hiv_deaths_all_cause:,.0f} -> {r.hiv_deaths_hiv_deleted:,.0f}"
              f"   (UNAIDS {r.unaids_deaths:,.0f})")
        print(f"  population       : {r.pop_all_cause:,.0f} -> {r.pop_hiv_deleted:,.0f}")
        print(f"  implied AIDS deaths from deleted rates: {r.implied_aids_deaths:,.0f}"
              f"   (UNAIDS {r.unaids_deaths:,.0f})")

    for pset in PARAM_SETS:
        plot_ab(df, attr, targets, FIG_DIR / f"mortality_ab_{pset}.png", pset)
    plot_aids_share(construction, FIG_DIR / "aids_share_by_age.png")
    plot_mortality_curves(construction, FIG_DIR / "mortality_curves_by_age.png")
    plot_implied(implied, targets, FIG_DIR / "implied_vs_unaids.png")
    print(f"\nWrote outputs to {OUT_DIR} and figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
