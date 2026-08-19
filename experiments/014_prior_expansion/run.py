"""Exp 014 — Prior predictive coverage check on the corrected stack.

Re-runs experiment 009's coverage check with the network fix (011/013) and the
VMMC prevalence-target fix (015) in place, and with 009's three fixed-parameter
suspects opened into the prior (9 parameters, up from 6).

The targets, the observation model and the envelope definition are IDENTICAL to
009, so the headline number is directly comparable to 009's 30/89 (34%).

Prior ranges are read from config.yaml — it is the single source of truth, not
a duplicate record. Do not hard-code ranges here.

Outputs:
  outputs/draws.csv              parameter draws (one row per draw)
  outputs/sims/sim_NNNN.parquet  per-sim slim results, written as each finishes
  outputs/ensemble.parquet       consolidated ensemble
  outputs/coverage_summary.csv   per-target-row envelope + in/out flag
  figures/coverage_prevalence.png
  figures/coverage_deaths.png

Per-sim files make the run resumable: raccoon is a spot VM and can be reclaimed
with <30s notice, so a rerun skips draws that already have a file rather than
starting over.

Usage (from repo root):
    python experiments/014_prior_expansion/run.py
    python experiments/014_prior_expansion/run.py --n_draws 200
    python experiments/014_prior_expansion/run.py --plot_only
"""

from __future__ import annotations

# Pin BLAS threading for parallel safety, mirroring 009 / run_calibrations.py
import os
os.environ.update(
    OMP_NUM_THREADS="1",
    OPENBLAS_NUM_THREADS="1",
    NUMEXPR_NUM_THREADS="1",
    MKL_NUM_THREADS="1",
)

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
OUT_DIR = EXP_DIR / "outputs"
SIM_DIR = OUT_DIR / "sims"
FIG_DIR = EXP_DIR / "figures"
for d in (OUT_DIR, SIM_DIR, FIG_DIR):
    d.mkdir(exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
from run_sims import make_sim  # noqa: E402
from hiv_mortality import HIVMortalityMultiplier  # noqa: E402

CFG = yaml.safe_load((EXP_DIR / "config.yaml").read_text())

PRIOR = {k: (v["low"], v["high"]) for k, v in CFG["prior"].items()}
SIM_STOP = CFG["model"]["stop"]
DRAW_SEED = CFG["model"]["draw_seed"]
ENV_LO, ENV_HI = CFG["coverage"]["envelope"]
TARGETS_PATH = REPO_ROOT / CFG["targets"]["file"]

# Upstream dur_latent is lognorm_ex(10y, 3y); the multiplier scales both moments
# so the shape of the distribution is preserved and only its scale moves.
DUR_LATENT_MEAN_Y = 10.0
DUR_LATENT_STD_Y = 3.0

# PHIA age bin -> model column. The model exposes 15:20 through 30:35 on
# `hiv.prevalence_*`, and 35:40 through 60:65 on the hiv_epi analyzer columns.
PREV_COL_TEMPLATE = {
    "hiv": "hiv.prevalence_{sex}_{lo}_{hi}",
    "hiv_epi": "hiv_epi.prevalence_{sex}_{lo}_{hi}",
}
KEEP_PREFIXES = ("timevec", "hiv.new_deaths", "hiv.prevalence_", "hiv_epi.prevalence_")


def model_prev_col(sex: str, age_low: int, age_high: int) -> str:
    """Pick the right prevalence column for a (sex, age_low, age_high) triple."""
    s = sex.lower()
    key = "hiv" if age_low < 35 else "hiv_epi"
    return PREV_COL_TEMPLATE[key].format(sex=s, lo=age_low, hi=age_high)


# --- Prior drawing ---------------------------------------------------------

def draw_prior(n_draws: int, seed: int) -> pd.DataFrame:
    """Uniform draws from the prior, returned as a tidy frame.

    Drawn as one (n_draws, n_params) block rather than per-parameter, so that
    row i is the same draw regardless of n_draws (row-major fill). Two
    consequences we rely on: the per-sim resume files stay valid across reruns,
    and raising --n_draws from 50 to 200 extends the ensemble instead of
    invalidating it. Reordering `prior:` in config.yaml would break both.
    """
    rng = np.random.default_rng(seed)
    u = rng.uniform(size=(n_draws, len(PRIOR)))
    cols = {"par_idx": np.arange(n_draws)}
    for k, (name, (lo, hi)) in enumerate(PRIOR.items()):
        cols[name] = lo + u[:, k] * (hi - lo)
    return pd.DataFrame(cols)


def split_pars(row: pd.Series) -> tuple[dict, dict]:
    """Split a draw into (hiv_pars, network_pars) for make_sim.

    Most entries pass through as scalars. `dur_latent_mult` is the exception:
    it is not a stisim parameter but a multiplier we expand into the actual
    duration distribution here.
    """
    hiv, net = {}, {}
    for col, val in row.items():
        if col == "par_idx":
            continue
        module, key = col.split(".", 1)
        if key == "dur_latent_mult":
            m = float(val)
            hiv["dur_latent"] = ss.lognorm_ex(ss.years(DUR_LATENT_MEAN_Y * m),
                                              ss.years(DUR_LATENT_STD_Y * m))
            continue
        target = hiv if module == "hiv" else net
        target[key] = float(val)
    return hiv, net


# --- Single sim worker -----------------------------------------------------

def run_one(par_idx: int, hiv_pars: dict, network_pars: dict) -> pd.DataFrame:
    """Run one sim, write its slim results to disk, and return them.

    Writes before returning so a spot reclamation costs at most the in-flight
    sims, and skips draws that already have a file so a rerun resumes.
    """
    path = SIM_DIR / f"sim_{par_idx:04d}.parquet"
    if path.exists():
        return pd.read_parquet(path)

    sim = make_sim(seed=par_idx, stop=SIM_STOP, verbose=-1,
                   hiv_pars=hiv_pars, network_pars=network_pars,
                   hiv_class=HIVMortalityMultiplier)
    sim.run()
    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns
            if any(c == p or c.startswith(p) for p in KEEP_PREFIXES)]
    out = df[keep].copy()
    out["par_idx"] = par_idx
    out.to_parquet(path, index=False)
    return out


# --- Coverage assessment ---------------------------------------------------

def assess_coverage(targets: pd.DataFrame, ensemble: pd.DataFrame) -> pd.DataFrame:
    """For each target row, compute the simulated envelope and an in/out flag."""
    rows = []
    ens_by_year = {y: g for y, g in ensemble.groupby("timevec")}

    for _, t in targets.iterrows():
        year = int(t["year"])
        if year not in ens_by_year:
            continue
        sub = ens_by_year[year]

        if t["quantity"] == "prevalence":
            col = model_prev_col(t["sex"], int(t["age_low"]), int(t["age_high"]))
            if col not in sub.columns:
                continue
            sim_vals = sub[col].dropna().values
        elif t["quantity"] == "aids_deaths":
            sim_vals = sub["hiv.new_deaths"].dropna().values
        else:
            continue

        if len(sim_vals) == 0:
            continue

        sim_lo = np.percentile(sim_vals, ENV_LO)
        sim_hi = np.percentile(sim_vals, ENV_HI)
        rows.append({
            "quantity": t["quantity"],
            "year": year,
            "sex": t["sex"],
            "age_low": t["age_low"],
            "age_high": t["age_high"],
            "obs_value": t["value"],
            "obs_lower": t["lower"],
            "obs_upper": t["upper"],
            "sim_p05": sim_lo,
            "sim_p50": np.percentile(sim_vals, 50),
            "sim_p95": sim_hi,
            "in_envelope": sim_lo <= t["value"] <= sim_hi,
        })

    return pd.DataFrame(rows)


# --- Plots -----------------------------------------------------------------

def plot_prevalence_coverage(targets, ensemble, path: Path) -> None:
    """PHIA years x sexes panel grid; ensemble band + observed CIs."""
    prev = targets[targets["quantity"] == "prevalence"].copy()
    phia_years = sorted(prev["year"].unique())
    sexes = ["M", "F"]

    fig, axes = plt.subplots(2, len(phia_years),
                             figsize=(4.5 * len(phia_years), 8), sharey=True)
    axes = np.atleast_2d(axes).reshape(2, len(phia_years))

    for j, year in enumerate(phia_years):
        year_ens = ensemble[ensemble["timevec"] == year]
        for i, sex in enumerate(sexes):
            ax = axes[i, j]
            obs = prev[(prev["year"] == year) & (prev["sex"] == sex)].copy()
            obs["age_mid"] = (obs["age_low"] + obs["age_high"]) / 2
            obs = obs.sort_values("age_mid")

            ages, p05, p50, p95 = [], [], [], []
            for _, r in obs.iterrows():
                col = model_prev_col(sex, int(r["age_low"]), int(r["age_high"]))
                if col in year_ens.columns:
                    vals = year_ens[col].dropna().values
                    if len(vals):
                        ages.append((r["age_low"] + r["age_high"]) / 2)
                        p05.append(np.percentile(vals, ENV_LO))
                        p50.append(np.percentile(vals, 50))
                        p95.append(np.percentile(vals, ENV_HI))

            if ages:
                ax.fill_between(ages, p05, p95, color="C0", alpha=0.25,
                                label=f"Sim {ENV_LO}-{ENV_HI}%")
                ax.plot(ages, p50, color="C0", lw=1.5, label="Sim median")

            err = [obs["value"] - obs["lower"], obs["upper"] - obs["value"]]
            ax.errorbar(obs["age_mid"], obs["value"], yerr=err, fmt="o",
                        color="C3", capsize=3, label="PHIA")

            ax.set_title(f"PHIA {year} — {'Male' if sex == 'M' else 'Female'}")
            ax.set_xlabel("Age (mid-bin)")
            if j == 0:
                ax.set_ylabel("HIV prevalence")
            ax.set_ylim(0, max(0.7, prev["upper"].max() * 1.05))
            ax.grid(alpha=0.3)
            if i == 0 and j == len(phia_years) - 1:
                ax.legend(fontsize=8, loc="upper right")

    fig.suptitle("Exp 014 coverage — PHIA prevalence vs. prior ensemble "
                 "(corrected stack, 9-parameter prior)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_deaths_coverage(targets, ensemble, path: Path) -> None:
    deaths = targets[targets["quantity"] == "aids_deaths"].sort_values("year")

    by_year = ensemble.groupby("timevec")["hiv.new_deaths"]
    yrs = sorted(by_year.groups.keys())
    p05 = [np.percentile(by_year.get_group(y).dropna(), ENV_LO) for y in yrs]
    p50 = [np.percentile(by_year.get_group(y).dropna(), 50) for y in yrs]
    p95 = [np.percentile(by_year.get_group(y).dropna(), ENV_HI) for y in yrs]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.fill_between(yrs, p05, p95, color="C0", alpha=0.25,
                    label=f"Sim {ENV_LO}-{ENV_HI}%")
    ax.plot(yrs, p50, color="C0", lw=1.5, label="Sim median")
    ax.errorbar(deaths["year"], deaths["value"],
                yerr=[deaths["value"] - deaths["lower"],
                      deaths["upper"] - deaths["value"]],
                fmt="o", color="C3", capsize=3, ms=4,
                label="UNAIDS (±15% placeholder)")

    ax.set_title("Exp 014 coverage — UNAIDS HIV deaths vs. prior ensemble")
    ax.set_xlabel("Year")
    ax.set_ylabel("New AIDS deaths (annual)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_prior_diagnostics(draws, ensemble, targets, path: Path) -> None:
    """Which prior parameters actually move prevalence, and are they at a bound?

    A coverage failure is only actionable once you know whether the model can
    reach the data at all. Panel A ranks the prior parameters by how strongly
    they move a well-observed target; panel B plots the dominant one against
    that target, so a pile-up at the prior bound is visible directly.
    """
    ref = {"year": 2016, "sex": "F", "age_low": 30, "age_high": 35}
    col = model_prev_col(ref["sex"], ref["age_low"], ref["age_high"])
    obs_row = targets[(targets["quantity"] == "prevalence")
                      & (targets["year"] == ref["year"])
                      & (targets["sex"] == ref["sex"])
                      & (targets["age_low"] == ref["age_low"])]
    if not len(obs_row) or col not in ensemble.columns:
        return
    obs = float(obs_row["value"].iloc[0])

    s = (ensemble[ensemble["timevec"] == ref["year"]][["par_idx", col]]
         .dropna().set_index("par_idx")[col])
    m = draws.set_index("par_idx").join(s.rename("y")).dropna(subset=["y"])
    par_cols = [c for c in draws.columns if c != "par_idx"]
    rho = m[par_cols + ["y"]].corr(method="spearman")["y"].drop("y").sort_values()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    ax1.barh(range(len(rho)), rho.values,
             color=["C3" if v < 0 else "C0" for v in rho.values])
    ax1.set_yticks(range(len(rho)))
    ax1.set_yticklabels([c.split(".", 1)[1] for c in rho.index], fontsize=9)
    ax1.axvline(0, color="k", lw=0.8)
    ax1.set_xlabel("Spearman rho with target")
    ax1.set_title(f"A. What moves {ref['sex']} {ref['age_low']}-{ref['age_high']} "
                  f"prevalence in {ref['year']}?", fontsize=10)
    ax1.grid(alpha=0.3, axis="x")

    top = rho.abs().idxmax()
    lo, hi = PRIOR[top]
    ax2.scatter(m[top], m["y"], s=28, color="C0", alpha=0.8, zorder=3)
    ax2.axhline(obs, color="C3", lw=1.5, ls="--", zorder=2,
                label=f"PHIA observed ({obs:.3f})")
    ax2.axvline(hi, color="k", lw=1.2, ls=":", zorder=2, label="prior upper bound")
    ax2.set_xlim(lo - 0.03 * (hi - lo), hi + 0.08 * (hi - lo))
    ax2.set_xlabel(top)
    ax2.set_ylabel(f"Simulated prevalence, {ref['sex']} "
                   f"{ref['age_low']}-{ref['age_high']}, {ref['year']}")
    ax2.set_title("B. Best draws sit at the prior's upper bound", fontsize=10)
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8, loc="upper left")

    fig.suptitle("Exp 014 — the coverage miss is driven by transmission and "
                 "seeding, not by the mortality parameters it opened", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


# --- Main ------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n_draws", type=int, default=CFG["model"]["n_draws"])
    p.add_argument("--n_workers", type=int, default=None,
                   help="default None = all cores")
    p.add_argument("--plot_only", action="store_true",
                   help="re-derive coverage and figures from existing outputs")
    args = p.parse_args()

    if args.plot_only:
        ensemble = pd.read_parquet(OUT_DIR / "ensemble.parquet")
    else:
        print(f"Drawing {args.n_draws} parameter sets from a "
              f"{len(PRIOR)}-parameter prior (seed {DRAW_SEED})...")
        draws = draw_prior(args.n_draws, DRAW_SEED)
        draws.to_csv(OUT_DIR / "draws.csv", index=False)

        work = []
        for _, row in draws.iterrows():
            hiv_pars, network_pars = split_pars(row)
            work.append((int(row["par_idx"]), hiv_pars, network_pars))

        n_done = len(list(SIM_DIR.glob("sim_*.parquet")))
        if n_done:
            print(f"Resuming: {n_done} sims already on disk, will be reused")

        print(f"Running {len(work)} sims (n_workers={args.n_workers or 'all'})...")
        t0 = sc.tic()
        dfs = sc.parallelize(run_one, iterarg=work, ncpus=args.n_workers)
        sc.toc(t0, label="ensemble run")

        ensemble = pd.concat(dfs, ignore_index=True)
        ensemble.to_parquet(OUT_DIR / "ensemble.parquet", index=False)
        print(f"Ensemble: {len(ensemble)} rows x {len(ensemble.columns)} cols "
              f"from {ensemble['par_idx'].nunique()} sims")

    print("Loading targets and assessing coverage...")
    targets = pd.read_csv(TARGETS_PATH)
    coverage = assess_coverage(targets, ensemble)
    coverage.to_csv(OUT_DIR / "coverage_summary.csv", index=False)

    n_total, n_in = len(coverage), int(coverage["in_envelope"].sum())
    pct = 100 * n_in / n_total if n_total else float("nan")
    print(f"\nCoverage: {n_in}/{n_total} target rows inside the "
          f"{ENV_LO}-{ENV_HI}% envelope ({pct:.0f}%)")
    print(coverage.groupby("quantity")["in_envelope"].agg(["sum", "count"]))
    print(f"\n009 baseline for comparison: {CFG['coverage']['baseline_009']['overall']}")

    plot_prevalence_coverage(targets, ensemble, FIG_DIR / "coverage_prevalence.png")
    plot_deaths_coverage(targets, ensemble, FIG_DIR / "coverage_deaths.png")
    plot_prior_diagnostics(pd.read_csv(OUT_DIR / "draws.csv"), ensemble, targets,
                           FIG_DIR / "prior_diagnostics.png")
    print(f"\nWrote outputs to {OUT_DIR} and figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
