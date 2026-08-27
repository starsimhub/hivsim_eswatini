"""Exp 018 — Adopt the 1.5.11 stack, then size the model.

Four parts, runnable independently (`--part adopt|size|repl|establish|all`).
Every part resumes from per-sim parquet in outputs/sims/, so an interrupted run
costs nothing.

  adopt      PrEP removal vs 017 arm C — the one adoption change with an
             unknown effect. The control is read from 017's outputs, not re-run.
  size       N in {5k, 10k, 20k, 50k}: per-stratum expected counts, CV, runtime.
  repl       CV at three parameter points -> replicates per parameter set.
  establish  beta_m2f x rel_init_prev grid -> which parts of the prior produce
             an epidemic at all.

Outputs
  outputs/sims/{part}__{tag}__{seed}.parquet   per-run, written as each finishes
  outputs/adopt.parquet, size.parquet, repl.parquet, establish.parquet
  outputs/expected_counts.csv    infected agents per PHIA stratum, by N
  outputs/replicates.csv         CV and recommended replicate count
  outputs/establishment.csv      fraction of seeds establishing, per grid cell
  outputs/runtime_scaling.csv    seconds per sim vs N
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
import yaml

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
EXP_017 = REPO_ROOT / "experiments" / "017_version_bump"
OUT_DIR, FIG_DIR = EXP_DIR / "outputs", EXP_DIR / "figures"
SIM_DIR = OUT_DIR / "sims"
for d in (OUT_DIR, FIG_DIR, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
from run_sims import make_sim  # noqa: E402
from analyzers import PopByAgeSex  # noqa: E402

CFG = yaml.safe_load((EXP_DIR / "config.yaml").read_text())
N_AGENTS_DEFAULT = 10_000
STOP_FULL = 2026

PARAM_SETS = {
    "default": {},
    "high_transmission": {"beta_m2f": 0.0139, "rel_init_prev": 0.49},
}
REPL_POINTS = {
    "plausible":        {"beta_m2f": 0.0139, "rel_init_prev": 0.49},
    "low_transmission": {"beta_m2f": 0.008,  "rel_init_prev": 0.2},
    "default":          {"beta_m2f": 0.01,   "rel_init_prev": 0.2},
}
N_SWEEP = [5_000, 10_000, 20_000, 50_000]
BETA_GRID = [0.006, 0.008, 0.010, 0.012, 0.014]
RIP_GRID = [0.1, 0.2, 0.35, 0.5]
ESTABLISH_STOP = 2006
ESTABLISH_YEAR = 2005
ESTABLISH_THRESHOLD = 0.05

KEEP = ("timevec", "hiv.new_deaths", "hiv.prevalence_15_49", "popagesex.")


# --- Runs --------------------------------------------------------------------

def _run(part, tag, seed, hiv_pars, n_agents, stop):
    """One sim; writes before returning so an interrupted batch is recoverable."""
    path = SIM_DIR / f"{part}__{tag}__{seed:03d}.parquet"
    if path.exists():
        return

    t0 = time.perf_counter()
    sim = make_sim(seed=seed, stop=stop, verbose=-1,
                   hiv_pars=dict(hiv_pars) or None,
                   analyzers=[PopByAgeSex()])
    sim.pars.n_agents = n_agents
    sim.run()
    elapsed = time.perf_counter() - t0

    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns if any(c == k or c.startswith(k) for k in KEEP)]
    out = df[keep].copy()
    out["part"], out["tag"], out["seed"] = part, tag, seed
    out["n_agents"], out["runtime_s"] = n_agents, elapsed
    out.to_parquet(path, index=False)


def work_adopt(n_seeds):
    return [("adopt", pset, s, PARAM_SETS[pset], N_AGENTS_DEFAULT, STOP_FULL)
            for pset in PARAM_SETS for s in range(n_seeds)]


def work_size(n_seeds):
    return [("size", f"n{n}", s, PARAM_SETS["high_transmission"], n, STOP_FULL)
            for n in N_SWEEP for s in range(n_seeds)]


def work_repl(n_seeds):
    return [("repl", name, s, pars, N_AGENTS_DEFAULT, STOP_FULL)
            for name, pars in REPL_POINTS.items() for s in range(n_seeds)]


def work_establish(n_seeds):
    return [("establish", f"b{b:.3f}_r{r:.2f}", s,
             {"beta_m2f": b, "rel_init_prev": r}, N_AGENTS_DEFAULT, ESTABLISH_STOP)
            for b in BETA_GRID for r in RIP_GRID for s in range(n_seeds)]


def load(part):
    files = sorted(SIM_DIR.glob(f"{part}__*.parquet"))
    if not files:
        return pd.DataFrame()
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    df.to_parquet(OUT_DIR / f"{part}.parquet", index=False)
    return df


# --- Analysis: expected counts per PHIA stratum ------------------------------

def expected_counts(df, label_col="n_agents"):
    """Infected AGENTS per PHIA target stratum — the rare-event floor check.

    The model's results are population-scaled, so a stratum's *agent* count is
    the scaled count divided by (total population / n_agents). Agent counts are
    what drive stochastic noise; the scaled numbers would flatter the model.
    """
    tg = pd.read_csv(REPO_ROOT / "calibration_data" / "prevalence_by_age_sex.csv")
    rows = []
    for label, g_all in df.groupby(label_col):
        n_agents = int(g_all["n_agents"].iloc[0])
        for _, t in tg.iterrows():
            yr, lo = int(t.Year), int(t["start age"])
            sex = "f" if t.Gender == 0 else "m"
            col = f"popagesex.n_alive_{sex}_{lo}_{lo + 5}"
            icol = f"popagesex.n_infected_{sex}_{lo}_{lo + 5}"
            g = g_all[g_all.timevec == yr]
            if col not in df.columns or not len(g):
                continue
            scale = g["popagesex.n_alive_total"].mean() / n_agents
            rows.append(dict(
                label=label, n_agents=n_agents, year=yr, sex=sex,
                age=f"{lo}-{lo + 5}", target_prev=t.NationalPrevalence,
                agents_in_stratum=g[col].mean() / scale,
                model_cases=(g[icol].mean() / scale) if icol in df.columns else np.nan,
                expected_cases_at_target=(g[col].mean() / scale) * t.NationalPrevalence,
            ))
    return pd.DataFrame(rows)


def cv_table(df, by):
    """Between-seed CV of trajectory-mean prevalence 15-49, 1995-2021."""
    w = df[(df.timevec >= 1995) & (df.timevec <= 2021)]
    rows = []
    for key, g in w.groupby(by):
        per_seed = g.groupby("seed")["hiv.prevalence_15_49"].mean()
        m, s = per_seed.mean(), per_seed.std(ddof=1)
        cv = s / m if m > 0 else np.inf
        rec = "3-5" if cv < 0.05 else ("10-20" if cv < 0.20 else "50+ or raise N")
        rows.append(dict(group=key, mean_prev=m, sd=s, cv=cv, n_seeds=len(per_seed),
                         recommended_replicates=rec))
    return pd.DataFrame(rows)


# --- Figures -----------------------------------------------------------------

def plot_adopt(cur, new, path):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for pset, ls in (("high_transmission", "-"), ("default", "--")):
        for d, c, lab in ((cur, "C0", "PrEP default ON (017 arm C)"),
                          (new, "C1", "PrEP removed (adopted)")):
            g = d[d.tag == pset] if "tag" in d.columns else d[d.pset == pset]
            if not len(g):
                continue
            for ax, col in zip(axes, ["hiv.prevalence_15_49", "hiv.new_deaths",
                                      "popagesex.n_alive_total"]):
                s = g.groupby("timevec")[col]
                ax.plot(s.median().index, s.median().values, ls, color=c, lw=2,
                        label=f"{lab}, {pset}")
    for ax, t in zip(axes, ["HIV prevalence 15-49", "AIDS deaths",
                            "Total population"]):
        ax.set_title(t, fontsize=11); ax.set_xlabel("year")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.suptitle("Exp 018 part 1 — removing the inherited PrEP default "
                 "(nobody chose it: 80% of FSW by 2025, ramp starting 2004)",
                 fontsize=12)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def plot_sizing(ec, rt, path):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for n, g in ec.groupby("n_agents"):
        v = np.sort(g.expected_cases_at_target.values)
        axes[0].plot(v, np.arange(1, len(v) + 1) / len(v), lw=2, label=f"N={n:,}")
    for thr, c in ((5, "C3"), (10, "C1"), (20, "C2")):
        axes[0].axvline(thr, color=c, ls=":", lw=1)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("expected infected agents in stratum")
    axes[0].set_ylabel("cumulative share of 54 PHIA strata")
    axes[0].set_title("A. Rare-event floor by N\n(dotted: 5, 10, 20)", fontsize=10)
    axes[0].legend(fontsize=8)

    axes[1].plot(rt.n_agents, rt["mean"], "o-", lw=2)
    lin = rt["mean"].iloc[0] * rt.n_agents / rt.n_agents.iloc[0]
    axes[1].plot(rt.n_agents, lin, "k--", lw=1, label="linear in N")
    axes[1].set_xscale("log"); axes[1].set_yscale("log")
    axes[1].set_xlabel("N agents"); axes[1].set_ylabel("seconds per sim")
    axes[1].set_title("B. Runtime scaling\n(above the dashed line = worse than linear)",
                      fontsize=10)
    axes[1].legend(fontsize=8)

    axes[2].plot(rt.n_agents, rt.cv * 100, "o-", lw=2)
    axes[2].axhline(5, color="C2", ls=":", lw=1)
    axes[2].axhline(20, color="C3", ls=":", lw=1)
    axes[2].set_xscale("log"); axes[2].set_xlabel("N agents")
    axes[2].set_ylabel("between-seed CV (%)")
    axes[2].set_title("C. Variance vs N\n(1/sqrt(N) if noise is demographic)", fontsize=10)
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle("Exp 018 part 2 — model sizing", fontsize=12)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def plot_establishment(est, path):
    piv = est.pivot(index="rel_init_prev", columns="beta_m2f", values="frac_established")
    fig, ax = plt.subplots(figsize=(7.5, 5))
    im = ax.imshow(piv.values, origin="lower", cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(piv.columns)), [f"{c:.3f}" for c in piv.columns])
    ax.set_yticks(range(len(piv.index)), [f"{r:.2f}" for r in piv.index])
    ax.set_xlabel("beta_m2f"); ax.set_ylabel("rel_init_prev")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            ax.text(j, i, f"{v:.0%}", ha="center", va="center",
                    color="white" if v < 0.6 else "black", fontsize=10)
    fig.colorbar(im, ax=ax, label="fraction of seeds establishing")
    ax.set_title("Exp 018 part 2c — where the epidemic establishes\n"
                 f"(prevalence 15-49 > {ESTABLISH_THRESHOLD:.0%} at {ESTABLISH_YEAR})",
                 fontsize=11)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


# --- Main --------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--part", default="all",
                   choices=["all", "adopt", "size", "repl", "establish"])
    p.add_argument("--n_seeds", type=int, default=10)
    p.add_argument("--n_workers", type=int, default=None)
    p.add_argument("--plot_only", action="store_true")
    args = p.parse_args()

    parts = ["adopt", "size", "repl", "establish"] if args.part == "all" else [args.part]
    builders = dict(adopt=work_adopt, size=work_size, repl=work_repl,
                    establish=work_establish)

    if not args.plot_only:
        for part in parts:
            n = 5 if part == "establish" else args.n_seeds
            work = builders[part](n)
            todo = [w for w in work
                    if not (SIM_DIR / f"{w[0]}__{w[1]}__{w[2]:03d}.parquet").exists()]
            print(f"\n[{part}] {len(work) - len(todo)}/{len(work)} on disk; "
                  f"running {len(todo)} (n_workers={args.n_workers or 'all'})")
            if todo:
                t0 = sc.tic()
                sc.parallelize(_run, iterarg=todo, ncpus=args.n_workers)
                sc.toc(t0, label=part)

    # --- Part 1: what did removing PrEP do? ---
    if "adopt" in parts:
        new = load("adopt")
        src = EXP_017 / "outputs" / "results.parquet"
        if len(new) and src.exists():
            cur = pd.read_parquet(src)
            cur = cur[cur.arm == "C_1_5_11_hiv_deleted"].rename(columns={"pset": "tag"})
            print("\n=== Part 1: PrEP removal (adopted) vs 017 arm C (PrEP on) ===")
            rows = []
            for pset in PARAM_SETS:
                for yr in (2005, 2015, 2021):
                    for m in ("hiv.prevalence_15_49", "hiv.new_deaths"):
                        a = cur[(cur.tag == pset) & (cur.timevec == yr)][m]
                        b = new[(new.tag == pset) & (new.timevec == yr)][m]
                        if not len(a) or not len(b):
                            continue
                        se = np.sqrt(a.std(ddof=1)**2/len(a) + b.std(ddof=1)**2/len(b))
                        z = (b.mean() - a.mean()) / se if se else np.nan
                        rows.append(dict(pset=pset, year=yr, metric=m,
                                         prep_on=a.mean(), prep_off=b.mean(),
                                         rel=(b.mean()/a.mean() - 1) if a.mean() else np.nan,
                                         z=z))
                        print(f"  {pset:<18} {yr} {m:<24} "
                              f"{a.mean():>10,.4f} -> {b.mean():>10,.4f} "
                              f"({rows[-1]['rel']:+.1%}, z={z:+.1f})")
            pd.DataFrame(rows).to_csv(OUT_DIR / "prep_removal.csv", index=False)
            plot_adopt(cur, new, FIG_DIR / "prep_removal.png")

    # --- Part 2a: population size ---
    if "size" in parts:
        size = load("size")
        if len(size):
            ec = expected_counts(size)
            ec.to_csv(OUT_DIR / "expected_counts.csv", index=False)
            rt = (size.groupby(["n_agents", "seed"])["runtime_s"].first()
                      .groupby("n_agents").agg(["mean", "std"]).reset_index())
            cv = cv_table(size, "n_agents").rename(columns={"group": "n_agents"})
            rt = rt.merge(cv[["n_agents", "cv"]], on="n_agents")
            rt.to_csv(OUT_DIR / "runtime_scaling.csv", index=False)

            print("\n=== Part 2a: population size ===")
            print(f"  {'N':>8} {'runtime':>9} {'CV':>7} {'<5 cases':>9} {'<10':>5} {'<20':>5}")
            for n, g in ec.groupby("n_agents"):
                r = rt[rt.n_agents == n].iloc[0]
                v = g.expected_cases_at_target
                print(f"  {n:>8,} {r['mean']:>8.0f}s {r['cv']:>6.1%} "
                      f"{(v < 5).sum():>9} {(v < 10).sum():>5} {(v < 20).sum():>5}"
                      f"   (of {len(v)} strata)")
            plot_sizing(ec, rt, FIG_DIR / "sizing.png")

    # --- Part 2b: replicate count ---
    if "repl" in parts:
        rep = load("repl")
        if len(rep):
            cv = cv_table(rep, "tag")
            cv.to_csv(OUT_DIR / "replicates.csv", index=False)
            print("\n=== Part 2b: replicate count ===")
            for _, r in cv.iterrows():
                print(f"  {r.group:<18} mean prev {r.mean_prev:.4f}  "
                      f"CV {r.cv:>6.1%}  -> {r.recommended_replicates} replicates")

    # --- Part 2c: establishment ---
    if "establish" in parts:
        est_raw = load("establish")
        if len(est_raw):
            g = est_raw[est_raw.timevec == ESTABLISH_YEAR]
            rows = []
            for tag, gg in g.groupby("tag"):
                b, r = tag.replace("b", "").split("_r")
                per_seed = gg.groupby("seed")["hiv.prevalence_15_49"].mean()
                rows.append(dict(beta_m2f=float(b), rel_init_prev=float(r),
                                 frac_established=float((per_seed > ESTABLISH_THRESHOLD).mean()),
                                 median_prev=float(per_seed.median()),
                                 n_seeds=len(per_seed)))
            est = pd.DataFrame(rows).sort_values(["rel_init_prev", "beta_m2f"])
            est.to_csv(OUT_DIR / "establishment.csv", index=False)
            print("\n=== Part 2c: establishment (fraction of seeds with "
                  f"prevalence > {ESTABLISH_THRESHOLD:.0%} at {ESTABLISH_YEAR}) ===")
            print(est.pivot(index="rel_init_prev", columns="beta_m2f",
                            values="frac_established").to_string())
            plot_establishment(est, FIG_DIR / "establishment.png")

    print(f"\nWrote outputs to {OUT_DIR} and figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
