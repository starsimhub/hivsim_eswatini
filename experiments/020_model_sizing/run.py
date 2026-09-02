"""Exp 020 — Size the model: population, replicates, and the rare-event floor.

Two parts, runnable independently (`--part size|repl|all`).

  size   N in {5k, 10k, 20k, 50k} x {high, low} transmission. Expected infected
         agents per PHIA stratum (corrected sex mapping), between-seed CV, wall
         time, peak RSS.
  repl   CV at three parameter points -> replicates per parameter set.

Memory is the binding constraint, not cores: the laptop has 12 logical cores but
~6 GB free of 34, and peak RSS above N=10,000 was unmeasured before this. Part A
runs N groups sequentially with a worker cap that falls as N rises, and records
peak RSS per sim so the question is answered rather than guessed. Every cell
writes its own parquet before returning, so an OOM costs only the cells in
flight.

Outputs
  outputs/sims/{part}__{tag}__{seed}.parquet   per-run, written as each finishes
  outputs/size.parquet, repl.parquet           concatenated
  outputs/expected_counts.csv   infected agents per PHIA stratum, by N and pset
  outputs/floor.csv             per (N, pset): thinnest stratum, n below 5/10
  outputs/scaling.csv           wall time and peak RSS vs N
  outputs/replicates.csv        CV and recommended replicate count per point
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
from analyzers import PopByAgeSex  # noqa: E402

STOP = 2026

PARAM_SETS = {
    # continuity with 016-019, and the favourable case for the count floor
    "high_transmission": {"beta_m2f": 0.0139, "rel_init_prev": 0.49},
    # just above 018's establishment threshold -- lowest prevalence, so the
    # floor bites hardest. This is the case that decides N, and 018 never ran it.
    "low_transmission":  {"beta_m2f": 0.008,  "rel_init_prev": 0.2},
}
REPL_POINTS = {
    "plausible":        {"beta_m2f": 0.0139, "rel_init_prev": 0.49},
    "low_transmission": {"beta_m2f": 0.008,  "rel_init_prev": 0.2},
    "default":          {"beta_m2f": 0.01,   "rel_init_prev": 0.2},
}
N_SWEEP = [5_000, 10_000, 20_000, 50_000]
N_REPL = 10_000

# Worker cap by N. Memory-driven, not core-driven: ~6 GB free and peak RSS
# unmeasured above 10k, so the large-N arms run few-at-a-time deliberately.
# Revised after measuring: peak RSS is 470/632/871 MB at N=5k/10k/20k, so
# ~1.5-2 GB at 50k. Three concurrent 50k sims fit in the ~6 GB available.
WORKER_CAP = {5_000: 10, 10_000: 8, 20_000: 4, 50_000: 3}

# PHIA Gender: 0 = Male, 1 = Female. Per experiments/008 and confirmed against
# the data -- Gender=1 gives 0.101 at 15-19 in 2007, SDHS 2006-07's figure for
# women. 018/run.py:144 had this inverted, which is half of why its floor
# estimate was wrong; the other half was computing it at high transmission only.
PHIA_SEX = {0: "m", 1: "f"}

KEEP = ("timevec", "hiv.new_deaths", "hiv.prevalence_15_49", "popagesex.")


def _peak_rss_mb():
    """Peak RSS of this process in MB. Windows has no resource.getrusage."""
    try:
        import psutil
        return psutil.Process().memory_info().peak_wset / 1e6
    except Exception:
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e3
        except Exception:
            return np.nan


# --- Runs --------------------------------------------------------------------

def _run(part, tag, seed, hiv_pars, n_agents):
    path = SIM_DIR / f"{part}__{tag}__{seed:03d}.parquet"
    if path.exists():
        return

    t0 = time.perf_counter()
    sim = make_sim(seed=seed, stop=STOP, verbose=-1,
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
    out["peak_rss_mb"] = _peak_rss_mb()
    out["pop_scale"] = float(out["popagesex.n_alive_total"].iloc[0]) / n_agents
    out.to_parquet(path, index=False)


def work_size(n_seeds):
    """Grouped by N so the caller can run groups sequentially under the cap."""
    groups = []
    for n in N_SWEEP:
        groups.append((n, [("size", f"n{n}_{pset}", s, PARAM_SETS[pset], n)
                           for pset in PARAM_SETS for s in range(n_seeds)]))
    return groups


def work_repl(n_seeds):
    return [("repl", name, s, pars, N_REPL)
            for name, pars in REPL_POINTS.items() for s in range(n_seeds)]


def load(part):
    files = sorted(SIM_DIR.glob(f"{part}__*.parquet"))
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


# --- Analysis ----------------------------------------------------------------

def expected_counts(df):
    """Expected infected AGENTS per PHIA target stratum, per (N, pset)."""
    tg = pd.read_csv(REPO_ROOT / "calibration_data" / "prevalence_by_age_sex.csv")
    rows = []
    for tag, g_all in df.groupby("tag"):
        n_agents = int(g_all["n_agents"].iloc[0])
        for _, t in tg.iterrows():
            yr, lo = int(t.Year), int(t["start age"])
            sex = PHIA_SEX[int(t.Gender)]
            icol = f"popagesex.n_infected_{sex}_{lo}_{lo + 5}"
            acol = f"popagesex.n_alive_{sex}_{lo}_{lo + 5}"
            if icol not in df.columns:
                continue
            g = g_all[g_all.timevec == yr]
            if not len(g):
                continue
            rows.append(dict(
                tag=tag, n_agents=n_agents, year=yr, sex=sex, age_low=lo,
                expected_agents=(g[icol] / g["pop_scale"]).mean(),
                model_prev=(g[icol] / g[acol].replace(0, np.nan)).mean(),
                phia_prev=t.NationalPrevalence,
            ))
    return pd.DataFrame(rows)


def floor_table(ec):
    rows = []
    for tag, g in ec.groupby("tag"):
        thin = g.nsmallest(1, "expected_agents").iloc[0]
        rows.append(dict(
            tag=tag, n_agents=int(g.n_agents.iloc[0]),
            min_expected=g.expected_agents.min(),
            thinnest=f"{thin.year} {thin.sex.upper()} {thin.age_low}-{thin.age_low + 5}",
            n_below_5=int((g.expected_agents < 5).sum()),
            n_below_10=int((g.expected_agents < 10).sum()),
            n_strata=len(g),
        ))
    return pd.DataFrame(rows).sort_values(["n_agents", "tag"])


def cv_table(df, by="tag"):
    rows = []
    for tag, g in df.groupby(by):
        tm = g.groupby("seed")["hiv.prevalence_15_49"].mean()
        cv = tm.std(ddof=1) / tm.mean() if tm.mean() else np.nan
        est = g[g.timevec == 2005].groupby("seed")["hiv.prevalence_15_49"].mean()
        band = ("3-5" if cv < 0.05 else "10-20" if cv < 0.20
                else "50+ or increase N")
        rows.append(dict(tag=tag, n_agents=int(g.n_agents.iloc[0]),
                         n_seeds=len(tm), mean_prev=tm.mean(), cv=cv,
                         replicates=band,
                         established=int((est > 0.05).sum()),
                         runtime_s=g.groupby("seed")["runtime_s"].first().mean(),
                         peak_rss_mb=g.groupby("seed")["peak_rss_mb"].first().mean()))
    return pd.DataFrame(rows).sort_values(["n_agents", "tag"])


def plot_scaling(sz, ec, path):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
    for pset, g in sz.groupby(sz.tag.str.split("_", n=1).str[1]):
        g = g.sort_values("n_agents")
        axes[0].plot(g.n_agents, g.runtime_s, "o-", label=pset)
        axes[1].plot(g.n_agents, g.cv * 100, "o-", label=pset)
    lin = sz.sort_values("n_agents")
    ref = lin.iloc[0]
    axes[0].plot(lin.n_agents, ref.runtime_s * lin.n_agents / ref.n_agents,
                 "k--", alpha=0.5, label="linear")
    axes[0].set(title="Wall time per sim", xlabel="N agents", ylabel="seconds",
                xscale="log", yscale="log")
    axes[1].set(title="Between-seed CV of mean prevalence", xlabel="N agents",
                ylabel="CV (%)", xscale="log")
    axes[1].axhspan(0, 5, color="green", alpha=0.08)
    axes[1].axhspan(5, 20, color="orange", alpha=0.08)

    for tag, g in ec.groupby("tag"):
        n = int(g.n_agents.iloc[0])
        pset = tag.split("_", 1)[1]
        axes[2].scatter([n] * len(g), g.expected_agents, s=8, alpha=0.5,
                        label=pset if n == ec.n_agents.min() else None)
    axes[2].axhline(5, color="r", ls="--", label="floor = 5 agents")
    axes[2].axhline(10, color="orange", ls=":", label="floor = 10")
    axes[2].set(title="Expected infected agents per PHIA stratum",
                xlabel="N agents", ylabel="expected agents",
                xscale="log", yscale="log")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Exp 020 — model sizing", y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


# --- Main --------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--part", default="all", choices=["all", "size", "repl"])
    p.add_argument("--n_seeds", type=int, default=10)
    p.add_argument("--max_n", type=int, default=None,
                   help="skip N above this, e.g. 20000 to defer the 50k arm")
    p.add_argument("--plot_only", action="store_true")
    args = p.parse_args()
    parts = ["size", "repl"] if args.part == "all" else [args.part]

    if not args.plot_only and "size" in parts:
        for n, work in work_size(args.n_seeds):
            if args.max_n and n > args.max_n:
                print(f"[size] skipping N={n} (--max_n {args.max_n})")
                continue
            todo = [w for w in work
                    if not (SIM_DIR / f"{w[0]}__{w[1]}__{w[2]:03d}.parquet").exists()]
            cap = WORKER_CAP.get(n, 4)
            print(f"[size] N={n}: {len(work) - len(todo)}/{len(work)} on disk; "
                  f"running {len(todo)} at {cap} workers")
            if todo:
                t0 = sc.tic()
                sc.parallelize(_run, iterarg=todo, ncpus=cap)
                sc.toc(t0, label=f"N={n}")

    if not args.plot_only and "repl" in parts:
        work = work_repl(args.n_seeds)
        todo = [w for w in work
                if not (SIM_DIR / f"{w[0]}__{w[1]}__{w[2]:03d}.parquet").exists()]
        print(f"[repl] {len(work) - len(todo)}/{len(work)} on disk; "
              f"running {len(todo)} at 8 workers")
        if todo:
            t0 = sc.tic()
            sc.parallelize(_run, iterarg=todo, ncpus=8)
            sc.toc(t0, label="repl")

    if "size" in parts:
        sz_raw = load("size")
        if len(sz_raw):
            sz_raw.to_parquet(OUT_DIR / "size.parquet", index=False)
            ec = expected_counts(sz_raw)
            ec.to_csv(OUT_DIR / "expected_counts.csv", index=False)
            fl = floor_table(ec)
            fl.to_csv(OUT_DIR / "floor.csv", index=False)
            cv = cv_table(sz_raw)
            cv.to_csv(OUT_DIR / "scaling.csv", index=False)
            print("\n=== Rare-event floor: expected infected agents per PHIA "
                  "stratum ===")
            print(fl.to_string(index=False))
            print("\n=== Scaling: CV, wall time, peak memory ===")
            print(cv.round(3).to_string(index=False))
            plot_scaling(cv, ec, FIG_DIR / "sizing.png")

    if "repl" in parts:
        rp = load("repl")
        if len(rp):
            rp.to_parquet(OUT_DIR / "repl.parquet", index=False)
            rt = cv_table(rp)
            rt.to_csv(OUT_DIR / "replicates.csv", index=False)
            print("\n=== Replicate count by parameter point (N = "
                  f"{N_REPL:,}) ===")
            print(rt.round(3).to_string(index=False))

    print(f"\nfigures -> {FIG_DIR}")


if __name__ == "__main__":
    main()
