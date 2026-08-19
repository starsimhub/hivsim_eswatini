"""Exp 015 — VMMC prevalence-target fix: broken (upstream) vs fixed (subclass).

Runs run_sims.make_sim unchanged except for the VMMC implementation, N seeds
each, and records per run:
  - circumcision coverage by 5-yr age bin over time (vs SHIMS3 2021 targets)
  - HIV prevalence 15-49 (overall, male, female) over time

Scalars + trajectories -> outputs/results.jsonl (one row per run, appended).
Figures -> figures/{vmmc_coverage_by_age,hiv_prevalence}.png.

Usage (from repo root):
    python experiments/015_vmmc_prevalence_target/run.py --n_seeds 10
    python experiments/015_vmmc_prevalence_target/run.py --plot_only
"""

import argparse
import json
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import starsim as ss
import stisim as sti

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from run_sims import make_sim

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))
from vmmc_prevalence_target import VMMCPrevalenceTarget

FIG = HERE / "figures"
OUT = HERE / "outputs"
JSONL = OUT / "results.jsonl"

# 11 age bins matching data/vmmc_coverage.csv
CSV_BINS = [(10, 15), (15, 20), (20, 25), (25, 30), (30, 35), (35, 40),
            (40, 45), (45, 50), (50, 55), (55, 60), (60, 65)]
# SHIMS3 2021 targets (p_vmmc) for the reference overlay
TARGET_2021 = {(10, 15): 0.85, (15, 20): 0.738, (20, 25): 0.604, (25, 30): 0.439,
               (30, 35): 0.345, (35, 40): 0.298, (40, 45): 0.325, (45, 50): 0.267,
               (50, 55): 0.252, (55, 60): 0.241, (60, 65): 0.166}

ARMS = {"broken_upstream": sti.VMMC, "fixed_prevalence_target": VMMCPrevalenceTarget}
COLORS = {"broken_upstream": "#d62728", "fixed_prevalence_target": "#1f77b4"}
SNAP_YEAR = 2021


def annual_last(years, vals):
    """Collapse a per-step series to one value per calendar year (last in year)."""
    yi = np.floor(np.asarray(years)).astype(int)
    uy = np.unique(yi)
    out = np.array([vals[np.where(yi == yr)[0][-1]] for yr in uy])
    return uy, out


class VMMCDiag(ss.Analyzer):
    """Circumcision prevalence by CSV age bin + HIV prevalence 15-49 by sex."""
    def __init__(self, **kw):
        super().__init__(**kw)
        self.years = []
        self.circ = {b: [] for b in CSV_BINS}
        self.prev_all, self.prev_m, self.prev_f = [], [], []

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result("_dummy", dtype=float, scale=False))

    def _find_vmmc(self):
        ivs = self.sim.interventions
        items = ivs.values() if hasattr(ivs, "values") else ivs
        for iv in items:
            if hasattr(iv, "circumcised"):
                return iv
        return None

    def step(self):
        sim = self.sim
        ppl = sim.people
        hiv = sim.diseases.hiv
        alive = ppl.alive
        self.years.append(float(sim.t.yearvec[sim.ti]))

        iv = self._find_vmmc()
        for lo, hi in CSV_BINS:
            males = alive & ppl.male & (ppl.age >= lo) & (ppl.age < hi)
            n = males.count()
            if n > 0 and iv is not None:
                self.circ[(lo, hi)].append(float((males & iv.circumcised).count() / n))
            else:
                self.circ[(lo, hi)].append(np.nan)

        for store, mask in [
            (self.prev_all, alive & (ppl.age >= 15) & (ppl.age < 50)),
            (self.prev_m,   alive & ppl.male & (ppl.age >= 15) & (ppl.age < 50)),
            (self.prev_f,   alive & ppl.female & (ppl.age >= 15) & (ppl.age < 50)),
        ]:
            store.append(float(np.mean(hiv.infected[mask])) if mask.count() > 0 else np.nan)


def run_one(arm, vmmc_class, seed):
    sim = make_sim(seed=seed, verbose=-1, vmmc_class=vmmc_class, analyzers=[VMMCDiag()])
    sim.run()
    # Retrieve the sim's own analyzer instance (starsim steps its copy, not the
    # object we passed in), matching how exp 013 reads analyzers back.
    anas = sim.analyzers.values() if hasattr(sim.analyzers, "values") else sim.analyzers
    diag = next(a for a in anas if isinstance(a, VMMCDiag))

    yrs, _ = annual_last(diag.years, diag.prev_all)
    circ_annual = {f"{lo}_{hi}": annual_last(diag.years, diag.circ[(lo, hi)])[1].tolist()
                   for lo, hi in CSV_BINS}
    return dict(
        arm=arm, seed=int(seed),
        years=[int(y) for y in yrs],
        prev_all=annual_last(diag.years, diag.prev_all)[1].tolist(),
        prev_m=annual_last(diag.years, diag.prev_m)[1].tolist(),
        prev_f=annual_last(diag.years, diag.prev_f)[1].tolist(),
        circ_by_bin=circ_annual,
    )


def run_all(n_seeds):
    FIG.mkdir(exist_ok=True)
    OUT.mkdir(exist_ok=True)
    (OUT / "version_stamp.json").write_text(json.dumps(
        {"starsim": ss.__version__, "stisim": sti.__version__, "n_seeds": n_seeds}, indent=2))
    if JSONL.exists():
        JSONL.unlink()
    for arm, cls in ARMS.items():
        for seed in range(1, n_seeds + 1):
            print(f"  {arm} seed {seed}/{n_seeds}...")
            res = run_one(arm, cls, seed)
            with JSONL.open("a") as f:
                f.write(json.dumps(res) + "\n")
    print("Runs complete ->", JSONL)


def load():
    return [json.loads(l) for l in JSONL.read_text().splitlines() if l.strip()]


def _mean_at(rows, key_fn, year):
    """Mean across rows of a per-year series, evaluated at `year`."""
    vals = []
    for r in rows:
        yrs = np.array(r["years"])
        i = int(np.argmin(np.abs(yrs - year)))
        v = key_fn(r)[i]
        if not np.isnan(v):
            vals.append(v)
    return np.mean(vals) if vals else np.nan


def plot():
    rows = load()

    # --- Figure 1: circumcision coverage by age at SNAP_YEAR vs SHIMS3 targets ---
    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(CSV_BINS))
    w = 0.38
    for i, arm in enumerate(ARMS):
        rs = [r for r in rows if r["arm"] == arm]
        means = [_mean_at(rs, lambda r, k=f"{lo}_{hi}": np.array(r["circ_by_bin"][k]), SNAP_YEAR) * 100
                 for lo, hi in CSV_BINS]
        ax.bar(x + (i - 0.5) * w, means, w, color=COLORS[arm], label=arm)
    ax.plot(x, [TARGET_2021[b] * 100 for b in CSV_BINS], "k*", markersize=13,
            label="SHIMS3 2021 target")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{lo}-{hi - 1}" for lo, hi in CSV_BINS], rotation=45)
    ax.set_ylabel("Circumcision coverage (%)")
    ax.set_xlabel("Age group")
    ax.set_title(f"Exp 015 — VMMC coverage by age ({SNAP_YEAR}): upstream vs prevalence-target fix")
    ax.legend()
    ax.set_ylim(0, 105)
    fig.tight_layout()
    fig.savefig(FIG / "vmmc_coverage_by_age.png", dpi=150, bbox_inches="tight")
    print("Saved", FIG / "vmmc_coverage_by_age.png")

    # --- Figure 2: HIV prevalence 15-49 (overall + male) over time ---
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    for ax, (field, title) in zip(axes, [("prev_all", "HIV prevalence 15-49 (all)"),
                                          ("prev_m", "HIV prevalence 15-49 (male)")]):
        for arm in ARMS:
            rs = [r for r in rows if r["arm"] == arm]
            yrs = np.array(rs[0]["years"])
            arr = np.array([r[field] for r in rs]) * 100
            med = np.nanmedian(arr, axis=0)
            lo = np.nanpercentile(arr, 5, axis=0)
            hi = np.nanpercentile(arr, 95, axis=0)
            ax.plot(yrs, med, color=COLORS[arm], lw=2, label=arm)
            ax.fill_between(yrs, lo, hi, color=COLORS[arm], alpha=0.15, linewidth=0)
        ax.set_title(title)
        ax.set_xlabel("Year")
        ax.set_ylabel("Prevalence (%)")
        ax.set_xlim(1985, 2031)
        ax.legend(fontsize=9)
    fig.suptitle("Exp 015 — HIV prevalence: broken (upstream) vs fixed (prevalence-target) VMMC",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "hiv_prevalence.png", dpi=150, bbox_inches="tight")
    print("Saved", FIG / "hiv_prevalence.png")

    # --- Summary table ---
    lines = ["arm,circ_15_49_2021_pct,prev_all_2021_pct,prev_m_2021_pct,prev_f_2021_pct"]
    print("\n  arm                         circ15-49%  prevAll%  prevM%  prevF%")
    for arm in ARMS:
        rs = [r for r in rows if r["arm"] == arm]
        # 15-49 circumcision = mean over the 7 bins spanning 15-49, pop-unweighted proxy
        c1549 = np.mean([_mean_at(rs, lambda r, k=f"{lo}_{hi}": np.array(r["circ_by_bin"][k]), 2021)
                         for lo, hi in CSV_BINS if lo >= 15 and hi <= 50]) * 100
        pa = _mean_at(rs, lambda r: np.array(r["prev_all"]), 2021) * 100
        pm = _mean_at(rs, lambda r: np.array(r["prev_m"]), 2021) * 100
        pf = _mean_at(rs, lambda r: np.array(r["prev_f"]), 2021) * 100
        print(f"  {arm:28s} {c1549:8.1f}  {pa:7.2f}  {pm:6.2f}  {pf:6.2f}")
        lines.append(f"{arm},{c1549:.2f},{pa:.2f},{pm:.2f},{pf:.2f}")
    (OUT / "summary_table.csv").write_text("\n".join(lines) + "\n")
    print("\nSaved", OUT / "summary_table.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_seeds", type=int, default=10)
    parser.add_argument("--plot_only", action="store_true")
    args = parser.parse_args()
    if not args.plot_only:
        run_all(args.n_seeds)
    plot()
