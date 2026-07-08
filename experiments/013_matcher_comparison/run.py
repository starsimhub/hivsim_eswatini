"""Exp 013 — matcher comparison: decompose the network prevalence shift.

Runs run_sims.make_sim unchanged except for the network `match_method`, across
three matchers, N seeds each, and records per run:
  - HIV prevalence 15-49 and new-infections 15-49 trajectories (hiv_epi)
  - partnership VOLUME: active MF edges over time, mean/total lifetime partners
  - realized male-female age gap by female age group (2020 snapshot)

Scalars + trajectories -> outputs/results.jsonl (one row per run, appended).
Comparison figure -> figures/matcher_comparison.png.

Usage (from repo root):
    python experiments/013_matcher_comparison/run.py --n_seeds 10
    python experiments/013_matcher_comparison/run.py --plot_only   # re-plot from jsonl
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
from plot_dashboard import PartnershipSnapshot, annual_mean, annual_last

HERE = pathlib.Path(__file__).parent
FIG = HERE / "figures"
OUT = HERE / "outputs"
JSONL = OUT / "results.jsonl"

# old default -> strict/untapered -> new default (see README for the contrast logic)
MATCHERS = ["sort_bisect", "kdtree_nn", "closest_age_tapered_seeking"]
COLORS = {"sort_bisect": "#7f7f7f", "kdtree_nn": "#1f77b4",
          "closest_age_tapered_seeking": "#d62728"}
AGE_GAP_BINS = [(15, 24), (25, 34), (35, 49)]
DHS_GAP = {(15, 24): 8.6, (25, 34): 7.4, (35, 49): 7.7}
SNAP_YEAR = 2020


class EdgeCount(ss.Analyzer):
    """Record the number of active edges in the MF/SW network each step."""
    def __init__(self, network="structuredsexual", **kw):
        super().__init__(**kw)
        self.network = network
        self.years = []
        self.n = []

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result("_dummy", dtype=float, scale=False))

    def step(self):
        net = self.sim.networks[self.network]
        self.years.append(float(self.sim.t.yearvec[self.sim.ti]))
        self.n.append(int((net.edges.dur > 1).sum()))


def annual_sum(yearvec, vals):
    yi = np.floor(yearvec).astype(int)
    uy = np.unique(yi)
    return uy, np.array([np.nansum(vals[yi == yr]) for yr in uy])


def run_one(matcher, seed):
    analyzers = [PartnershipSnapshot(year=SNAP_YEAR), EdgeCount()]
    sim = make_sim(seed=seed, verbose=-1,
                   network_pars={"match_method": matcher}, analyzers=analyzers)
    sim.run()

    yv = sim.t.yearvec
    epi = sim.analyzers["hiv_epi"]
    yrs, prev = annual_last(yv, np.array(epi.results["prevalence_15_49"]))
    _, newinf = annual_sum(yv, np.array(epi.results["new_infections_15_49"]))
    snap = sim.analyzers["partnershipsnapshot"]
    ec = sim.analyzers["edgecount"]
    _, nedges = annual_mean(np.array(ec.years), np.array(ec.n))

    nw = sim.networks.structuredsexual
    ppl = sim.people
    part = nw.participant & ppl.alive
    deb = part & (ppl.age >= nw.debut)
    lpF = np.array(nw.lifetime_partners[deb & ppl.female])
    lpM = np.array(nw.lifetime_partners[deb & ppl.male])

    gaps = {}
    n_active = 0
    if snap.f_ages is not None:
        n_active = int(len(snap.f_ages))
        gap = snap.m_ages - snap.f_ages
        for lo, hi in AGE_GAP_BINS:
            m = (snap.f_ages >= lo) & (snap.f_ages < hi + 1)
            gaps[f"{lo}_{hi}"] = float(np.nanmean(gap[m])) if m.any() else None

    return dict(
        matcher=matcher, seed=int(seed),
        years=[int(y) for y in yrs],
        prev_15_49=[float(x) for x in prev],
        new_inf_15_49=[float(x) for x in newinf],
        n_edges=[float(x) for x in nedges],
        n_edges_years=[int(y) for y in np.unique(np.floor(np.array(ec.years)).astype(int))],
        mean_lp_f=float(lpF.mean()), mean_lp_m=float(lpM.mean()),
        total_pairs_f=float(lpF.sum()),
        n_active_pairs_2020=n_active,
        age_gaps=gaps,
    )


def run_all(n_seeds):
    FIG.mkdir(exist_ok=True)
    OUT.mkdir(exist_ok=True)
    (OUT / "version_stamp.json").write_text(json.dumps(
        {"starsim": ss.__version__, "stisim": sti.__version__, "n_seeds": n_seeds}, indent=2))
    if JSONL.exists():
        JSONL.unlink()  # fresh full run
    for matcher in MATCHERS:
        for seed in range(1, n_seeds + 1):
            print(f"  {matcher} seed {seed}/{n_seeds}...")
            res = run_one(matcher, seed)
            with JSONL.open("a") as f:
                f.write(json.dumps(res) + "\n")
    print("Runs complete ->", JSONL)


def load():
    rows = [json.loads(l) for l in JSONL.read_text().splitlines() if l.strip()]
    return rows


def _band(ax, x, arr2d, color, label):
    med = np.median(arr2d, axis=0)
    lo = np.percentile(arr2d, 5, axis=0)
    hi = np.percentile(arr2d, 95, axis=0)
    ax.plot(x, med, color=color, lw=2, label=label)
    ax.fill_between(x, lo, hi, color=color, alpha=0.15, linewidth=0)


def plot():
    rows = load()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    a, b, c, d = axes.ravel()

    for matcher in MATCHERS:
        rs = [r for r in rows if r["matcher"] == matcher]
        if not rs:
            continue
        col = COLORS[matcher]
        yrs = np.array(rs[0]["years"])
        prev = np.array([r["prev_15_49"] for r in rs]) * 100
        newi = np.array([r["new_inf_15_49"] for r in rs])
        ey = np.array(rs[0]["n_edges_years"])
        ne = np.array([r["n_edges"] for r in rs])
        _band(a, yrs, prev, col, matcher)
        _band(b, yrs, newi, col, matcher)
        _band(c, ey, ne, col, matcher)

    a.set_title("(A) HIV prevalence 15-49"); a.set_ylabel("Prevalence (%)"); a.legend(fontsize=8); a.set_xlim(1985, 2031)
    b.set_title("(B) New infections 15-49 per year"); b.set_ylabel("New infections"); b.set_xlim(1985, 2031)
    c.set_title("(C) Partnership volume: active edges"); c.set_ylabel("Active edges"); c.set_xlim(1985, 2031)

    # (D) age gap by female group, grouped bars
    x = np.arange(len(AGE_GAP_BINS))
    w = 0.25
    for i, matcher in enumerate(MATCHERS):
        rs = [r for r in rows if r["matcher"] == matcher]
        means = []
        for lo, hi in AGE_GAP_BINS:
            vals = [r["age_gaps"].get(f"{lo}_{hi}") for r in rs if r["age_gaps"].get(f"{lo}_{hi}") is not None]
            means.append(np.mean(vals) if vals else np.nan)
        d.bar(x + (i - 1) * w, means, w, color=COLORS[matcher], label=matcher)
    d.plot(x, [DHS_GAP[b_] for b_ in AGE_GAP_BINS], "k*", markersize=12, label="DHS Eswatini")
    d.set_title("(D) Realized male-female age gap (2020)")
    d.set_ylabel("Age gap (years)")
    d.set_xticks(x); d.set_xticklabels([f"{lo}-{hi}" for lo, hi in AGE_GAP_BINS])
    d.legend(fontsize=8)

    fig.suptitle("Exp 013 — matcher comparison (10 seeds, starsim 3.5.0 / stisim 1.5.8)", fontsize=14)
    fig.tight_layout()
    out = FIG / "matcher_comparison.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("Saved", out)

    # Summary table -> outputs/summary_table.csv + stdout
    lines = ["matcher,final_prev_2021_pct,mean_lp_f,total_pairs_f,n_active_2020,gap_15_24,gap_25_34,gap_35_49"]
    yrs = np.array(rows[0]["years"])
    i2021 = int(np.argmin(np.abs(yrs - 2021)))
    print("\n  matcher                        prev2021%  mean_lpF  totalPairsF  nActive2020  gaps(15-24/25-34/35-49)")
    for matcher in MATCHERS:
        rs = [r for r in rows if r["matcher"] == matcher]
        prev = np.mean([r["prev_15_49"][i2021] for r in rs]) * 100
        mlp = np.mean([r["mean_lp_f"] for r in rs])
        tot = np.mean([r["total_pairs_f"] for r in rs])
        nact = np.mean([r["n_active_pairs_2020"] for r in rs])
        g = [np.mean([r["age_gaps"].get(f"{lo}_{hi}") for r in rs
                      if r["age_gaps"].get(f"{lo}_{hi}") is not None]) for lo, hi in AGE_GAP_BINS]
        print(f"  {matcher:30s} {prev:7.1f}   {mlp:7.2f}   {tot:10.0f}   {nact:9.0f}   {g[0]:.2f}/{g[1]:.2f}/{g[2]:.2f}")
        lines.append(f"{matcher},{prev:.2f},{mlp:.3f},{tot:.0f},{nact:.0f},{g[0]:.3f},{g[1]:.3f},{g[2]:.3f}")
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
