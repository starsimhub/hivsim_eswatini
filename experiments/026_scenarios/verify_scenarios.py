"""What each arm actually does, stated explicitly and then verified from output.

Two jobs:

1. **Spec.** Print, for every arm, exactly what changes relative to baseline --
   the intervention added, the coverage table overridden, the ramp years. Read
   off ARMS/ARM_OVERRIDES in run.py rather than retyped, so the description
   cannot drift from what ran.

2. **Verification.** Plot the REALISED coverage from the simulation output, not
   the requested coverage. These are different things and the gap is where bugs
   live: a coverage table that stops before the sim end forward-fills, an
   eligibility mask that is smaller than intended, a stock target that never
   reaches its value because the pool is too small. Two arm-definition errors
   in this experiment were caught this way rather than by anything failing.

Usage (from the repo ROOT)
  python experiments/026_scenarios/verify_scenarios.py
"""

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import run as R          # noqa: E402  -- ARMS, ARM_OVERRIDES, CONTRASTS
import scenarios as S    # noqa: E402

OUT, FIG = HERE / "outputs", HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)

ORDER = ["baseline", "prep_agyw", "prep_agyw_risk", "prep_fsw",
         "art_95", "both", "art_95_early", "prep_at_high_art"]


def spec_table():
    """What each arm changes vs baseline, derived from run.py's own definitions."""
    art_tbl, vls_tbl = S.baseline_tables()
    base_art = art_tbl[art_tbl.Year == art_tbl.Year.max()].p_art
    rows = []
    for arm in ORDER:
        if arm not in R.ARMS:
            continue
        prep_spec, casc = R.ARMS[arm]
        ov = R.ARM_OVERRIDES.get(arm, {})
        prep = "-"
        if prep_spec is not None:
            cov, elig = prep_spec
            prep = (f"LEN {cov:.0%} of {elig}, "
                    f"{R.SCEN_START}-{R.SCEN_REACH}, "
                    f"eff {S.LEN_EFF}, q{S.LEN_DUR_MONTHS}mo")
        cascade = "-"
        if casc is not None:
            art_t, _ = casc
            s = ov.get("cascade_start", R.SCEN_START)
            r = ov.get("cascade_reach", R.SCEN_REACH)
            cascade = (f"ART coverage -> {art_t:.0%} by {r} "
                       f"(ramp from {s}; baseline "
                       f"{base_art.min():.2f}-{base_art.max():.2f})")
        rows.append(dict(arm=arm, prep_vs_baseline=prep,
                         cascade_vs_baseline=cascade))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "scenario_spec.csv", index=False)
    return t


def load():
    files = sorted(glob.glob(str(OUT / "sims" / "*.parquet")))
    if not files:
        sys.exit("no outputs/sims/*.parquet")
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


def derive(d):
    d = d.copy()
    ni = [c for c in d.columns
          if c.startswith("popagesex.new_infections_")
          and c.split("_")[-3] in ("f", "m")]
    inf = [c for c in d.columns
           if c.startswith("popagesex.n_infected_")
           and c.split("_")[-3] in ("f", "m")]
    alv = [c.replace("n_infected", "n_alive") for c in inf]
    d["new_inf"] = d[ni].sum(axis=1)
    d["n_inf"] = d[inf].sum(axis=1)
    d["n_alive"] = d[alv].sum(axis=1)
    d["susc"] = d.n_alive - d.n_inf
    d["incidence"] = 100 * d.new_inf / d.susc.replace(0, np.nan)

    # PrEP coverage among HIV-negative AGYW -- the denominator the AGYW arms
    # actually target. Reported so a "30% coverage" arm can be checked against
    # the 30% it was asked for, rather than assumed.
    agyw_a = [f"popagesex.n_alive_f_{lo}_{lo+5}" for lo in (15, 20)]
    agyw_i = [f"popagesex.n_infected_f_{lo}_{lo+5}" for lo in (15, 20)]
    have = [c for c in agyw_a + agyw_i if c in d.columns]
    if len(have) == 4:
        d["agyw_neg"] = d[agyw_a].sum(axis=1) - d[agyw_i].sum(axis=1)
        d["prep_cov_agyw"] = d["hiv.n_on_prep"] / d.agyw_neg.replace(0, np.nan)
    return d


def verify_table(d):
    """Realised values at key years -- the numbers behind the figure."""
    rows = []
    for arm in [a for a in ORDER if a in set(d.arm)]:
        g = d[d.arm == arm].groupby("timevec")
        r = {"arm": arm}
        for y in (2021, 2025, 2030, 2040):
            if y not in g.groups:
                continue
            sub = d[(d.arm == arm) & (d.timevec == y)]
            r[f"art_{y}"] = round(float(sub["hiv.p_on_art"].mean()), 3)
            r[f"prep_n_{y}"] = int(sub["hiv.n_on_prep"].mean())
        rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "scenario_verification.csv", index=False)
    return t


def plots(d):
    arms = [a for a in ORDER if a in set(d.arm)]
    cmap = plt.get_cmap("tab10")
    col = {a: cmap(i) for i, a in enumerate(arms)}
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    def band(ax, series_name, ylab, title, pct=False):
        for a in arms:
            g = d[d.arm == a].groupby("timevec")[series_name]
            m, sd = g.mean(), g.std(ddof=1).fillna(0)
            k = 100 if pct else 1
            ax.plot(m.index, m.values * k, lw=2, color=col[a], label=a)
            ax.fill_between(m.index, (m - sd) * k, (m + sd) * k,
                            color=col[a], alpha=0.10, lw=0)
        ax.axvline(R.SCEN_START, ls=":", c="grey", lw=1.2)
        ax.set_xlim(2015, 2040); ax.set_xlabel("year"); ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=10); ax.grid(alpha=0.3)

    ax = axes[0, 0]
    band(ax, "hiv.p_on_art", "ART coverage (of PLHIV)",
         "ART coverage -- cascade arms should separate here", pct=True)
    ax.axhline(95, ls="--", c="k", lw=1)
    ax.annotate("95% target", (2016, 95.5), fontsize=7, color="k")
    ax.legend(fontsize=7, loc="lower right")

    ax = axes[0, 1]
    band(ax, "hiv.n_on_prep", "people on PrEP",
         "PrEP uptake -- non-PrEP arms must sit flat at zero")
    ax.legend(fontsize=7, loc="upper left")

    ax = axes[1, 0]
    if "prep_cov_agyw" in d.columns:
        band(ax, "prep_cov_agyw", "% of HIV-negative AGYW on PrEP",
             "Realised AGYW coverage vs the 30% requested", pct=True)
        ax.axhline(30, ls="--", c="k", lw=1)
        ax.annotate("30% requested", (2016, 31), fontsize=7)
        ax.legend(fontsize=7, loc="upper left")

    ax = axes[1, 1]
    band(ax, "incidence", "HIV incidence (% per year)",
         "Incidence -- the outcome the contrasts are built from")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=7, loc="upper right")

    fig.suptitle("Exp 026 scenario verification -- REALISED coverage from the "
                 "simulation, not requested coverage\n"
                 "dotted line = scenario start (2026); bands = 1 SD over 10 seeds",
                 y=1.0, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "scenario_verification.png", dpi=130,
                bbox_inches="tight")
    plt.close(fig)


def main():
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 70)
    print("=== What each arm changes, relative to baseline ===")
    print(spec_table().to_string(index=False))

    d = derive(load())
    print(f"\n=== Realised values ({d.arm.nunique()} arms x "
          f"{d.seed.nunique()} seeds) ===")
    print(verify_table(d).to_string(index=False))

    print("\n=== Checks ===")
    ok = True
    for a in [x for x in ORDER if x in set(d.arm)]:
        sub = d[d.arm == a]
        has_prep = R.ARMS[a][0] is not None
        pk = float(sub["hiv.n_on_prep"].max())
        if has_prep and pk <= 0:
            print(f"  FAIL {a}: PrEP arm but nobody ever on PrEP"); ok = False
        if not has_prep and pk > 0:
            print(f"  FAIL {a}: no-PrEP arm but {pk:.0f} on PrEP"); ok = False
        pre = sub[sub.timevec < R.SCEN_START]
        if has_prep and float(pre["hiv.n_on_prep"].max()) > 0:
            print(f"  FAIL {a}: PrEP before {R.SCEN_START}"); ok = False
    # arms must be identical before the earliest scenario change
    first = min([R.ARM_OVERRIDES.get(a, {}).get("cascade_start", R.SCEN_START)
                 for a in R.ARMS])
    pre = d[d.timevec < first].groupby(["arm", "seed"]).new_inf.sum().unstack(0)
    spread = (pre.max(axis=1) - pre.min(axis=1)).max()
    print(f"  pre-{first} infections identical across arms: "
          f"max spread {spread:.0f} {'OK' if spread == 0 else 'MISMATCH'}")
    if spread != 0:
        ok = False
    print("  all checks passed" if ok else "  SOME CHECKS FAILED")

    plots(d)
    print("\nwrote outputs/scenario_spec.csv, scenario_verification.csv")
    print("      figures/scenario_verification.png")


if __name__ == "__main__":
    main()
