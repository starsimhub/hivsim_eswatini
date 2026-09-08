"""Scorecard and figures for exp 026.

Reads outputs/sims/*.parquet and produces the numbers the abstract needs.

Primary outcome: cumulative new infections over the scenario window (2026-2040)
and infections averted against `baseline`. Counts are already at real
population scale -- starsim Results carry scale=True by default, which is why
the capacity check read n_on_art = 229,726 rather than a 10,000-agent count.

The headline quantity is the PrEP INCREMENT at high ART coverage, taken with
the ART ramp held constant: `both` minus `art_95`. See contrasts() -- every
reported difference is named in run.py::CONTRASTS with its comparator, because
the two errors in the first version of this experiment were both differences
between arms that differed in more than one thing.

Uncertainty: seeds only. One parameter set means these intervals are Monte
Carlo error, NOT parameter uncertainty. Reported as a seed range, never as a
credible interval.
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
OUT, FIG = HERE / "outputs", HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)

WINDOW = (2026, 2040)
ORDER = ["baseline", "prep_agyw", "prep_agyw_risk", "prep_fsw",
         "art_95", "both", "art_95_early", "prep_at_high_art"]
LABEL = {"baseline": "Baseline",
         "prep_agyw": "LEN 30% AGYW",
         "prep_agyw_risk": "LEN 30% AGYW (higher-risk)",
         "prep_fsw": "LEN 60% FSW",
         "art_95": "ART coverage to 95%",
         "both": "LEN + ART 95% (concurrent)",
         "art_95_early": "ART 95% by 2025 (early)",
         "prep_at_high_art": "LEN added at ART 95% (early)"}


def load():
    files = sorted(glob.glob(str(OUT / "sims" / "*.parquet")))
    if not files:
        sys.exit("no outputs/sims/*.parquet -- run run.py first")
    d = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    print(f"loaded {len(files)} runs: {d.arm.nunique()} arms x "
          f"{d.seed.nunique()} seeds")
    return d


def _band_cols(d, stem, sex=None):
    """popagesex.<stem>_<sex>_<lo>_<hi> columns, excluding aggregates."""
    out = []
    for c in d.columns:
        if not c.startswith(f"popagesex.{stem}_"):
            continue
        parts = c.split("_")
        if len(parts) < 3:
            continue
        s = parts[-3]
        if s not in ("f", "m"):
            continue
        if sex is not None and s != sex:
            continue
        out.append(c)
    return out


def add_derived(d):
    d = d.copy()
    d["new_inf"] = d[_band_cols(d, "new_infections")].sum(axis=1)
    for s in ("f", "m"):
        d[f"new_inf_{s}"] = d[_band_cols(d, "new_infections", s)].sum(axis=1)
    inf = _band_cols(d, "n_infected")
    alv = [c.replace("n_infected", "n_alive") for c in inf]
    d["prev_all"] = d[inf].sum(axis=1) / d[alv].sum(axis=1)
    for col in ("cascade.p_vls_f", "cascade.p_vls_m", "cascade.p_vls"):
        if col in d.columns:
            d["pvls"] = d[col]
            break
    if "hiv.new_deaths" in d.columns:
        d["deaths"] = d["hiv.new_deaths"]
    return d


def scorecard(d):
    lo, hi = WINDOW
    w = d[(d.timevec >= lo) & (d.timevec <= hi)]
    per_seed = (w.groupby(["arm", "seed"])
                  .agg(cum_inf=("new_inf", "sum"),
                       cum_inf_f=("new_inf_f", "sum"),
                       cum_inf_m=("new_inf_m", "sum"),
                       py_prep=("hiv.n_on_prep", "sum"))
                  .reset_index())
    base = per_seed[per_seed.arm == "baseline"].set_index("seed")

    rows = []
    for arm in [a for a in ORDER if a in set(per_seed.arm)]:
        g = per_seed[per_seed.arm == arm].set_index("seed")
        # Paired by SEED against baseline. Pairing matters: seed-to-seed
        # variation in cumulative infections is larger than several of the arm
        # effects, so an unpaired difference of means would be much noisier
        # than the design actually is.
        av = (base.cum_inf - g.cum_inf).dropna()
        rows.append(dict(
            arm=arm, label=LABEL.get(arm, arm),
            cum_inf_mean=g.cum_inf.mean(), cum_inf_sd=g.cum_inf.std(ddof=1),
            averted_mean=av.mean(), averted_sd=av.std(ddof=1),
            averted_min=av.min(), averted_max=av.max(),
            pct_averted=100 * av.mean() / base.cum_inf.mean(),
            female_share=100 * g.cum_inf_f.mean() / g.cum_inf.mean(),
            py_prep=g.py_prep.mean(),
            # Person-years on PrEP per infection averted -- the efficiency
            # metric, and the only fair way to compare arms that cover
            # different NUMBERS of people. `prep_agyw_risk` is 30% of a
            # smaller group, so it covers 0.39x the person-years of
            # `prep_agyw` and averts less in total while being MORE efficient
            # per person. Comparing totals alone would read that backwards.
            py_per_averted=(g.py_prep.mean() / av.mean()
                            if av.mean() > 0 else np.nan)))
    sc = pd.DataFrame(rows)
    sc.to_csv(OUT / "scorecard.csv", index=False)
    per_seed.to_csv(OUT / "per_seed.csv", index=False)
    return sc, per_seed


def headline(per_seed):
    """What PrEP adds on top of ART at 95%, with ART TIMING HELD CONSTANT.

    CORRECTED. The obvious comparison, prep_at_high_art - art_95, is
    CONFOUNDED: prep_at_high_art was defined with its cascade starting in 2020
    rather than 2026, so it also carries six extra years of ART scale-up. Its
    ART coverage is 0.948 at 2025 against 0.918 in art_95, and the difference
    came out at 8,757 -- of which only part is PrEP.

    `both` has the SAME 2026-2030 ART ramp as `art_95` and adds PrEP on top, so
    `both - art_95` isolates the PrEP increment. That is 5,445, not 8,757.

    prep_at_high_art is still reported, as the "earlier and harder cascade
    push, plus PrEP" arm it actually is -- but it is not the headline.
    """
    have = set(per_seed.arm)
    if not {"both", "art_95"} <= have:
        return None
    a = per_seed[per_seed.arm == "art_95"].set_index("seed").cum_inf
    b = per_seed[per_seed.arm == "both"].set_index("seed").cum_inf
    diff = (a - b).dropna()
    base = per_seed[per_seed.arm == "baseline"].set_index("seed").cum_inf
    res = dict(mean=diff.mean(), sd=diff.std(ddof=1),
               lo=diff.min(), hi=diff.max(),
               pct_of_baseline=100 * diff.mean() / base.mean(),
               pct_of_remaining=100 * diff.mean() / a.mean(),
               n_seeds=len(diff),
               n_seeds_positive=int((diff > 0).sum()))
    pd.Series(res).to_csv(OUT / "headline.csv")
    return res


def contrasts(per_seed):
    """Every named contrast from run.py::CONTRASTS, paired by seed.

    Reported from a named dict rather than assembled ad hoc, because the two
    errors in the first version of this experiment were both differences taken
    between arms that differed in more than one thing. Naming each contrast and
    its comparator makes that visible instead of implicit.
    """
    sys.path.insert(0, str(HERE))
    import run as R

    rows = []
    for name, (treat, comp) in R.CONTRASTS.items():
        have = set(per_seed.arm)
        if not {treat, comp} <= have:
            rows.append(dict(contrast=name, treat=treat, comparator=comp,
                             note="MISSING ARM"))
            continue
        t = per_seed[per_seed.arm == treat].set_index("seed").cum_inf
        c = per_seed[per_seed.arm == comp].set_index("seed").cum_inf
        diff = (c - t).dropna()
        rows.append(dict(
            contrast=name, treat=treat, comparator=comp,
            averted=diff.mean(), sd=diff.std(ddof=1),
            lo=diff.min(), hi=diff.max(),
            pct_of_comparator=100 * diff.mean() / c.mean(),
            seeds_positive=f"{int((diff > 0).sum())}/{len(diff)}"))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "contrasts.csv", index=False)
    return t


def plots(d, sc, per_seed):
    lo, hi = WINDOW
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8),
                             gridspec_kw=dict(width_ratios=[1.15, 1.35, 1]))

    # 1. infections averted, with the seed spread shown as points
    ax = axes[0]
    s = sc[sc.arm != "baseline"]
    y = np.arange(len(s))
    ax.barh(y, s.averted_mean, color="#4a7fb5", zorder=2)
    for i, arm in enumerate(s.arm):
        base = per_seed[per_seed.arm == "baseline"].set_index("seed").cum_inf
        g = per_seed[per_seed.arm == arm].set_index("seed").cum_inf
        ax.plot((base - g).dropna(), np.full(len(g), i), "o", ms=3,
                color="k", alpha=0.45, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(s.label, fontsize=8)
    ax.invert_yaxis(); ax.axvline(0, c="k", lw=1)
    ax.set_xlabel(f"infections averted, {lo}-{hi}")
    ax.set_title("Infections averted vs baseline\n(dots = individual seeds)",
                 fontsize=10)
    ax.grid(axis="x", alpha=0.3)

    # 2. incidence trajectories
    ax = axes[1]
    inf = _band_cols(d, "n_infected")
    alv = [c.replace("n_infected", "n_alive") for c in inf]
    d2 = d.copy()
    d2["susc"] = d2[alv].sum(axis=1) - d2[inf].sum(axis=1)
    d2["inc"] = 100 * d2.new_inf / d2.susc
    cmap = plt.get_cmap("tab10")
    for i, arm in enumerate([a for a in ORDER if a in set(d2.arm)]):
        g = d2[d2.arm == arm].groupby("timevec").inc
        m, sd = g.mean(), g.std(ddof=1).fillna(0)
        ax.plot(m.index, m.values, lw=2, color=cmap(i), label=LABEL.get(arm, arm))
        ax.fill_between(m.index, m - sd, m + sd, color=cmap(i), alpha=0.12, lw=0)
    ax.axvline(lo, ls=":", c="grey")
    ax.set_xlim(2015, hi); ax.set_ylim(bottom=0)
    ax.set_xlabel("year"); ax.set_ylabel("HIV incidence (% per year)")
    ax.set_title("Incidence by arm (mean ± 1 SD over seeds)", fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)

    # 3. the headline subtraction
    ax = axes[2]
    have = set(per_seed.arm)
    if {"prep_at_high_art", "art_95"} <= have:
        a = per_seed[per_seed.arm == "art_95"].set_index("seed").cum_inf
        b = per_seed[per_seed.arm == "prep_at_high_art"].set_index("seed").cum_inf
        diff = (a - b).dropna()
        ax.axhline(0, c="k", lw=1)
        ax.plot(np.zeros(len(diff)) + np.random.uniform(-.05, .05, len(diff)),
                diff, "o", ms=7, color="#c0392b", alpha=0.8)
        ax.plot([-.18, .18], [diff.mean()] * 2, "-", lw=3, color="#c0392b")
        ax.set_xlim(-.5, .5); ax.set_xticks([])
        ax.set_ylabel(f"infections averted, {lo}-{hi}")
        ax.set_title("What LEN adds AFTER ART reaches 95%\n"
                     "(each dot one seed; bar = mean)", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "scenarios.png", dpi=130)
    plt.close(fig)


def main():
    d = add_derived(load())
    sc, per_seed = scorecard(d)
    pd.set_option("display.width", 200)
    print(f"\n=== Cumulative new infections {WINDOW[0]}-{WINDOW[1]} "
          f"(real population scale) ===")
    print(sc[["label", "cum_inf_mean", "cum_inf_sd", "averted_mean",
              "averted_sd", "pct_averted", "female_share"]]
          .round(1).to_string(index=False))

    eff = sc[(sc.py_prep > 0)]
    if len(eff):
        print("\n=== PrEP efficiency (person-years on PrEP per infection averted) ===")
        print("    The fair comparison between arms that cover different numbers")
        print("    of people. Lower is better.")
        print(eff[["label", "py_prep", "averted_mean", "py_per_averted"]]
              .round(1).to_string(index=False))

    ct = contrasts(per_seed)
    print("\n=== Every named contrast, paired by seed ===")
    print("    'averted' = comparator minus treat, so positive means the treat")
    print('    arm has FEWER infections.')
    print(ct.round(1).to_string(index=False))

    h = headline(per_seed)
    if h:
        print("\n=== HEADLINE: what LEN adds ON TOP of ART 95% ===")
        print(f"  infections averted: {h['mean']:,.0f} "
              f"(seed range {h['lo']:,.0f} to {h['hi']:,.0f}, "
              f"SD {h['sd']:,.0f}, n={h['n_seeds']})")
        print(f"  = {h['pct_of_baseline']:.1f}% of baseline infections")
        print(f"  = {h['pct_of_remaining']:.1f}% of those REMAINING once "
              f"ART is at 95%")
        print(f"  positive in {h['n_seeds_positive']}/{h['n_seeds']} seeds")
        conf = per_seed[per_seed.arm == "prep_at_high_art"]
        if len(conf):
            a = per_seed[per_seed.arm == "art_95"].set_index("seed").cum_inf
            c = conf.set_index("seed").cum_inf
            print(f"\n  For contrast, prep_at_high_art - art_95 = "
                  f"{(a - c).mean():,.0f} -- but that arm's cascade starts in "
                  f"2020, so it carries six extra years of ART scale-up "
                  f"(coverage 0.948 vs 0.918 at 2025) and is NOT the PrEP "
                  f"increment.")
        if h["n_seeds_positive"] < h["n_seeds"] * 0.9:
            print("  NOTE: not consistently positive across seeds -- the effect "
                  "is within Monte Carlo noise at this seed count.")
    plots(d, sc, per_seed)
    print(f"\nwrote outputs/scorecard.csv, per_seed.csv, headline.csv")
    print(f"      figures/scenarios.png")


if __name__ == "__main__":
    main()
