"""Scorecard, cascade reachability audit and figures for exp 031.

026's analyse.py plus the measurement it never took. Three sections are new:

`gate()`      -- baseline cumulative infections against 026's 63,162. v1.4 and
                 v1.5 are recorded in their tags as near-inert, so the testing
                 split should be the ONLY difference between 026 and 031. If it
                 is not, no contrast below is trustworthy and the run stops
                 here.

`reachability()` -- per arm, stratum and year: the ART target the arm actually
                 asked for, the awareness ceiling, and the coverage delivered.
                 029 established that awareness is a hard ceiling (stisim fills
                 a stratified target only from `diagnosed & ~on_art`) and that
                 art_95 sat above it in three strata; 026 never checked, so its
                 headline number is for a smaller intervention than its label.
                 This is the section the experiment exists for.

`prevalence()` -- the standard figure, from the baseline arm.

Everything else -- scorecard, contrasts, headline, scenarios.png -- is 026's,
copied rather than re-derived so the two are directly comparable.

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
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))

from cascade_analysis import aggregate, conditionals, ART_BINS  # noqa: E402
import standard_figures as SF                                   # noqa: E402

OUT, FIG = HERE / "outputs", HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)

WINDOW = (2026, 2040)

# 026's headline numbers, for the gate and for the SUMMARY's comparison column.
# Source: experiments/026_scenarios/SUMMARY.md. Recorded here as literals so
# this script fails loudly if the new run disagrees, rather than quietly
# reporting a different number under the same label.
REF_026 = dict(baseline_cum_inf=63162, art_95_averted=16525,
               prep_fsw_averted=8771, prep_increment=4449)

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

BINLAB = {"[15,25)": "15-24", "[25,35)": "25-34",
          "[35,45)": "35-44", "[45,100)": "45+"}
SEXLAB = {"f": "Women", "m": "Men"}


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
        # Paired by SEED against baseline. Seed-to-seed variation in cumulative
        # infections is larger than several of the arm effects, so an unpaired
        # difference of means would be much noisier than the design really is.
        av = (base.cum_inf - g.cum_inf).dropna()
        rows.append(dict(
            arm=arm, label=LABEL.get(arm, arm),
            cum_inf_mean=g.cum_inf.mean(), cum_inf_sd=g.cum_inf.std(ddof=1),
            averted_mean=av.mean(), averted_sd=av.std(ddof=1),
            averted_min=av.min(), averted_max=av.max(),
            pct_averted=100 * av.mean() / base.cum_inf.mean(),
            female_share=100 * g.cum_inf_f.mean() / g.cum_inf.mean(),
            py_prep=g.py_prep.mean(),
            py_per_averted=(g.py_prep.mean() / av.mean()
                            if av.mean() > 0 else np.nan)))
    sc = pd.DataFrame(rows)
    sc.to_csv(OUT / "scorecard.csv", index=False)
    per_seed.to_csv(OUT / "per_seed.csv", index=False)
    return sc, per_seed


def gate(per_seed):
    """Does the baseline still reproduce 026's 63,162?

    026 ran at model-v1.3. Between then and v1.6 sit v1.4 (acute-phase CD4
    floor) and v1.5 (non-negative CD4 invariant), both recorded in their tags as
    near-inert -- v1.4's tag states an identical seeded run gives cum_infections
    82943 either way -- plus the testing split, which IS meant to change things.

    The testing split changes who is DIAGNOSED. Under baseline it should barely
    touch who is INFECTED, because ART coverage is forced to its data target
    either way; 030 obs 6 found baseline coverage survives in 215 of 216
    strata-years. So the baseline arm is where "nothing else moved" gets
    checked. A large move here means something other than testing changed, and
    every contrast below inherits the problem.
    """
    b = per_seed[per_seed.arm == "baseline"].cum_inf
    ref = REF_026["baseline_cum_inf"]
    z = (b.mean() - ref) / b.std(ddof=1) if b.std(ddof=1) > 0 else np.nan
    ok = bool(b.min() <= ref <= b.max())
    res = dict(mean=b.mean(), sd=b.std(ddof=1), lo=b.min(), hi=b.max(),
               ref_026=ref, diff=b.mean() - ref,
               pct_diff=100 * (b.mean() - ref) / ref,
               sd_units=z, ref_inside_seed_range=ok, n_seeds=len(b))
    pd.Series(res).to_csv(OUT / "gate.csv")
    return res


def reachability(d):
    """Target vs awareness ceiling vs achieved, per arm, stratum and year.

    The primary output. Format follows 030's, extended across arms.

    Each arm asks for a DIFFERENT ART target -- that is what the cascade arms
    are -- so the target series is rebuilt per arm from run.make_arm_kwargs(),
    the same call the run itself used. Reading one shared target table would
    compare every arm against the baseline's target and silently report the
    cascade arms as enormously over-delivering.
    """
    import run as R

    art_cond = conditionals(
        aggregate(d, ART_BINS, group_cols=("seed", "arm")),
        keys=("arm", "year", "sex", "bin"))
    piv = (art_cond[art_cond.measure.isin(["aware_of_all_plhiv",
                                           "on_art_of_all_plhiv"])]
           .pivot_table(index=["arm", "year", "sex", "bin"],
                        columns="measure", values="model").reset_index())

    years = np.sort(piv.year.unique())
    out = []
    for arm in sorted(piv.arm.unique()):
        tgt = R.make_arm_kwargs(arm).get("art_coverage")
        if tgt is None:                       # no cascade spec -> the baseline table
            tgt = pd.read_csv(REPO / "data" / "art_coverage.csv")
        tgt = tgt.copy()
        # Gender 0 = FEMALE here (stisim's convention, per
        # art_coverage_construction.py). This is the OPPOSITE of
        # calibration_data/art_coverage_by_age_sex.csv. Reversing it silently
        # swaps the sexes and still produces a plausible figure (029 obs 8).
        tgt["sex"] = tgt.Gender.map({0: "f", 1: "m"})
        tgt = tgt.rename(columns={"Year": "year", "AgeBin": "bin",
                                  "p_art": "target"})
        # The table is tabulated at a handful of years and stisim INTERPOLATES
        # between them, so merging on tabulated years alone would check a few
        # dozen strata-years instead of every year the model ran -- and the
        # binding constraint can fall in an interpolated year. np.interp holds
        # flat outside the tabulated range, which is stisim's forward-fill.
        for (sex, b), g in tgt.groupby(["sex", "bin"]):
            g = g.sort_values("year")
            out.append(pd.DataFrame({
                "arm": arm, "year": years, "sex": sex, "bin": b,
                "target": np.interp(years, g.year.values, g.target.values),
                "interpolated": ~np.isin(years, g.year.values)}))
    tgt_annual = pd.concat(out, ignore_index=True)

    reach = piv.merge(tgt_annual, on=["arm", "year", "sex", "bin"], how="inner")
    reach["shortfall_pp"] = (reach.on_art_of_all_plhiv - reach.target) * 100
    reach["ceiling_slack_pp"] = (reach.aware_of_all_plhiv - reach.target) * 100
    reach = reach.sort_values(["arm", "sex", "bin", "year"])
    reach.to_csv(OUT / "art_reachability.csv", index=False)
    return reach


def delivered(reach, arm="art_95", year=2030):
    """What the cascade arm actually delivered, per stratum -- the answer.

    Separates the two ways an arm can fall short, which 026 could not
    distinguish and which have different implications:

    `ceiling_slack_pp` < 0  -- the target sits ABOVE awareness. Structurally
                               unreachable; no amount of linkage effort fills
                               it. This is 029 obs 2.
    `shortfall_pp`     < 0  -- the target was not met. If the ceiling had slack,
                               this is ordinary under-delivery; if it did not,
                               the ceiling is the explanation.
    """
    g = reach[(reach.arm == arm) & (reach.year == year)].copy()
    g["label"] = g.sex.map(SEXLAB) + " " + g["bin"].map(BINLAB).fillna(g["bin"])
    g["reason"] = np.where(g.ceiling_slack_pp < 0, "ceiling",
                           np.where(g.shortfall_pp < -0.5, "under-delivered",
                                    "met"))
    cols = ["label", "target", "aware_of_all_plhiv", "on_art_of_all_plhiv",
            "ceiling_slack_pp", "shortfall_pp", "reason"]
    g[cols].to_csv(OUT / f"delivered_{arm}_{year}.csv", index=False)
    return g[cols]


def cascade_decomposition(d, arm="baseline", years=(2021, 2030, 2040)):
    """Is Eswatini at 95-95-95, and what does the aggregate hide?

    Added after the main analysis, from the same outputs. The UNAIDS target is
    a POPULATION AVERAGE, so it can be met while individual strata sit far
    below it -- and whether that happens depends entirely on where the PLHIV
    stock sits. In Eswatini the stock is concentrated in older women, who have
    had decades to be diagnosed and treated, so they carry the average.

    Writes two things: the aggregate product against the 0.857 target
    (0.95^3), and the per-stratum decomposition weighted by PLHIV share.
    """
    from cascade_analysis import SHIMS_BINS

    sub = d[d.arm == arm]
    rows = []
    agg = conditionals(aggregate(sub, {"15+": (15, 200)}, group_cols=("seed",)),
                       keys=("year", "bin"))
    for y in years:
        p = agg[agg.year == y].set_index("measure").model
        a, l, v = (p["aware_of_all_plhiv"], p["on_art_given_aware"],
                   p["vls_given_art"])
        rows.append(dict(year=y, aware=a, on_art_given_aware=l, vls_given_art=v,
                         suppressed_of_plhiv=a * l * v, target=0.95 ** 3,
                         meets_95_95_95=bool(a * l * v >= 0.95 ** 3)))
    ag = pd.DataFrame(rows)
    ag.to_csv(OUT / "cascade_aggregate.csv", index=False)

    # Per stratum, with the PLHIV weight that explains the aggregate.
    lg = aggregate(sub, SHIMS_BINS, group_cols=("seed",))
    y = years[1]
    c = conditionals(lg, keys=("year", "sex", "bin"))
    p = (c[(c.year == y) & c.measure.isin(["aware_of_all_plhiv",
                                           "on_art_given_aware",
                                           "vls_given_art"])]
         .pivot_table(index=["sex", "bin"], columns="measure",
                      values="model").reset_index())
    w = lg[lg.year == y].groupby(["sex", "bin"], as_index=False).n_infected.mean()
    m = p.merge(w, on=["sex", "bin"])
    m = m[m["bin"].isin(["15-24", "25-34", "35-49", "50+"])].copy()
    m["share_of_plhiv"] = 100 * m.n_infected / m.n_infected.sum()
    m["suppressed_of_plhiv"] = (m.aware_of_all_plhiv * m.on_art_given_aware
                                * m.vls_given_art)
    m["year"] = y
    m = m.sort_values("suppressed_of_plhiv", ascending=False)
    m.to_csv(OUT / "cascade_by_stratum.csv", index=False)
    return ag, m


def contrasts(per_seed):
    """Every named contrast from run.py::CONTRASTS, paired by seed.

    From a named dict rather than assembled ad hoc, because the two errors in
    the first version of 026 were both differences taken between arms that
    differed in more than one thing.
    """
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


def headline(per_seed):
    """What PrEP adds on top of ART at 95%, with ART TIMING HELD CONSTANT.

    `both` has the same 2026-2030 ART ramp as `art_95` and adds PrEP on top, so
    `both - art_95` isolates the PrEP increment. `prep_at_high_art - art_95` is
    confounded -- that arm's cascade starts in 2020, so it carries six extra
    years of ART scale-up. See 026's SUMMARY.
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
               n_seeds=len(diff), n_seeds_positive=int((diff > 0).sum()),
               ref_026=REF_026["prep_increment"])
    pd.Series(res).to_csv(OUT / "headline.csv")
    return res


def fig_scenarios(d, sc, per_seed):
    lo, hi = WINDOW
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8),
                             gridspec_kw=dict(width_ratios=[1.15, 1.35, 1]))

    ax = axes[0]
    s = sc[sc.arm != "baseline"]
    y = np.arange(len(s))
    ax.barh(y, s.averted_mean, color="#4a7fb5", zorder=2)
    base = per_seed[per_seed.arm == "baseline"].set_index("seed").cum_inf
    for i, arm in enumerate(s.arm):
        g = per_seed[per_seed.arm == arm].set_index("seed").cum_inf
        ax.plot((base - g).dropna(), np.full(len(g), i), "o", ms=3,
                color="k", alpha=0.45, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(s.label, fontsize=8)
    ax.invert_yaxis(); ax.axvline(0, c="k", lw=1)
    ax.set_xlabel(f"infections averted, {lo}-{hi}")
    ax.set_title("Infections averted vs baseline\n(dots = individual seeds)",
                 fontsize=10)
    ax.grid(axis="x", alpha=0.3)

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
    ax.set_title("Incidence by arm (mean +/- 1 SD over seeds)", fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)

    ax = axes[2]
    if {"both", "art_95"} <= set(per_seed.arm):
        a = per_seed[per_seed.arm == "art_95"].set_index("seed").cum_inf
        b = per_seed[per_seed.arm == "both"].set_index("seed").cum_inf
        diff = (a - b).dropna()
        rng = np.random.default_rng(0)
        ax.axhline(0, c="k", lw=1)
        ax.plot(rng.uniform(-.05, .05, len(diff)), diff, "o", ms=7,
                color="#c0392b", alpha=0.8)
        ax.plot([-.18, .18], [diff.mean()] * 2, "-", lw=3, color="#c0392b")
        ax.set_xlim(-.5, .5); ax.set_xticks([])
        ax.set_ylabel(f"infections averted, {lo}-{hi}")
        ax.set_title("What LEN adds AFTER ART reaches 95%\n"
                     "(each dot one seed; bar = mean)", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "scenarios.png", dpi=130)
    plt.close(fig)


def fig_target_vs_achieved(reach, arm="art_95"):
    """The figure the experiment exists for: did the cascade arm get there?"""
    g = reach[reach.arm == arm]
    bins = [b for b in ART_BINS if b in set(g["bin"])]
    fig, axes = plt.subplots(2, len(bins), figsize=(3.5 * len(bins), 6.6),
                             sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    for r, sex in enumerate(("f", "m")):
        for c, b in enumerate(bins):
            ax = axes[r, c]
            s = g[(g.sex == sex) & (g["bin"] == b)].sort_values("year")
            ax.fill_between(s.year, s.target, s.aware_of_all_plhiv,
                            where=s.aware_of_all_plhiv >= s.target,
                            color="#7fbf7b", alpha=0.30, lw=0,
                            label="headroom under the ceiling")
            ax.fill_between(s.year, s.target, s.aware_of_all_plhiv,
                            where=s.aware_of_all_plhiv < s.target,
                            color="#c0392b", alpha=0.30, lw=0,
                            label="target above the ceiling")
            ax.plot(s.year, s.aware_of_all_plhiv, color="#2c6fbb", lw=1.6,
                    label="awareness (the ceiling)")
            ax.plot(s.year, s.target, color="k", lw=1.4, ls="--",
                    label="ART target asked for")
            ax.plot(s.year, s.on_art_of_all_plhiv, color="#c0392b", lw=2.0,
                    label="ART coverage delivered")
            ax.set_xlim(2005, s.year.max())
            ax.set_ylim(0, 1.02)
            ax.grid(alpha=0.25)
            if r == 0:
                ax.set_title(BINLAB.get(b, b), fontsize=10)
            if c == 0:
                ax.set_ylabel(f"{SEXLAB[sex]}\nshare of PLHIV", fontsize=9)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.04))
    fig.suptitle(f"Did the '{LABEL.get(arm, arm)}' arm get where it was asked "
                 f"to go?", fontsize=12)
    fig.text(0.5, -0.085,
             "Awareness is a hard ceiling: stisim fills a stratified ART target "
             "only from agents already diagnosed. Where the dashed target rises "
             "above the blue\nceiling (red shading) the target cannot fill, with "
             "no warning and no error -- the arm silently delivers less than its "
             "label claims.",
             ha="center", fontsize=8, color="#6b6b6b")
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.savefig(FIG / "art_target_vs_achieved.png", dpi=140,
                bbox_inches="tight")
    plt.close(fig)


def prevalence(d):
    """The standard figure, from the baseline arm.

    Required of every experiment (CLAUDE.md). 030 skipped it on the argument
    that testing changes who is diagnosed, not who is infected. That argument
    does not carry here: the arms change ART coverage, and v1.6 changes the
    ceiling on it, so prevalence can move.
    """
    b = d[d.arm == "baseline"].copy()
    try:
        return SF.plot_prevalence_fit(b, "031 baseline (model-v1.6)",
                                      FIG / "prevalence_fit_vs_phia.png")
    except Exception as e:                       # noqa: BLE001
        print(f"  ! standard prevalence figure failed: {e}")
        return None


def main():
    d = add_derived(load())
    pd.set_option("display.width", 220)

    # --- the gate -----------------------------------------------------------
    sc, per_seed = scorecard(d)
    g = gate(per_seed)
    print("\n=== GATE: does the baseline still reproduce 026? ===")
    print(f"  031 baseline cumulative infections {WINDOW[0]}-{WINDOW[1]}: "
          f"{g['mean']:,.0f} (seed range {g['lo']:,.0f} to {g['hi']:,.0f}, "
          f"SD {g['sd']:,.0f}, n={g['n_seeds']})")
    print(f"  026 reported:                        {g['ref_026']:,.0f}")
    print(f"  difference: {g['diff']:+,.0f} ({g['pct_diff']:+.1f}%, "
          f"{g['sd_units']:+.1f} SD)")
    print(f"  026's value inside 031's seed range: "
          f"{'YES' if g['ref_inside_seed_range'] else 'NO'}")
    if not g["ref_inside_seed_range"]:
        print("  -> Read the contrasts below with care. The tags record v1.4 "
              "and v1.5 as near-inert,\n     so a move here needs explaining "
              "before any number goes in the abstract.")

    print(f"\n=== Cumulative new infections {WINDOW[0]}-{WINDOW[1]} ===")
    print(sc[["label", "cum_inf_mean", "cum_inf_sd", "averted_mean",
              "averted_sd", "pct_averted", "female_share"]]
          .round(1).to_string(index=False))

    eff = sc[sc.py_prep > 0]
    if len(eff):
        print("\n=== PrEP efficiency (person-years on PrEP per infection "
              "averted; lower is better) ===")
        print(eff[["label", "py_prep", "averted_mean", "py_per_averted"]]
              .round(1).to_string(index=False))

    ct = contrasts(per_seed)
    print("\n=== Every named contrast, paired by seed ===")
    print("    'averted' = comparator minus treat, so positive means the treat "
          "arm has FEWER infections.")
    print(ct.round(1).to_string(index=False))

    h = headline(per_seed)
    if h:
        print("\n=== HEADLINE: what LEN adds ON TOP of ART 95% ===")
        print(f"  infections averted: {h['mean']:,.0f} "
              f"(seed range {h['lo']:,.0f} to {h['hi']:,.0f}, "
              f"SD {h['sd']:,.0f}, n={h['n_seeds']})")
        print(f"  = {h['pct_of_remaining']:.1f}% of those REMAINING once ART "
              f"is at 95%")
        print(f"  positive in {h['n_seeds_positive']}/{h['n_seeds']} seeds")
        print(f"  026 reported {h['ref_026']:,.0f}")

    # --- the reason the experiment exists -----------------------------------
    reach = reachability(d)
    print("\n=== ART reachability: target vs ceiling vs delivered ===")
    for arm in [a for a in ORDER if a in set(reach.arm)]:
        s = reach[(reach.arm == arm) & (reach.year.between(2026, 2040))]
        if not len(s):
            continue
        worst = s.loc[s.shortfall_pp.idxmin()]
        tight = s.loc[s.ceiling_slack_pp.idxmin()]
        print(f"  {LABEL.get(arm, arm):32} worst shortfall "
              f"{worst.shortfall_pp:+6.2f} pp ({worst.sex} {worst['bin']} "
              f"{int(worst.year)}) | tightest ceiling "
              f"{tight.ceiling_slack_pp:+6.2f} pp | unreachable strata-years "
              f"{int((s.ceiling_slack_pp < 0).sum())}/{len(s)}")

    if "art_95" in set(reach.arm):
        print("\n=== What 'ART coverage to 95%' actually delivered, by "
              "stratum, 2030 ===")
        print(delivered(reach, "art_95", 2030).round(3).to_string(index=False))

    # --- is Eswatini already at 95-95-95, and what does the average hide? ---
    ag, strat = cascade_decomposition(d)
    print("\n=== Aggregate cascade, baseline, both sexes 15+ ===")
    print("    UNAIDS 95-95-95 is a population AVERAGE; the product target is "
          "0.857.")
    print(ag.round(3).to_string(index=False))
    print(f"\n=== What the average hides: by stratum at {int(strat.year.iloc[0])}"
          f" ===")
    print(strat[["sex", "bin", "share_of_plhiv", "aware_of_all_plhiv",
                 "on_art_given_aware", "vls_given_art", "suppressed_of_plhiv"]]
          .round(3).to_string(index=False))
    lag = strat[strat.suppressed_of_plhiv < 0.857]
    print(f"  strata below the 95-95-95 product: {len(lag)} of {len(strat)}, "
          f"holding {lag.share_of_plhiv.sum():.1f}% of PLHIV")
    print(f"  spread: {strat.suppressed_of_plhiv.min():.3f} to "
          f"{strat.suppressed_of_plhiv.max():.3f}")

    fig_scenarios(d, sc, per_seed)
    if "art_95" in set(reach.arm):
        fig_target_vs_achieved(reach, "art_95")
    pv = prevalence(d)
    if pv:
        print(f"\nstandard prevalence fit: {pv}")

    print("\nwrote outputs/{gate,scorecard,per_seed,contrasts,headline,"
          "art_reachability,delivered_art_95_2030}.csv")
    print("      figures/{scenarios,art_target_vs_achieved,"
          "prevalence_fit_vs_phia}.png")


if __name__ == "__main__":
    main()
