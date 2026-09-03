"""Exp 023 — Prior sensitivity sweep over nine candidate wave-1 parameters.

400 draws, 1 replicate each, N = 10,000. Deliberately NOT the coverage check:
024 runs that at N = 20,000 with 10 replicates per 020's sizing. This asks the
smaller, cheaper question -- does the prior move the targets, and can the data
tell the parameters apart?

  beta_m2f         0.008-0.025  log   per-act M->F risk, reference band
  rel_beta_f2m     0.15-0.60    log   male:female direction ratio
  s_f_young        0.8-3.0      log   women 15-24 susceptibility multiplier
  rel_init_prev    0.1-0.5      lin   1985 seed prevalence
  age_gap_shift    -2 to +3 yr  lin   additive, all 9 age_diff_pars means
  age_gap_sd_mult  0.6-1.8      log   multiplicative, all 9 SDs (assortativity)
  prop_f0          0.45-0.85    lin   female low-risk share
  prop_m0          0.40-0.80    lin   male low-risk share
  conc_mult        0.5-2.0      log   on f1/f2/m1/m2 concurrency

The 25-34 susceptibility band is anchored at 1.0 for both sexes -- a reference
point, not a biological claim. Without it, halving beta_m2f and doubling every
multiplier gives an identical model.

Outputs
  outputs/sims/draw_{i}.parquet   per-draw, written as each finishes
  outputs/draws.csv               the prior sample
  outputs/ensemble.parquet        concatenated trajectories
  outputs/summary.csv             per-draw target statistics
  outputs/sensitivity.csv         Spearman rho, every parameter x every target
  outputs/confounding.csv         effect-signature correlation between parameters
  outputs/orthogonality.csv       prior-draw correlations (uninformative -- see
                                  effect_signature_confounding for why)
  outputs/diagnostics.csv         the specific questions from the README
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
from scipy.stats import spearmanr

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
OUT_DIR, FIG_DIR = EXP_DIR / "outputs", EXP_DIR / "figures"
SIM_DIR = OUT_DIR / "sims"
for d in (OUT_DIR, FIG_DIR, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
from run_sims import make_sim  # noqa: E402
from analyzers import PopByAgeSex, Cascade  # noqa: E402
from standard_figures import load_targets, fit_by_stratum  # noqa: E402

N_AGENTS = 10_000
STOP = 2026
N_DRAWS = 400
DRAW_SEED = 20260902
SIM_SEED = 1          # one replicate per draw; the draw index varies, not the seed

# (low, high, 'log' or 'lin')
PRIOR = {
    # --- transmission and susceptibility ---
    "beta_m2f":       (0.008, 0.025, "log"),   # per-act M->F, reference band
    "rel_beta_f2m":   (0.15,  0.60,  "log"),   # male:female direction ratio
    "s_f_young":      (0.8,   3.0,   "log"),   # women 15-24 susceptibility
    # --- seeding ---
    "rel_init_prev":  (0.1,   0.5,   "lin"),
    # --- mixing: the mechanism seven experiments implicated and none tested ---
    "age_gap_shift":  (-2.0,  3.0,   "lin"),   # additive, on all 9 age_diff_pars means
    "age_gap_sd_mult":(0.6,   1.8,   "log"),   # multiplicative, on all 9 SDs
    # --- risk structure ---
    "prop_f0":        (0.45,  0.85,  "lin"),   # female low-risk share
    "prop_m0":        (0.40,  0.80,  "lin"),   # male low-risk share
    "conc_mult":      (0.5,   2.0,   "log"),   # on f1/f2/m1/m2 concurrency
}

# stisim defaults for age_diff_pars, (mean, sd) by female age group x risk group.
# age_gap_shift and age_gap_sd_mult are single scalars over all nine, which
# preserves the relative structure -- teens partner with slightly younger men
# than adults, higher-risk groups with smaller and tighter gaps -- while giving
# the mixing hypothesis one degree of freedom each instead of eighteen.
AGE_DIFF_BASE = {
    "teens": [(7, 3), (6, 3), (5, 1)],
    "young": [(8, 3), (7, 3), (5, 2)],
    "adult": [(8, 3), (7, 3), (5, 2)],
}
# make_sim's concurrency values, scaled by conc_mult. f0/m0 stay at the stisim
# default of ~0 (group 0 is strictly monogamous by construction).
CONC_BASE = {"f1_conc": 0.15, "f2_conc": 0.25, "m1_conc": 0.15, "m2_conc": 0.5}

BANDS = [(15, 25), (25, 35), (35, 50)]
INCIDENCE_TARGETS = REPO_ROOT / "calibration_data" / "incidence_by_age_sex.csv"
ESTABLISH_YEAR, ESTABLISH_THRESHOLD = 2005, 0.05

KEEP = ("timevec", "hiv.new_deaths", "hiv.prevalence_15_49", "popagesex.",
        "cascade.p_vls_")


def sample_prior(n=N_DRAWS, seed=DRAW_SEED):
    rng = np.random.default_rng(seed)
    u = rng.uniform(size=(n, len(PRIOR)))
    cols = {"draw": np.arange(n)}
    for k, (name, (lo, hi, scale)) in enumerate(PRIOR.items()):
        if scale == "log":
            cols[name] = np.exp(np.log(lo) + u[:, k] * (np.log(hi) - np.log(lo)))
        else:
            cols[name] = lo + u[:, k] * (hi - lo)
    return pd.DataFrame(cols)


def _run(draw, pars):
    path = SIM_DIR / f"draw_{draw:04d}.parquet"
    if path.exists():
        return

    # s_f_young enters through rel_sus_age. The 25-34 and 35+ female bands and
    # all male bands stay at 1.0 -- the anchor that makes beta_m2f identifiable.
    rel_sus_age = [(15, 25, 'f', float(pars["s_f_young"])),
                   (25, 50, 'f', 1.0),
                   (15, 50, 'm', 1.0)]
    hiv_pars = dict(beta_m2f=float(pars["beta_m2f"]),
                    rel_beta_f2m=float(pars["rel_beta_f2m"]),
                    rel_init_prev=float(pars["rel_init_prev"]),
                    rel_sus_age=rel_sus_age)

    shift, sd_mult = float(pars["age_gap_shift"]), float(pars["age_gap_sd_mult"])
    age_diff_pars = {
        grp: [(max(m + shift, 1.0), max(s * sd_mult, 0.2)) for m, s in vals]
        for grp, vals in AGE_DIFF_BASE.items()}
    cm = float(pars["conc_mult"])
    network_pars = dict(
        prop_f0=float(pars["prop_f0"]), prop_m0=float(pars["prop_m0"]),
        age_diff_pars=age_diff_pars,
        **{k: v * cm for k, v in CONC_BASE.items()})

    t0 = time.perf_counter()
    sim = make_sim(seed=SIM_SEED, stop=STOP, verbose=-1, hiv_pars=hiv_pars,
                   network_pars=network_pars,
                   analyzers=[PopByAgeSex(), Cascade()])
    sim.pars.n_agents = N_AGENTS
    sim.run()
    elapsed = time.perf_counter() - t0

    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns if any(c == k or c.startswith(k) for k in KEEP)]
    out = df[keep].copy()
    out["draw"], out["runtime_s"] = draw, elapsed
    for k in PRIOR:
        out[k] = float(pars[k])
    out.to_parquet(path, index=False)


def load():
    files = sorted(SIM_DIR.glob("draw_*.parquet"))
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


# --- Summary statistics per draw ---------------------------------------------

def model_incidence(g, sex, lo, hi, year):
    """Annual incidence per 100 PY for one band, from PopByAgeSex."""
    sub = g[g.timevec == year]
    if not len(sub):
        return np.nan
    # PopByAgeSex bins start at multiples of 5, so a target band whose edges
    # are not multiples of 5 has to be snapped. SHIMS1 2011 reports 18-49; the
    # model is read as 15-49. That includes 15-17-year-olds, who have very low
    # incidence, so it biases the MODEL slightly DOWN relative to an 18-49
    # measure -- i.e. the comparison is mildly conservative, not flattering.
    lo = (lo // 5) * 5
    hi = -(-hi // 5) * 5
    inf = alive = infected = 0.0
    for b in range(lo, hi, 5):
        n_col = f"popagesex.new_infections_{sex}_{b}_{b + 5}"
        a_col = f"popagesex.n_alive_{sex}_{b}_{b + 5}"
        i_col = f"popagesex.n_infected_{sex}_{b}_{b + 5}"
        if n_col not in g.columns:
            return np.nan
        inf += float(sub[n_col].iloc[0])
        alive += float(sub[a_col].iloc[0])
        infected += float(sub[i_col].iloc[0])
    susc = alive - infected
    return 100.0 * inf / susc if susc > 0 else np.nan


def summarise(df, tg):
    """One row per draw: the target-facing statistics sensitivity is measured on."""
    inc_tg = pd.read_csv(INCIDENCE_TARGETS)
    rows = []
    for draw, g in df.groupby("draw"):
        r = {"draw": int(draw)}
        for k in PRIOR:
            r[k] = float(g[k].iloc[0])

        est = g[g.timevec == ESTABLISH_YEAR]["hiv.prevalence_15_49"]
        r["established"] = bool(len(est) and est.iloc[0] > ESTABLISH_THRESHOLD)

        # Prevalence: level, and the shape statistics 014's single-stratum
        # ranking would have scored at ~0
        st = fit_by_stratum(g, tg=tg)
        if len(st):
            st["band"] = pd.cut(st.age_low, [14, 24, 34, 44, 64],
                                labels=["15-24", "25-34", "35-44", "45-64"])
            piv = st.pivot_table(index="band", columns="sex", values="model",
                                 observed=True)
            for band in piv.index:
                for sex in ("f", "m"):
                    if sex in piv.columns:
                        r[f"prev_{sex}_{band}"] = piv.loc[band, sex]
            r["prev_level"] = st.model.mean()
            fy = piv.loc["15-24", "f"] if "f" in piv else np.nan
            fo = piv.loc["35-44", "f"] if "f" in piv else np.nan
            mm = piv["m"].mean() if "m" in piv else np.nan
            ff = piv["f"].mean() if "f" in piv else np.nan
            r["prev_young_old_f"] = fy / fo if fo else np.nan
            r["prev_fm_ratio"] = ff / mm if mm else np.nan
            r["prev_mae"] = (st.model - st.phia).abs().mean()

        # Incidence, on the target bands
        for _, t in inc_tg.iterrows():
            key = f"inc_{t.sex}_{int(t.age_low)}_{int(t.age_high)}_{int(t.year)}"
            r[key] = model_incidence(g, t.sex, int(t.age_low), int(t.age_high),
                                     int(t.year))
        # F:M ratio from whichever bands that round actually publishes -- 2011
        # is a single 18-50 band, 2016 is three bands. Averaging over the year's
        # own keys avoids assuming a band structure that may not exist.
        for year in (2011, 2016):
            vals = {}
            for sex in ("m", "f"):
                keys = [k for k in r
                        if k.startswith(f"inc_{sex}_") and k.endswith(f"_{year}")]
                got = [r[k] for k in keys if r.get(k) is not None
                       and np.isfinite(r[k])]
                vals[sex] = float(np.mean(got)) if got else np.nan
            r[f"inc_fm_ratio_{year}"] = (vals["f"] / vals["m"]
                                         if vals["m"] and np.isfinite(vals["m"])
                                         else np.nan)

        d = g.groupby("timevec")["hiv.new_deaths"].mean()
        r["peak_deaths"] = d.max() if len(d) else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


# --- Sensitivity --------------------------------------------------------------

def sensitivity(s):
    ok = s[s.established]
    targets = [c for c in s.columns
               if c.startswith(("prev_", "inc_")) or c == "peak_deaths"]
    rows = []
    for par in PRIOR:
        for tgt in targets:
            v = ok[[par, tgt]].dropna()
            if len(v) < 10:
                continue
            rho = spearmanr(v[par], v[tgt]).statistic
            rows.append(dict(parameter=par, target=tgt, rho=rho,
                             abs_rho=abs(rho), n=len(v)))
    if not rows:   # too few established draws to correlate anything
        return pd.DataFrame(columns=["parameter", "target", "rho", "abs_rho", "n"])
    return pd.DataFrame(rows).sort_values("abs_rho", ascending=False)


def effect_signature_confounding(sens):
    """Which parameters DO THE SAME THING to the model.

    The pairwise correlation of prior draws is uninformative here: the design is
    an independent uniform sample, so those correlations are ~0 by construction
    and confirm only that the sampler works. The question that matters is
    whether two parameters have the same effect signature across targets -- if
    they do, the data cannot tell them apart however they were sampled.
    """
    piv = sens.pivot_table(index="target", columns="parameter",
                           values="rho").dropna()
    C = piv.corr(method="spearman")
    rows = []
    names = list(C.columns)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            r = C.loc[a, b]
            rows.append(dict(par_a=a, par_b=b, signature_r=r, abs_r=abs(r),
                             n_targets=len(piv),
                             flag="CONFOUNDED" if abs(r) > 0.8 else ""))
    return pd.DataFrame(rows).sort_values("abs_r", ascending=False)


def orthogonality(s):
    ok = s[s.established]
    rows = []
    names = list(PRIOR)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            rho = spearmanr(ok[a], ok[b]).statistic
            rows.append(dict(par_a=a, par_b=b, rho=rho, abs_rho=abs(rho),
                             flag="NON-IDENTIFIABLE" if abs(rho) > 0.8 else ""))
    return pd.DataFrame(rows).sort_values("abs_rho", ascending=False)


def plot_sensitivity(sens, s, path):
    fig, axes = plt.subplots(1, 3, figsize=(19, 5))

    ax = axes[0]
    fam = {"prev_f": "prevalence, women", "prev_m": "prevalence, men",
           "inc_f": "incidence, women", "inc_m": "incidence, men",
           "peak_deaths": "peak deaths"}
    best = []
    for par in PRIOR:
        sub = sens[sens.parameter == par]
        for pref, lbl in fam.items():
            f = sub[sub.target.str.startswith(pref)]
            best.append(dict(parameter=par, family=lbl,
                             max_abs_rho=f.abs_rho.max() if len(f) else np.nan))
    b = pd.DataFrame(best).pivot(index="parameter", columns="family",
                                 values="max_abs_rho")
    im = ax.imshow(b.values, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set(xticks=range(len(b.columns)), yticks=range(len(b.index)),
           title="max |Spearman rho| by parameter and target family")
    ax.set_xticklabels(b.columns, rotation=20, ha="right", fontsize=8)
    ax.set_yticklabels(b.index, fontsize=9)
    for i in range(len(b.index)):
        for j in range(len(b.columns)):
            v = b.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="w" if v < 0.6 else "k", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.04)

    ok = s[s.established]
    ax = axes[1]
    if "prev_f_15-24" in ok:
        ax.scatter(ok.s_f_young, ok["prev_f_15-24"], s=14, alpha=0.7,
                   label="women 15-24")
        if "prev_f_35-44" in ok:
            ax.scatter(ok.s_f_young, ok["prev_f_35-44"], s=14, alpha=0.7,
                       label="women 35-44")
    ax.set(xscale="log", xlabel="s_f_young", ylabel="model prevalence",
           title="Does s_f_young move young women\nspecifically, or everything?")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[2]
    if "inc_fm_ratio_2016" in ok:
        ax.scatter(ok.rel_beta_f2m, ok["inc_fm_ratio_2016"], s=14, alpha=0.7)
    for y, lbl in ((1.903, "observed 2011"), (2.035, "observed 2016")):
        ax.axhline(y, ls="--", lw=1, color="k")
        ax.annotate(lbl, (ax.get_xlim()[0], y), fontsize=7, va="bottom")
    ax.set(xscale="log", xlabel="rel_beta_f2m",
           ylabel="model F:M incidence ratio",
           title="Can rel_beta_f2m reach the\nobserved sex ratio?")
    ax.grid(alpha=0.3)

    fig.suptitle(f"Exp 023 — prior sensitivity over {len(PRIOR)} candidate parameters", y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_draws", type=int, default=N_DRAWS)
    p.add_argument("--n_workers", type=int, default=None)
    p.add_argument("--plot_only", action="store_true")
    args = p.parse_args()

    draws = sample_prior(args.n_draws)
    draws.to_csv(OUT_DIR / "draws.csv", index=False)
    print(f"prior sample: {len(draws)} draws over {list(PRIOR)}")
    print(draws.describe().loc[["min", "50%", "max"]].round(4).to_string())

    if not args.plot_only:
        todo = [(int(r.draw), r._asdict() if hasattr(r, "_asdict") else dict(r))
                for r in draws.itertuples(index=False)
                if not (SIM_DIR / f"draw_{int(r.draw):04d}.parquet").exists()]
        print(f"\n{len(draws) - len(todo)}/{len(draws)} on disk; running {len(todo)}")
        if todo:
            t0 = sc.tic()
            sc.parallelize(_run, iterarg=todo, ncpus=args.n_workers)
            sc.toc(t0, label="023")

    df = load()
    if not len(df):
        print("no results")
        return
    df.to_parquet(OUT_DIR / "ensemble.parquet", index=False)
    tg = load_targets()
    s = summarise(df, tg)
    s.to_csv(OUT_DIR / "summary.csv", index=False)

    n_est = int(s.established.sum())
    print(f"\n=== Established epidemics: {n_est}/{len(s)} "
          f"({100 * n_est / len(s):.0f}%) ===")
    if n_est < len(s):
        dead = s[~s.established]
        print(f"  failed draws had beta_m2f {dead.beta_m2f.min():.4f}-"
              f"{dead.beta_m2f.max():.4f}, rel_init_prev "
              f"{dead.rel_init_prev.min():.3f}-{dead.rel_init_prev.max():.3f}")

    sens = sensitivity(s)
    sens.to_csv(OUT_DIR / "sensitivity.csv", index=False)
    print("\n=== Strongest target for each parameter ===")
    for par in PRIOR:
        sub = sens[sens.parameter == par]
        if len(sub):
            top = sub.iloc[0]
            print(f"  {par:14s} rho = {top.rho:+.2f}  vs  {top.target}")

    orth = orthogonality(s)
    orth.to_csv(OUT_DIR / "orthogonality.csv", index=False)
    conf = effect_signature_confounding(sens)
    conf.to_csv(OUT_DIR / "confounding.csv", index=False)
    print("\n=== Which parameters do the same thing? (|r| > 0.8 = confounded) ===")
    print("(prior-draw correlations are ~0 by construction and say nothing;")
    print(" this compares effect signatures across targets instead)")
    print(conf.head(6).round(3).to_string(index=False))

    print("\n=== The README's specific questions ===")
    ok = s[s.established]
    diag = []
    if "prev_f_15-24" in ok and "prev_f_35-44" in ok:
        r_young = spearmanr(ok.s_f_young, ok["prev_f_15-24"]).statistic
        r_old = spearmanr(ok.s_f_young, ok["prev_f_35-44"]).statistic
        diag.append(dict(question="s_f_young moves young women specifically?",
                         answer=f"rho young {r_young:+.2f} vs older {r_old:+.2f}",
                         verdict="specific" if abs(r_young) - abs(r_old) > 0.2
                                 else "acts like a level parameter"))
    if "prev_m_25-34" in ok and "prev_m_35-44" in ok:
        for par in ("age_gap_shift", "age_gap_sd_mult", "prop_m0", "conc_mult"):
            r25 = spearmanr(ok[par], ok["prev_m_25-34"]).statistic
            r35 = spearmanr(ok[par], ok["prev_m_35-44"]).statistic
            diag.append(dict(
                question=f"{par} reaches the male age profile?",
                answer=f"rho vs men 25-34 {r25:+.2f}, vs men 35-44 {r35:+.2f}",
                verdict=("differential" if abs(r25 - r35) > 0.2
                         else "acts uniformly on men")))
    if "inc_fm_ratio_2016" in ok:
        rng = ok["inc_fm_ratio_2016"].dropna()
        diag.append(dict(question="F:M incidence ratio reachable (obs 1.90-2.04)?",
                         answer=f"model spans {rng.min():.2f}-{rng.max():.2f}",
                         verdict="reachable" if rng.min() <= 2.0 <= rng.max()
                                 else "NOT reachable"))
    d = pd.DataFrame(diag)
    if len(d):
        d.to_csv(OUT_DIR / "diagnostics.csv", index=False)
        print(d.to_string(index=False))

    plot_sensitivity(sens, s, FIG_DIR / "sensitivity.png")
    print(f"\nfigures -> {FIG_DIR}")


if __name__ == "__main__":
    main()
