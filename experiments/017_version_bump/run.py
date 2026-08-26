"""Exp 017 — Version bump: stisim 1.5.8 -> 1.5.11, starsim 3.5.0 -> 3.5.2.

Three arms at fixed parameters and seeds; each ADJACENT PAIR differs by one
thing:

    A (1.5.8  / 3.5.0, all-cause mortality)
    B (1.5.11 / 3.5.2, all-cause mortality)   A->B isolates the version
    C (1.5.11 / 3.5.2, HIV-deleted mortality) B->C isolates the data

stisim and starsim are editable installs from a single checkout each, so the
version is a property of the ENVIRONMENT, not of the run. This script detects
the installed version and runs only the arms that match it, so a forgotten
checkout move cannot silently corrupt the comparison.

    # on the current stack
    python experiments/017_version_bump/run.py

    # then move the checkouts and re-run; arm A resumes from disk
    git -C ../../star_sim/stisim  checkout v1.5.11
    git -C ../../star_sim/starsim checkout v3.5.2
    python experiments/017_version_bump/run.py

Outputs
  outputs/data_hiv_deleted/                  alternate datafolder (arm C)
  outputs/sims/{arm}_{pset}_{seed}.parquet   per-run results, as each finishes
  outputs/results.parquet                    consolidated
  outputs/scorecard.csv                      headline numbers per arm/year
  outputs/death_attribution.csv              deaths among PLHIV, by module
  outputs/reproduction_check.csv             arm A vs 016 arm A, arm C vs 016 arm B
  outputs/runtime.csv                        seconds per sim, per arm
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
import starsim as ss
import stisim as sti
import yaml

EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
EXP_016 = REPO_ROOT / "experiments" / "016_double_counted_mortality"
OUT_DIR, FIG_DIR = EXP_DIR / "outputs", EXP_DIR / "figures"
SIM_DIR = OUT_DIR / "sims"
for d in (OUT_DIR, FIG_DIR, SIM_DIR):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(EXP_016))
from run_sims import make_sim  # noqa: E402
# 016 is closed and its SUMMARY immutable, so importing its construction module
# is read-only and avoids duplicating ~350 lines. Only the pure module is
# imported, never 016's run.py (which has module-level side effects).
from mortality_construction import (  # noqa: E402
    build_hiv_deleted, deleted_fraction, make_datafolder,
)

CFG = yaml.safe_load((EXP_DIR / "config.yaml").read_text())
SIM_STOP = CFG["model"]["stop"]
DEATHS_FILE = "eswatini_deaths.csv"
ALT_DATAFOLDER = OUT_DIR / "data_hiv_deleted"

# --- Arms and version gating -------------------------------------------------

# arm -> (required stisim version, HIV-deleted mortality?, extra HIV pars)
#
# Arm D isolates `rel_death_f`, the mechanism observation 6 proposes for the
# prevalence rise. It is a NEW parameter in 1.5.11 -- zero occurrences in 1.5.8
# and 1.5.10 -- giving women a 26% HIV-mortality reduction that did not exist
# before. Setting it to 1.0 removes the sex differential while leaving every
# other 1.5.11 change in place, so D vs B is the clean measurement of it.
ARMS = {
    "A_1_5_8":              ("1.5.8",  False, {}),
    "B_1_5_11":             ("1.5.11", False, {}),
    "C_1_5_11_hiv_deleted": ("1.5.11", True,  {}),
    "D_1_5_11_no_reldeathf": ("1.5.11", False, {"rel_death_f": 1.0}),
}

# 016's two sets, deliberately unchanged: arm A then reproduces 016's arm A and
# arm C reproduces 016's arm B, which is the control the whole experiment leans
# on. `high_transmission` is 014's best-fitting draw on the two parameters that
# actually drive the epidemic; at defaults the model often fails to establish.
PARAM_SETS = {
    "default": {},
    "high_transmission": {"beta_m2f": 0.0139, "rel_init_prev": 0.49},
}

# 016's arm names, for the reproduction check.
REPRO_PAIRS = {"A_1_5_8": "all_cause", "C_1_5_11_hiv_deleted": "hiv_deleted"}


def installed():
    return sti.__version__, ss.__version__


def arms_runnable():
    """Arms whose required stisim version matches what is actually installed."""
    sti_ver, _ = installed()
    return [a for a, (req, *_) in ARMS.items() if req == sti_ver]


def vmmc_class_for_version():
    """Which VMMC to use on the installed stack.

    Upstream `sti.VMMC` gained prevalence/stock-target semantics in 1.5.9 — the
    behaviour the in-repo `VMMCPrevalenceTarget` subclass was written to supply,
    and 1.5.11 additionally moved circumcision state onto HIV, which breaks the
    subclass outright. So on 1.5.11 we use upstream directly; on 1.5.8 we keep
    the subclass, because upstream there applies coverage as a per-step hazard
    and overshoots to ~100% (exp 015).

    Passing the class explicitly means NO repo edits are needed to run this
    experiment — `interventions.py` still defaults to the subclass. Actually
    deleting `vmmc.py` is an adoption step for 018, once arm B has confirmed
    the two are equivalent (metric 5).
    """
    sti_ver, _ = installed()
    if sti_ver.startswith("1.5.8"):
        return None  # interventions.py default = VMMCPrevalenceTarget
    return sti.VMMC


# --- Analyzers ---------------------------------------------------------------

class PopByAgeSex(ss.Analyzer):
    """Population and infection counts for the mortality attribution.

    Redefined here rather than imported from 016's run.py (which has
    module-level side effects). Kept byte-identical in behaviour so the
    reproduction check compares like with like.
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


class CircByAge(ss.Analyzer):
    """Circumcision prevalence by age band — the VMMC equivalence check.

    Version-agnostic by necessity: `circumcised` lives on the VMMC intervention
    in 1.5.8 and on HIV in 1.5.11. Reading the state directly (rather than a
    results column) is what makes arm A and arm B comparable at all.

    This is a required metric, not a nice-to-have: 015 established that
    prevalence-target vs hazard semantics roughly DOUBLE male HIV prevalence, so
    an upstream VMMC difference would surface in arm B looking like a mortality
    effect.
    """
    AGE_BANDS = [(15, 20), (20, 25), (25, 30), (30, 35), (35, 50), (50, 100)]

    def __init__(self, *args, name="circ", **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name

    def init_results(self):
        super().init_results()
        res = [ss.Result(f"p_circ_{lo}_{hi}", dtype=float, scale=False)
               for lo, hi in self.AGE_BANDS]
        res.append(ss.Result("p_circ_15_49", dtype=float, scale=False))
        self.define_results(*res)

    def _circumcised(self):
        hiv = self.sim.diseases.hiv
        if hasattr(hiv, "circumcised"):           # 1.5.11+
            return hiv.circumcised
        iv = self.sim.interventions.get("vmmc")   # 1.5.8
        return getattr(iv, "circumcised", None)

    def step(self):
        ti, ppl = self.ti, self.sim.people
        circ = self._circumcised()
        if circ is None:
            return
        for lo, hi in self.AGE_BANDS:
            men = ppl.alive & ppl.male & (ppl.age >= lo) & (ppl.age < hi)
            if men.count() > 0:
                self.results[f"p_circ_{lo}_{hi}"][ti] = float(np.mean(circ[men]))
        men = ppl.alive & ppl.male & (ppl.age >= 15) & (ppl.age < 50)
        if men.count() > 0:
            self.results["p_circ_15_49"][ti] = float(np.mean(circ[men]))


class PrepUptake(ss.Analyzer):
    """PrEP uptake — the confounder check (metric 8).

    `sti.Prep()` is called bare in interventions.py, and in BOTH versions
    coverage=None falls back to a built-in ramp to 80% of FSW by 2025 — a
    parameter nobody in this project chose. 1.5.11 reimplemented the mechanism
    (per-agent probability -> stock target with duration/retention), so part of
    any A->B change is PrEP rather than mortality.

    Needs an analyzer rather than a results column: 1.5.8's `Prep` defines no
    results at all, and `on_prep` moved from the intervention onto HIV in
    1.5.11.

    `on_prep` is NOT comparable across the two versions. In 1.5.8,
    `Prep.step()` computes the newly-covered agents and applies the `rel_sus`
    reduction but never writes `self.on_prep[...] = True` — the state is
    declared and left dead. Two consequences: `n_on_prep` reads 0 for the whole
    run, and because `~self.on_prep` is therefore always all-True while
    `rel_sus` is reset each timestep, coverage is re-drawn from scratch every
    month. Nobody is ever continuously protected. 1.5.11 fixes this by tracking
    enrollment on HIV with an explicit course duration.

    So `p_fsw_protected` is the metric that actually compares: the share of
    alive, uninfected female FSW carrying a susceptibility reduction right now.
    Females only, because VMMC also writes `rel_sus` — but only for males, so
    among women the reduction is PrEP's alone. `n_on_prep` is kept anyway,
    since 0-in-A vs nonzero-in-B is the direct evidence for the dead state.
    """
    def __init__(self, *args, name="prepuptake", **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result("n_on_prep", dtype=int, scale=True),
            ss.Result("p_fsw_on_prep", dtype=float, scale=False),
            ss.Result("p_fsw_protected", dtype=float, scale=False),
            ss.Result("n_fsw", dtype=int, scale=True),
        )

    def _on_prep(self):
        hiv = self.sim.diseases.hiv
        if hasattr(hiv, "on_prep"):                 # 1.5.11+
            return hiv.on_prep
        for iv in self.sim.interventions():         # 1.5.8
            if hasattr(iv, "on_prep"):
                return iv.on_prep
        return None

    def step(self):
        ti = self.ti
        ppl, hiv = self.sim.people, self.sim.diseases.hiv
        alive = ppl.alive

        on_prep = self._on_prep()
        if on_prep is not None:
            self.results["n_on_prep"][ti] = (alive & on_prep).count()

        net = self.sim.networks.get("structuredsexual")
        fsw = getattr(net, "fsw", None) if net is not None else None
        if fsw is None:
            return

        fsw_alive = alive & fsw
        self.results["n_fsw"][ti] = fsw_alive.count()
        if fsw_alive.count() == 0:
            return
        if on_prep is not None:
            self.results["p_fsw_on_prep"][ti] = float(np.mean(on_prep[fsw_alive]))

        # The cross-version measure: currently carrying a susceptibility
        # reduction. Restricted to uninfected FSW (rel_sus is meaningless once
        # infected) and implicitly female, since `fsw` is female-only.
        at_risk = fsw_alive & ~hiv.infected
        if at_risk.count() > 0:
            self.results["p_fsw_protected"][ti] = float(np.mean(hiv.rel_sus[at_risk] < 1))


# --- Runs --------------------------------------------------------------------

KEEP = ("timevec", "hiv.new_deaths", "hiv.prevalence_15_49",
        "popagesex.", "circ.", "prepuptake.")


def run_one(arm: str, pset: str, seed: int) -> pd.DataFrame:
    """One run; writes before returning so a crash mid-batch is recoverable.

    Resume matters more than usual here: once the checkouts move to 1.5.11, arm
    A cannot be reproduced without moving them back, so an existing file is
    never overwritten.
    """
    path = SIM_DIR / f"{arm}_{pset}_{seed:03d}.parquet"
    if path.exists():
        return pd.read_parquet(path)

    _, hiv_deleted, extra_pars = ARMS[arm]
    datafolder = str(ALT_DATAFOLDER) if hiv_deleted else None
    sti_ver, ss_ver = installed()

    t0 = time.perf_counter()
    sim = make_sim(seed=seed, stop=SIM_STOP, verbose=-1, datafolder=datafolder,
                   vmmc_class=vmmc_class_for_version(),
                   hiv_pars={**PARAM_SETS[pset], **extra_pars} or None,
                   analyzers=[PopByAgeSex(), CircByAge(), PrepUptake()])
    sim.run()
    elapsed = time.perf_counter() - t0

    df = sim.to_df(resample="year", use_years=True, sep=".")
    keep = [c for c in df.columns if any(c == k or c.startswith(k) for k in KEEP)]
    out = df[keep].copy()
    out["arm"], out["pset"], out["seed"] = arm, pset, seed
    out["stisim"], out["starsim"], out["runtime_s"] = sti_ver, ss_ver, elapsed
    out.to_parquet(path, index=False)
    return out


# --- Metrics -----------------------------------------------------------------

def unaids_deaths() -> pd.DataFrame:
    """UNAIDS Spectrum annual AIDS deaths — the target, from the correct file.

    NOT data/eswatini_deaths.csv, which is background mortality rates; that
    mis-citation was corrected by 016.
    """
    d = pd.read_csv(REPO_ROOT / "data" / "eswatini_hiv_calib.csv",
                    usecols=["time", "hiv.new_deaths"]).dropna()
    return d.rename(columns={"time": "year", "hiv.new_deaths": "unaids_deaths"})


def attribute_deaths(df: pd.DataFrame) -> pd.DataFrame:
    """Split deaths among PLHIV into HIV-module vs background-module.

    Stock identity, since the HIV module is the only route out of the infected
    state apart from death:
        deaths_among_plhiv[t] = n_infected[t-1] + new_infections[t] - n_infected[t]
    Same caveat as 016: migration is on, so agents entering or leaving while
    infected also move this quantity and are not separated out.
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


def reproduction_check(df: pd.DataFrame) -> pd.DataFrame:
    """Metric 1 — do our arms land on 016's published numbers?

    Arm A should reproduce 016's arm A (both 1.5.8, all-cause) and arm C should
    reproduce 016's arm B (both HIV-deleted) up to the version change. If arm A
    does NOT match, something other than the version differs and the rest of
    the experiment is uninterpretable — that is the whole point of running it.

    Returns an empty frame if 016's results are not on disk (outputs/ is
    gitignored, so a fresh clone will not have them).
    """
    src = EXP_016 / "outputs" / "results.parquet"
    if not src.exists():
        print(f"  [repro] 016 results not found at {src} — skipping check")
        return pd.DataFrame()

    old = pd.read_parquet(src)
    rows = []
    for new_arm, old_arm in REPRO_PAIRS.items():
        if new_arm not in set(df.arm):
            continue
        for pset in PARAM_SETS:
            for year in (2005, 2021):
                n = df[(df.arm == new_arm) & (df.pset == pset) & (df.timevec == year)]
                o = old[(old.arm == old_arm) & (old.pset == pset) & (old.timevec == year)]
                if not len(n) or not len(o):
                    continue
                for metric in ("hiv.prevalence_15_49", "hiv.new_deaths"):
                    nv, ov = n[metric].mean(), o[metric].mean()
                    # Two-sample z on the between-seed spread: a difference
                    # small relative to seed noise is not a difference. Must be
                    # two-sample, not (new - mean_old)/se_old — at low n_seeds
                    # the new arm's own spread dominates, and treating it as a
                    # point estimate manufactures large z from pure seed
                    # lottery. 016 measured default-parameter prevalence
                    # ranging 0.004-0.166 across 10 seeds, so this matters.
                    se_n = n[metric].std(ddof=1) / np.sqrt(len(n)) if len(n) > 1 else np.nan
                    se_o = o[metric].std(ddof=1) / np.sqrt(len(o)) if len(o) > 1 else np.nan
                    se = np.sqrt(np.nansum([se_n ** 2, se_o ** 2]))
                    rows.append({
                        "arm_017": new_arm, "arm_016": old_arm, "pset": pset,
                        "year": year, "metric": metric,
                        "v017": nv, "v016": ov, "n_017": len(n), "n_016": len(o),
                        "abs_diff": nv - ov,
                        "rel_diff": (nv - ov) / ov if ov else np.nan,
                        "se_017": se_n, "se_016": se_o,
                        "z": (nv - ov) / se if se else np.nan,
                    })
    return pd.DataFrame(rows)


# --- Figures -----------------------------------------------------------------

ARM_COLORS = {"A_1_5_8": "C3", "B_1_5_11": "C0", "C_1_5_11_hiv_deleted": "C2",
              "D_1_5_11_no_reldeathf": "C4"}
ARM_LABELS = {"A_1_5_8": "A: 1.5.8, all-cause",
              "B_1_5_11": "B: 1.5.11, all-cause",
              "C_1_5_11_hiv_deleted": "C: 1.5.11, HIV-deleted",
              "D_1_5_11_no_reldeathf": "D: 1.5.11, rel_death_f=1"}


def _band(ax, g, col, color, label):
    """Median with a min-max band across seeds."""
    s = g.groupby("timevec")[col]
    med, lo, hi = s.median(), s.min(), s.max()
    ax.plot(med.index, med.values, color=color, lw=2, label=label)
    ax.fill_between(med.index, lo.values, hi.values, color=color, alpha=0.15, lw=0)


def plot_arms(df, attr, targets, path, pset):
    """Four-panel comparison, same layout as 016's A/B figure."""
    d = df[df.pset == pset]
    a = attr[attr.pset == pset]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    for arm in ARMS:
        g = d[d.arm == arm]
        if not len(g):
            continue
        c, lab = ARM_COLORS[arm], ARM_LABELS[arm]
        _band(axes[0, 0], g, "hiv.prevalence_15_49", c, lab)
        _band(axes[0, 1], g, "hiv.new_deaths", c, lab)
        _band(axes[1, 0], g, "popagesex.n_alive_total", c, lab)
        ga = a[a.arm == arm]
        if len(ga):
            s = ga.groupby("year")["background_deaths_plhiv"]
            axes[1, 1].plot(s.median().index, s.median().values, color=c, lw=2, label=lab)

    axes[0, 1].plot(targets.year, targets.unaids_deaths, "k--", lw=2, label="UNAIDS")
    for ax, title, ylab in (
        (axes[0, 0], "A. HIV prevalence 15-49", "prevalence"),
        (axes[0, 1], "B. AIDS deaths (HIV module)", "deaths/year"),
        (axes[1, 0], "C. Total population", "people"),
        (axes[1, 1], "D. Deaths among PLHIV from the background module", "deaths/year"),
    ):
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(ylab)
        ax.set_xlabel("year")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.suptitle(f"Exp 017 — version bump, {pset} parameters "
                 f"(A->B = version, B->C = mortality data)", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_circumcision(df, path, pset):
    """Metric 5 — is upstream 1.5.11 VMMC equivalent to our subclass?"""
    d = df[df.pset == pset]
    bands = [c for c in df.columns if c.startswith("circ.p_circ_")
             and not c.endswith("15_49")]
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharey=True)
    try:
        tgt = pd.read_csv(REPO_ROOT / "data" / "vmmc_coverage.csv")
    except Exception:
        tgt = None

    for ax, col in zip(axes.flat, bands + ["circ.p_circ_15_49"]):
        for arm in ARMS:
            g = d[d.arm == arm]
            if len(g):
                _band(ax, g, col, ARM_COLORS[arm], ARM_LABELS[arm])
        ax.set_title(col.replace("circ.p_circ_", "ages "), fontsize=10)
        ax.set_xlabel("year")
        ax.grid(alpha=0.3)
    axes.flat[0].set_ylabel("proportion circumcised")
    axes.flat[0].legend(fontsize=7)
    for ax in axes.flat[len(bands) + 1:]:
        ax.axis("off")

    note = "SHIMS3 targets in data/vmmc_coverage.csv" if tgt is not None else ""
    fig.suptitle(f"Exp 017 — circumcision by age, {pset} parameters. "
                 f"Arms A and B should coincide if upstream VMMC == our subclass. {note}",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_prep(df, path):
    """Metric 8 — the PrEP confounder. If A and B diverge, part of any
    prevalence difference is the PrEP rewrite, not mortality."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    d = df[df.pset == "high_transmission"]
    for arm in ARMS:
        g = d[d.arm == arm]
        if len(g):
            _band(axes[0], g, "prepuptake.n_on_prep", ARM_COLORS[arm], ARM_LABELS[arm])
            _band(axes[1], g, "prepuptake.p_fsw_on_prep", ARM_COLORS[arm], ARM_LABELS[arm])
            _band(axes[2], g, "prepuptake.p_fsw_protected", ARM_COLORS[arm], ARM_LABELS[arm])
    axes[0].set_title("Number on PrEP (state)", fontsize=11)
    axes[1].set_title("Share of FSW with on_prep set\n"
                      "(dead state in 1.5.8 — reads 0)", fontsize=10)
    axes[2].set_title("Share of uninfected FSW protected\n"
                      "(rel_sus < 1 — the comparable measure)", fontsize=10)
    for ax in (axes[1], axes[2]):
        ax.axhline(0.5, color="k", ls=":", lw=1)
        ax.axhline(0.8, color="k", ls=":", lw=1)
    for ax in axes:
        ax.set_xlabel("year")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("Exp 017 — inherited PrEP default (nobody chose this: "
                 "50% of FSW by 2015, 80% by 2025), high-transmission parameters",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# --- Main --------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_seeds", type=int, default=CFG["model"]["n_seeds"])
    p.add_argument("--n_workers", type=int, default=None)
    p.add_argument("--plot_only", action="store_true")
    args = p.parse_args()

    sti_ver, ss_ver = installed()
    runnable = arms_runnable()
    print(f"Installed: stisim {sti_ver}, starsim {ss_ver}")
    print(f"Arms runnable on this stack: {runnable or '(none)'}")
    if not runnable and not args.plot_only:
        expected = sorted({v for v, _ in ARMS.values()})
        raise SystemExit(
            f"\nNo arm matches stisim {sti_ver}. This experiment expects one of "
            f"{expected}.\nMove the editable checkouts to a tagged version, e.g.\n"
            f"  git -C ../../star_sim/stisim  checkout v1.5.11\n"
            f"  git -C ../../star_sim/starsim checkout v3.5.2\n"
        )

    print("Building HIV-deleted mortality (arm C)...")
    all_cause = pd.read_csv(REPO_ROOT / "data" / DEATHS_FILE)
    hiv_deleted = build_hiv_deleted(all_cause)
    construction = deleted_fraction(all_cause, hiv_deleted)
    construction.to_csv(OUT_DIR / "mortality_construction.csv", index=False)
    make_datafolder(REPO_ROOT / "data", ALT_DATAFOLDER, hiv_deleted, DEATHS_FILE)
    n_changed = int((construction.deleted_rate > 1e-12).sum())
    print(f"  {n_changed} of {len(construction)} rows lowered")

    if not args.plot_only:
        work = [(arm, pset, seed) for arm in runnable for pset in PARAM_SETS
                for seed in range(args.n_seeds)]
        todo = [w for w in work
                if not (SIM_DIR / f"{w[0]}_{w[1]}_{w[2]:03d}.parquet").exists()]
        print(f"{len(work) - len(todo)} of {len(work)} already on disk; "
              f"running {len(todo)} (n_workers={args.n_workers or 'all'})...")
        if todo:
            t0 = sc.tic()
            sc.parallelize(run_one, iterarg=todo, ncpus=args.n_workers)
            sc.toc(t0, label="017 run")

    # Consolidate everything on disk, across however many stacks have been run.
    files = sorted(SIM_DIR.glob("*.parquet"))
    if not files:
        raise SystemExit("No results on disk yet.")
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    df.to_parquet(OUT_DIR / "results.parquet", index=False)
    present = sorted(set(df.arm))
    print(f"\nArms with results on disk: {present}")

    rt = (df.groupby(["arm", "stisim", "starsim"])["runtime_s"]
            .agg(["mean", "min", "max", "count"]).reset_index())
    rt.to_csv(OUT_DIR / "runtime.csv", index=False)

    targets = unaids_deaths()
    attr = attribute_deaths(df)
    attr.to_csv(OUT_DIR / "death_attribution.csv", index=False)

    repro = reproduction_check(df)
    if len(repro):
        repro.to_csv(OUT_DIR / "reproduction_check.csv", index=False)

    rows = []
    for pset in PARAM_SETS:
        for year in (2005, 2021):
            r = {"pset": pset, "year": year}
            for arm in present:
                g = df[(df.arm == arm) & (df.pset == pset) & (df.timevec == year)]
                if not len(g):
                    continue
                r[f"prev_15_49_{arm}"] = g["hiv.prevalence_15_49"].mean()
                r[f"hiv_deaths_{arm}"] = g["hiv.new_deaths"].mean()
                r[f"pop_{arm}"] = g["popagesex.n_alive_total"].mean()
                r[f"p_circ_15_49_{arm}"] = g["circ.p_circ_15_49"].mean()
                r[f"p_fsw_prep_{arm}"] = g["prepuptake.p_fsw_on_prep"].mean()
                r[f"p_fsw_protected_{arm}"] = g["prepuptake.p_fsw_protected"].mean()
            t = targets[targets.year == year]
            r["unaids_deaths"] = float(t.unaids_deaths.iloc[0]) if len(t) else np.nan
            rows.append(r)
    score = pd.DataFrame(rows)
    score.to_csv(OUT_DIR / "scorecard.csv", index=False)

    print("\n=== Scorecard ===")
    for _, r in score.iterrows():
        print(f"\n{r.pset}, {int(r.year)}  (UNAIDS deaths {r.unaids_deaths:,.0f})")
        for arm in present:
            if f"prev_15_49_{arm}" not in r or pd.isna(r.get(f"prev_15_49_{arm}")):
                continue
            print(f"  {ARM_LABELS[arm]:<28} prev {r[f'prev_15_49_{arm}']:.4f}"
                  f"   deaths {r[f'hiv_deaths_{arm}']:>8,.0f}"
                  f"   pop {r[f'pop_{arm}']:>10,.0f}"
                  f"   circ {r[f'p_circ_15_49_{arm}']:.3f}"
                  f"   fsw-prot {r[f'p_fsw_protected_{arm}']:.3f}"
                  f"   on_prep {r[f'p_fsw_prep_{arm}']:.3f}")

    if len(repro):
        print("\n=== Reproduction check vs 016 ===")
        print("(z = two-sample, on the pooled between-seed standard error)")
        for _, r in repro.iterrows():
            flag = "  <-- CHECK" if abs(r.z) > 3 else ""
            print(f"  {r.arm_017:<24} {r.pset:<18} {int(r.year)} "
                  f"{r.metric:<24} {r.v017:>10,.4f} vs {r.v016:>10,.4f} "
                  f"({r.rel_diff:+.1%}, z={r.z:+.1f}){flag}")
    elif "A_1_5_8" in present:
        print("\n[repro] 016 results unavailable — reproduction check skipped")

    if len(present) < len(ARMS):
        missing = [a for a in ARMS if a not in present]
        print(f"\nStill to run: {missing}")
        print("Move the editable checkouts and re-run; arms already on disk resume.")

    for pset in PARAM_SETS:
        plot_arms(df, attr, targets, FIG_DIR / f"arms_{pset}.png", pset)
        plot_circumcision(df, FIG_DIR / f"circumcision_{pset}.png", pset)
    plot_prep(df, FIG_DIR / "prep_uptake.png")
    print(f"\nWrote outputs to {OUT_DIR} and figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
