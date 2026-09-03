"""Wave-1 figures. Separate from run.py so figures can be rebuilt without compute.

Two figures, answering two different questions:

  prevalence_fit_vs_phia.png  -- does the PRIOR ENSEMBLE cover PHIA? This is the
      coverage check (workflow step 3) that 009 and 014 both failed. Uses
      kind='ensemble' so the band is the 5-95% over the 1000 draws, not over
      seeds.

  best_joint_point.png        -- can a SINGLE parameter point fit everything at
      once? Covering the data draw-by-draw is necessary but not sufficient; if
      every target is reachable but only by mutually exclusive points, the model
      still cannot fit. This is the figure that distinguishes those two cases.

Both call standard_figures.plot_prevalence_fit so the format cannot drift from
the rest of the project (see CLAUDE.md).
"""
import sys, hashlib, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))
import standard_figures as sf

OUT, FIG = HERE / "outputs", HERE / "figures"
FIG.mkdir(exist_ok=True)
BOUND_KEYS = ["age_gap_shift", "log_age_gap_sd_mult", "log_beta_m2f",
              "log_rel_beta_f2m", "log_s_f_young", "prop_f0", "prop_m0"]


def point_key(row):
    """Must match run.py.point_key exactly, or the join silently returns nothing."""
    vec = ",".join(f"{float(row[k]):.10g}" for k in sorted(BOUND_KEYS))
    return hashlib.sha1(vec.encode()).hexdigest()[:16]


def main():
    ens = pd.read_parquet(OUT / "ensemble.parquet")
    D = pd.read_csv(OUT / "design_scored.csv")
    tg = sf.load_targets()

    # --- 1. ensemble coverage ---
    sf.plot_prevalence_fit(ens, "024 wave-1 prior ensemble (1000 draws)",
                           FIG / "prevalence_fit_vs_phia.png",
                           kind="ensemble", tg=tg,
                           stamp="HM wave 1, N=10,000, 1 replicate/point")

    sf.plot_incidence_fit(ens, "024 wave-1 prior ensemble (1000 draws)",
                          FIG / "incidence_fit_vs_shims.png",
                          kind="ensemble",
                          stamp="HM wave 1, N=10,000, 1 replicate/point")

    # --- 2. the best joint point ---
    best = D.sort_values(["n_targets_within_3sigma", "impl_tierA"],
                         ascending=[False, True]).iloc[0]
    key = point_key(best)
    sub = ens[ens.point_key == key]
    if not len(sub):
        raise SystemExit(f"point_key {key} not in ensemble -- hash mismatch")
    lab = (f"024 best joint point ({int(best.n_targets_within_3sigma)}/48 within 3σ)  "
           f"β={np.exp(best.log_beta_m2f):.4f}, rel_f2m={np.exp(best.log_rel_beta_f2m):.3f}, "
           f"s_f_young={np.exp(best.log_s_f_young):.2f}")
    sf.plot_prevalence_fit(sub, lab, FIG / "best_joint_point.png",
                           kind="arm", tg=tg, stamp="single draw, 1 replicate")
    sf.plot_incidence_fit(sub, lab, FIG / "best_joint_point_incidence.png",
                          kind="arm", stamp="single draw, 1 replicate")

    # --- 3. what the wave constrained, and what it did not ---
    S = pd.read_csv(OUT / "sim_results.csv"); o = pd.read_csv(OUT / "observations.csv")
    o = o[o.feature.isin(S.columns)]
    Z = pd.DataFrame({r.feature: (S[r.feature] - r.mean) / r.sigma for r in o.itertuples()})
    best_z = Z.abs().min()                       # closest any single draw gets
    joint_z = Z.loc[best.name].abs()             # the joint point's own residual

    fig, axes = plt.subplots(2, 1, figsize=(14, 9))
    order = joint_z.sort_values().index
    x = np.arange(len(order))
    ax = axes[0]
    ax.axhspan(0, 3, color="#e8f4e8", zorder=0)
    ax.bar(x - 0.2, best_z[order], 0.4, label="closest any single draw gets", color="#4a7fb5")
    ax.bar(x + 0.2, joint_z[order], 0.4, label="the best JOINT draw", color="#c0392b")
    ax.axhline(3, ls="--", c="k", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_ylabel("|z| = |sim − target| / σ"); ax.legend(fontsize=8)
    n_over = int((joint_z > 3).sum())
    ax.set_title(f"Every target is individually reachable, and one draw reaches "
                 f"{int(best.n_targets_within_3sigma)} of {len(joint_z)} jointly.\n"
                 + ("No target exceeds 3σ at that draw." if not n_over else
                    f"The {n_over} red bar(s) above 3σ are the entire residual defect."),
                 fontsize=10)

    ax = axes[1]
    disc = np.array([0.0, 0.005, 0.010, 0.015, 0.020])
    prev = [f for f in Z if f.startswith("prev_15_49") and "_all_" not in f]
    ci = {r.feature: np.sqrt(max(r.sigma**2 - 0.02**2, 1e-12))
          for r in o.itertuples() if r.feature in prev}
    frac, bestv = [], []
    for dsc in disc:
        zz = pd.DataFrame({f: (S[f] - float(o.loc[o.feature == f, "mean"].iloc[0]))
                           / np.sqrt(ci[f]**2 + dsc**2) for f in prev})
        m = zz.abs().max(axis=1)
        frac.append((m < 4).mean() * 100); bestv.append(m.min())
    ax.plot(disc, bestv, "o-", color="#c0392b", label="best draw's max |z| over the six PHIA 15-49 targets")
    ax.axhline(4, ls="--", c="k", lw=1, label="implausibility threshold 4")
    ax.set_xlabel("assumed model discrepancy σ_disc on prevalence")
    ax.set_ylabel("max |z| of the best draw"); ax.legend(fontsize=8)
    ax.set_title("Even at σ_disc = 0 the best draw sits at 3.98σ — the −4.3 pp deficit of 016–022\n"
                 "was a parameter-value deficit, not the structural bias σ_disc was set to absorb.",
                 fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "target_residuals.png", dpi=130); plt.close(fig)

    json.dump({"best_design_row": int(best.name),
               "n_targets_within_3sigma": int(best.n_targets_within_3sigma),
               "impl_tierA": float(best.impl_tierA),
               "params_natural": {"beta_m2f": float(np.exp(best.log_beta_m2f)),
                                  "rel_beta_f2m": float(np.exp(best.log_rel_beta_f2m)),
                                  "s_f_young": float(np.exp(best.log_s_f_young)),
                                  "age_gap_shift": float(best.age_gap_shift),
                                  "age_gap_sd_mult": float(np.exp(best.log_age_gap_sd_mult)),
                                  "prop_f0": float(best.prop_f0), "prop_m0": float(best.prop_m0)},
               "misses_over_3sigma": {k: round(float(v), 2)
                                      for k, v in joint_z[joint_z > 3].items()},
               "sigma_disc_scan": dict(zip([f"{d:.3f}" for d in disc],
                                           [round(float(b), 2) for b in bestv]))},
              open(OUT / "best_point.json", "w"), indent=2)
    print("figures ->", FIG)
    print("best point:", json.load(open(OUT / "best_point.json"))["params_natural"])


if __name__ == "__main__":
    main()
