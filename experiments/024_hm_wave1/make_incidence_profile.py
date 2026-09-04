"""The continuous (5-year band) modelled incidence age profile, 2011 and 2016.

Separate from make_figures.py because it answers a different question: not "does
the model match the banded targets" but "what age profile does the model actually
produce". The banded figure cannot show that -- SHIMS publishes three coarse
bands, and 2011 has no age detail at all.
"""
import sys, hashlib
from pathlib import Path
import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))
import standard_figures as sf

OUT, FIG = HERE / "outputs", HERE / "figures"
BK = ["age_gap_shift", "log_age_gap_sd_mult", "log_beta_m2f", "log_rel_beta_f2m",
      "log_s_f_young", "prop_f0", "prop_m0"]


def point_key(row):
    vec = ",".join(f"{float(row[k]):.10g}" for k in sorted(BK))
    return hashlib.sha1(vec.encode()).hexdigest()[:16]


ens = pd.read_parquet(OUT / "ensemble.parquet")
D = pd.read_csv(OUT / "design_scored.csv")
best = D.sort_values(["n_targets_within_3sigma", "impl_tierA"],
                     ascending=[False, True]).iloc[0]

sf.plot_incidence_age_profile(
    ens, "024 wave-1 prior ensemble (1000 draws)",
    FIG / "incidence_age_profile.png", kind="ensemble",
    stamp="median and 5-95% of draws, N=10,000")
sf.plot_incidence_age_profile(
    ens[ens.point_key == point_key(best)],
    "024 best joint point (48/48 within 3σ)",
    FIG / "incidence_age_profile_best.png", kind="arm",
    stamp="single draw, 1 replicate")
print("wrote incidence_age_profile.png and incidence_age_profile_best.png")
