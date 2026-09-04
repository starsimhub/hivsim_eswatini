"""Re-derive tier-C features and observations from the cached ensemble.

Wave 1's emulated feature (prev_15_49_all_mean, ages 15-50) is unaffected by the
tier_c_bands() fix, so the emulator, the NROY and the 87.2% all stand as run.
Only the 2011 top band changes, and both sides of it can be recomputed from
outputs/ensemble.parquet without touching the simulator. This script does that
and rewrites sim_results.csv / observations.csv in place.
"""
import sys, hashlib
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent.parent))
import run as R

OUT = HERE / "outputs"
BK = sorted(R.BOUNDS)
key = lambda r: hashlib.sha1(",".join(f"{float(r[k]):.10g}" for k in BK).encode()).hexdigest()[:16]

ens = pd.read_parquet(OUT / "ensemble.parquet")
D = pd.read_csv(OUT / "design.csv")
groups = dict(tuple(ens.groupby("point_key")))

rows, missing = [], 0
for _, r in D.iterrows():
    g = groups.get(key(r))
    if g is None:
        rows.append({}); missing += 1; continue
    rows.append(R.summarise_point(g))
S = pd.DataFrame(rows)
print(f"re-derived {len(S)} points ({missing} missing), {len(S.columns)} features")

obs, prov = R.build_observations()
prov = prov[prov.feature.isin(S.columns)]
S.to_csv(OUT / "sim_results.csv", index=False)
prov.to_csv(OUT / "observations.csv", index=False)

Z = pd.DataFrame({r.feature: (S[r.feature] - r.mean) / r.sigma for r in prov.itertuples()})
A = [f for f in Z if dict(zip(prov.feature, prov.tier))[f] == "A"]
n_ok = (Z.abs() < 3).sum(axis=1)
rank = pd.DataFrame({"n_ok": n_ok, "impA": Z[A].abs().max(axis=1)}).sort_values(
    ["n_ok", "impA"], ascending=[False, True])
best = rank.index[0]
print(f"\nbest joint point: design row {best}, {int(n_ok[best])}/{Z.shape[1]} within 3 sigma")
w = Z.loc[best].abs()
print("misses:", {k: round(float(v), 2) for k, v in w[w > 3].items()} or "none")
print("worst 5:", {k: round(float(v), 2) for k, v in w.sort_values(ascending=False).head(5).items()})
print(f"\ntargets-hit quantiles: " +
      str({q: int(np.percentile(n_ok, q)) for q in (50, 75, 90, 95, 99, 100)}))
print("2011 top-band z at the best point:",
      {f: round(float(Z.loc[best, f]), 2) for f in Z if "2011" in f and "_45_" in f})
D.assign(n_targets_within_3sigma=n_ok, impl_tierA=Z[A].abs().max(axis=1)
         ).to_csv(OUT / "design_scored.csv", index=False)
print("\nrewrote sim_results.csv, observations.csv, design_scored.csv")
