"""Shared helpers for turning `analyzers.CascadeByAge` counts into cascade rates.

Promoted from experiments/029_cascade_audit/analyse.py when exp 030 needed the
same aggregation. 029's copy is left in place: it is closed, and its code is the
record of what actually ran there.

The analyzer stores 5-year bands so that any grouping is a downstream sum --
SHIMS3 reports 15-24 / 25-34 / 35-49 / 50+, while data/art_coverage.csv drives
the model on [15,25) / [25,35) / [35,45) / [45,100), and those disagree above 35.
"""

import numpy as np
import pandas as pd

BANDS = [(a, a + 5) for a in range(15, 80, 5)] + [(80, 200)]

SHIMS_BINS = {"15-24": (15, 25), "25-34": (25, 35), "35-49": (35, 50),
              "50+": (50, 200), "15-49": (15, 50), "15+": (15, 200)}
PLOT_BINS = ["15-24", "25-34", "35-49", "50+"]

# data/art_coverage.csv's strata. [45,100) is read as 45+ -- no agents above 100.
ART_BINS = {"[15,25)": (15, 25), "[25,35)": (25, 35),
            "[35,45)": (35, 45), "[45,100)": (45, 200)}

STATES = ("n_infected", "n_diagnosed", "n_on_art", "n_effective_art")

# measure -> (numerator, denominator)
MEASURES = {
    "aware_of_all_plhiv":  ("n_diagnosed", "n_infected"),
    "on_art_given_aware":  ("n_on_art", "n_diagnosed"),
    "vls_given_art":       ("n_effective_art", "n_on_art"),
    "on_art_of_all_plhiv": ("n_on_art", "n_infected"),
    "vls_of_all_plhiv":    ("n_effective_art", "n_infected"),
}

# House colours, carried from standard_figures.py. Validated: adjacent-pair
# CVD dE 21.8 (protan), 28.8 normal vision.
COLS = {"f": "#c0392b", "m": "#2c6fbb"}
SEXLAB = {"f": "Women", "m": "Men"}
INK, MUTED = "#222222", "#6b6b6b"


def aggregate(df, bins, group_cols=("seed",)):
    """Sum 5-year bands into `bins`.

    `group_cols` are carried through so a scan can keep (k_m, k_f) alongside
    the seed. Returns long: group_cols + year, sex, bin, and the four counts.
    """
    group_cols = list(group_cols)
    rows = []
    for name, (lo_b, hi_b) in bins.items():
        members = [(lo, hi) for lo, hi in BANDS if lo >= lo_b and hi <= hi_b]
        if not members:
            raise ValueError(f"no 5-year bands inside {name}")
        for sex in ("f", "m"):
            acc = {}
            for state in STATES:
                cols = [f"cascadeage.{state}_{sex}_{lo}_{hi}" for lo, hi in members]
                missing = [c for c in cols if c not in df.columns]
                if missing:
                    raise KeyError(f"missing analyzer columns: {missing[:3]}")
                acc[state] = df[cols].sum(axis=1)
            sub = pd.DataFrame(acc)
            sub["year"] = df["timevec"].values
            for g in group_cols:
                sub[g] = df[g].values
            sub["sex"] = sex
            sub["bin"] = name
            rows.append(sub)
    return pd.concat(rows, ignore_index=True)


def conditionals(long, keys=("year", "sex", "bin")):
    """Seed-pooled cascade rates, plus the across-seed spread.

    `keys` is the explicit grouping -- everything NOT named there (normally just
    `seed`) is pooled over. Pass e.g. ("k_m", "k_f", "year", "sex", "bin") for a
    scan.

    The point estimate pools numerator and denominator across seeds BEFORE
    dividing. Averaging per-seed ratios would bias small strata, where a seed
    with few PLHIV would carry the same weight as one with many.
    """
    keys = list(keys)
    pooled = long.groupby(keys, as_index=False)[list(STATES)].sum()
    out = []
    for measure, (num, den) in MEASURES.items():
        p = pooled.copy()
        p["measure"] = measure
        p["model"] = p[num] / p[den].replace(0, np.nan)
        per = long.copy()
        per["r"] = per[num] / per[den].replace(0, np.nan)
        sd = (per.groupby(keys)["r"]
                 .agg(model_sd="std", n_seeds="count").reset_index())
        merged = (p.merge(sd, on=keys, how="left")
                   .rename(columns={num: "numerator", den: "denominator"}))
        out.append(merged[keys + ["measure", "model", "model_sd", "n_seeds",
                                  "numerator", "denominator"]])
    return pd.concat(out, ignore_index=True)


def load_awareness_targets(repo, year=2021):
    """The 8 SHIMS3 awareness strata, tidy."""
    t = pd.read_csv(repo / "data" / "eswatini_cascade_95s.csv")
    # The file's own `survey` column is the survey NAME; move it aside before
    # `value` takes that slot, or two columns end up called `survey`.
    t = (t[t.sex.isin(["f", "m"]) & t.age.isin(PLOT_BINS)
           & (t.measure == "aware_of_all_plhiv") & (t.year == year)]
         .rename(columns={"survey": "survey_name", "age": "bin",
                          "value": "survey"}))
    return t[["sex", "bin", "measure", "survey", "n_unweighted",
              "small_denominator"]]
