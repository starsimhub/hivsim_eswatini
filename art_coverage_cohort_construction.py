"""Build data/art_coverage_cohorts.csv: forced ART coverage on exp 027's age cohorts.

Exp 027 step 2d (DTC_check/context_dtc.md, 2026-11-19: "subdivide testing and
treatment rates by age cohorts (0-24), (25-34), (35-49), (50+)"; no 90-90-90
data exist under 15, so 0-14 take the 15-24 values).

`data/art_coverage.csv` (used by experiments up to 027c) has bins [15,25),
[25,35), [35,45), [45,100) and nothing under 15, so children were never given
a forced ART target. This file uses [0,25), [25,35), [35,50), [50,100) -- the
SHIMS 90-90-90 table's own bands -- so the ART input and the cascade targets
describe the same people.

ART coverage among ALL PLHIV = (% diagnosed) x (% on ART among diagnosed),
from DTC_check/SHIMS_Datasets/SHIMS_909090_coverage.csv (SHIMS2 Table 10.3.B,
SHIMS3), by year:

  2004  0 everywhere (as data/art_coverage.csv: no ART before the programme).
  2011  SHIMS1, carried over from data/art_coverage.csv and REMAPPED onto the
        new bands: [15,25)->[0,25), [25,35)->[25,35), [35,45)->[35,50),
        [45,100)->[50,100). An approximation: SHIMS1 is not re-tabulated.
  2016  SHIMS2 (fieldwork 2016-17; labelled 2016 like the other SHIMS2-based
        inputs, so make_interventions' shims2_year moves it to 2017).
        50+ is not published: it is backed out of the 15+ and 15-49 rows,
        weighting by the table's UNWEIGHTED n (an approximation; survey
        weights differ). Men 50+: 88.5% diagnosed, 95.0% of them on ART;
        women 93.2% / 91.0%.
  2021  SHIMS3, all four bands published.

Gender uses stisim's convention: 0 = female, 1 = male.

Run from the repo root:  python art_coverage_cohort_construction.py
"""

import pandas as pd

SRC = "DTC_check/SHIMS_Datasets/SHIMS_909090_coverage.csv"
OLD = "data/art_coverage.csv"
OUT = "data/art_coverage_cohorts.csv"

COHORTS = {"15-24": "[0,25)", "25-34": "[25,35)", "35-49": "[35,50)", "50+": "[50,100)"}
OLD_TO_NEW = {"[15,25)": "[0,25)", "[25,35)": "[25,35)", "[35,45)": "[35,50)", "[45,100)": "[50,100)"}
GENDER = {"f": 0, "m": 1}
SURVEY_YEAR = {2017: 2016, 2021: 2021}  # SHIMS2 rows are written as 2016 (see docstring)


def shims_table():
    s = pd.read_csv(SRC, encoding="utf-8-sig")
    s["sex"] = s.Gender.str.lower().map({"female": "f", "male": "m"})
    s["p"] = s.Percentage / 100
    return s


def backout_50plus(s, year, metric, sex):
    """50+ from the 15+ and 15-49 rows, weighting by unweighted n."""
    g = s[(s.Year == year) & (s.Metric == metric) & (s.sex == sex)].set_index("Age group")
    n_all, n_young = g.Number["15+"], g.Number["15-49"]
    return (n_all * g.p["15+"] - n_young * g.p["15-49"]) / (n_all - n_young)


def build():
    s = shims_table()
    rows = []
    for sex, gcode in GENDER.items():
        for band in COHORTS.values():
            rows.append(dict(Year=2004, Gender=gcode, AgeBin=band, p_art=0.0, source="no ART"))
    old = pd.read_csv(OLD)
    for _, r in old[old.Year == 2011].iterrows():
        rows.append(dict(Year=2011, Gender=r.Gender, AgeBin=OLD_TO_NEW[r.AgeBin],
                         p_art=r.p_art, source=f"SHIMS1 via {OLD} {r.AgeBin}"))
    for year, label in SURVEY_YEAR.items():
        for sex, gcode in GENDER.items():
            g = s[(s.Year == year) & (s.sex == sex)]
            for grp, band in COHORTS.items():
                if grp == "50+" and not ((g["Age group"] == "50+").any()):
                    dx = backout_50plus(s, year, "Diagnosed", sex)
                    art = backout_50plus(s, year, "Diagnosed_OnART", sex)
                    src = "SHIMS2, 50+ backed out of 15+ and 15-49 (unweighted n)"
                else:
                    pick = lambda m: g[(g.Metric == m) & (g["Age group"] == grp)].p.iloc[0]  # noqa: E731
                    dx, art = pick("Diagnosed"), pick("Diagnosed_OnART")
                    src = f"{'SHIMS2' if year == 2017 else 'SHIMS3'} {grp}"
                rows.append(dict(Year=label, Gender=gcode, AgeBin=band,
                                 p_art=round(dx * art, 4), source=src))
    return pd.DataFrame(rows).sort_values(["Year", "Gender", "AgeBin"]).reset_index(drop=True)


if __name__ == "__main__":
    df = build()
    df.drop(columns="source").to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(df)} rows")
    print(df.to_string(index=False))
