"""Extract the SHIMS3 2021 conditional and overall 95-95-95 cascade by age and sex.

Run `python cascade_construction.py` to regenerate `data/eswatini_cascade_95s.csv`.

Why this file exists
--------------------
Exp 029 needed the cascade decomposed into its three conditional steps by age
and sex. The repo had no tabulated awareness data at all -- `art_coverage.csv`
carries on-ART among *all* PLHIV (unconditional), and `eswatini_vls.csv` carries
suppression-given-ART only as a 15+ aggregate by sex. The first 95 was missing
entirely, which meant the step that turns out to bind hardest had no target.

It was in the SHIMS3 report the whole time, in exactly the bins we wanted.

Source: data/241123_SHIMS_ENG_RR3_Final-1.pdf
  Table 9.1.A (p.79, PDF page 81) -- OVERALL percentages: each of the three
      indicators expressed as a share of all PLHIV.
  Table 9.1.B (p.80, PDF page 82) -- CONDITIONAL percentages: aware | PLHIV,
      on ART | aware, suppressed | on ART. This is the decomposition.

Both are "self-reported and antiretroviral biomarker data" -- i.e. awareness and
treatment status are self-report, corrected upward by a detectable ARV in blood.
Tables 9.2.A/B are the alternative basis (self-report adjusted by VL<200) and are
NOT used here; mixing the two bases would break the multiplicative identity below.

No confidence intervals
-----------------------
Neither table publishes CIs -- only the percentage and the unweighted
denominator. `n_unweighted` is carried so a binomial interval can be formed
downstream if needed, but that would ignore the survey weighting and design
effect, so it is left to the consumer rather than baked in here.

Values in parentheses in the report flag a small/unreliable denominator; those
are carried with `small_denominator=True`. One cell is affected: men 15-24 VLS
among those on treatment, (86.5), n=48.

The identity that validates the extraction
------------------------------------------
Conditional percentages must multiply to the overall ones:

    aware|PLHIV x onART|aware            = onART|PLHIV     (9.1.B -> 9.1.A)
    aware|PLHIV x onART|aware x VLS|onART = VLS|PLHIV      (9.1.B -> 9.1.A)

This is checked at the end and the script raises if any cell is off by more than
0.6 percentage points (the tolerance absorbs the report's own rounding to 1dp).
It is the only real guard against a mis-parse, since the PDF text layer gives a
flat token stream with no table structure.

The number that motivated exp 029
---------------------------------
Men 25-34: 65.0% on ART among all PLHIV -- the value exp 026 built its `art_95`
arm around -- decomposes as 74.8% aware x 86.8% on-ART-given-diagnosed. Two
thirds of that gap is the FIRST 95, not the second.
"""

import re
import fitz
import pandas as pd

SRC = "data/241123_SHIMS_ENG_RR3_Final-1.pdf"
OUT = "data/eswatini_cascade_95s.csv"

PAGE_OVERALL = 80      # 0-indexed; Table 9.1.A
PAGE_CONDITIONAL = 81  # 0-indexed; Table 9.1.B

AGES = ["15-24", "25-34", "35-49", "50+", "15-49", "15+"]

# Block header -> (measure, basis). Order matters: blocks appear in this order
# on the page and are sliced between successive headers.
BLOCKS_OVERALL = [
    ("Diagnosed", "aware_of_all_plhiv"),
    ("On Treatment", "on_art_of_all_plhiv"),
    ("Viral Load Suppression (VLS) on Treatment", "vls_of_all_plhiv"),
]
BLOCKS_CONDITIONAL = [
    ("Diagnosed", "aware_of_all_plhiv"),
    ("On Treatment Among Those Diagnosed", "on_art_given_aware"),
    ("Viral Load Suppression (VLS) Among Those on Treatment", "vls_given_art"),
]

NUM = re.compile(r"\(?([\d]+\.[\d])\)?|^([\d,]+)$")


def _tokens(page_text):
    return [ln.strip() for ln in page_text.split("\n") if ln.strip()]


def _parse_block(toks, start, end):
    """Rows run: age, m_pct, m_n, f_pct, f_n, t_pct, t_n -- in that order."""
    rows = []
    i = start
    while i < end:
        if toks[i] in AGES:
            age = toks[i]
            vals = toks[i + 1:i + 7]
            if len(vals) < 6:
                break
            rec = {"age": age}
            for sex, (p_tok, n_tok) in zip(("m", "f", "all"),
                                           [(vals[0], vals[1]), (vals[2], vals[3]),
                                            (vals[4], vals[5])]):
                small = p_tok.startswith("(")
                pct = float(p_tok.strip("()*"))
                n = int(n_tok.replace(",", ""))
                rec[sex] = (pct / 100.0, n, small)
            rows.append(rec)
            i += 7
        else:
            i += 1
    return rows


def _extract(page_text, blocks):
    toks = _tokens(page_text)
    # Locate each block header, then slice to the next header (or end of page).
    idx = []
    for header, _ in blocks:
        found = next((j for j, t in enumerate(toks) if t == header), None)
        if found is None:
            raise ValueError(f"block header not found: {header!r}")
        idx.append(found)
    idx.append(len(toks))

    out = []
    for k, (header, measure) in enumerate(blocks):
        for rec in _parse_block(toks, idx[k], idx[k + 1]):
            for sex in ("m", "f", "all"):
                pct, n, small = rec[sex]
                out.append(dict(year=2021, survey="SHIMS3 2021", sex=sex,
                                age=rec["age"], measure=measure, value=pct,
                                n_unweighted=n, small_denominator=small,
                                source_table=("Table 9.1.A" if blocks is BLOCKS_OVERALL
                                              else "Table 9.1.B")))
    return out


def build():
    doc = fitz.open(SRC)
    rows = _extract(doc[PAGE_OVERALL].get_text(), BLOCKS_OVERALL)
    rows += _extract(doc[PAGE_CONDITIONAL].get_text(), BLOCKS_CONDITIONAL)
    df = pd.DataFrame(rows)

    # 'aware_of_all_plhiv' is identical in both tables (the first 95 is already
    # conditional on nothing); drop the duplicate from 9.1.B.
    df = df.drop_duplicates(subset=["year", "sex", "age", "measure"], keep="first")
    return df.sort_values(["measure", "sex", "age"]).reset_index(drop=True)


def check(df, tol=0.006):
    """Conditional percentages must multiply to the overall ones."""
    p = df.pivot_table(index=["sex", "age"], columns="measure", values="value")
    problems = []
    lhs = p["aware_of_all_plhiv"] * p["on_art_given_aware"]
    for key, (got, want) in zip(p.index, zip(lhs, p["on_art_of_all_plhiv"])):
        if abs(got - want) > tol:
            problems.append(f"on_art {key}: {got:.4f} vs {want:.4f}")
    lhs2 = lhs * p["vls_given_art"]
    for key, (got, want) in zip(p.index, zip(lhs2, p["vls_of_all_plhiv"])):
        if abs(got - want) > tol:
            problems.append(f"vls {key}: {got:.4f} vs {want:.4f}")
    if problems:
        raise ValueError("cascade identity failed -- likely a mis-parse:\n  "
                         + "\n  ".join(problems))
    return True


if __name__ == "__main__":
    df = build()
    check(df)
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(df)} rows, {df.measure.nunique()} measures")
    print(df.pivot_table(index=["sex", "age"], columns="measure",
                         values="value").round(3).to_string())
