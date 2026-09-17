# Exp 029 — Where does the treatment cascade actually stand, step by conditional step, by age and sex

**Question.** Before drawing conclusions from 026's scenario ordering or 028's
identifiability verdict, what does the model's care cascade look like decomposed
into its three conditional steps — **diagnosed | infected**, **on ART |
diagnosed**, **suppressed | on ART** — for each of 15-24, 25-34, 35-49 by sex?

026 concluded the cascade is the larger prize and built its `art_95` arm as a
*coverage* push, on the strength of an aggregate observation (men 25-34 at 0.650
ART coverage). That arm is only meaningful if the model can actually reach 95%.
In stisim, `ART` with a coverage target fills a stratum **only from
`diagnosed & ~on_art` agents in that stratum** (`hiv_interventions.py:573`).
If testing has not diagnosed enough people, the target is silently unmet — the
run does not warn. So the arm 026 costed could be measuring an impossibility,
and no existing output can tell us either way.

**Why this cannot be answered from what we already have.** The `Cascade`
analyzer (`analyzers.py:237`) records `infected`, `on_art`, `on_effective_art`
and `on_nonsupp_art` — but **not `diagnosed`** — and only for 15-49 and 15+ by
sex. There is no age-banded cascade output anywhere in `experiments/`. The
per-agent file `026/outputs/capacity_agents_baseline.csv` has diagnosis and ART
*timings* but carries neither age nor sex, so it cannot be re-binned.

**And it cannot be fully answered from data either.** The conditional
decomposition is not observed at this stratification:

| step | observed by age x sex? |
|---|---|
| diagnosed \| infected | **no** — no tabulated awareness data in the repo; only inside `data/241123_SHIMS_ENG_RR3_Final-1.pdf` |
| on ART \| diagnosed | **no** — `art_coverage_by_age_sex.csv` is on-ART \| **PLHIV**, unconditional |
| suppressed \| on ART | **no** — `eswatini_vls.csv` has `vls_given_art` only as a 15+ aggregate by sex (f 0.959, m 0.967) |

What *is* observed by age and sex is the two unconditional levels: on-ART | PLHIV
(3 bins x sex x 2011/2016/2021) and VLS | PLHIV (5-year bands x sex, 2016/2021).
So this experiment reports the conditional decomposition as a **model
diagnostic**, anchored against the two unconditional levels that have data. The
first 95 is model output with no target to check it against — that is a finding
about the target set, not a defect of the run, and it is recorded as such.

**Plan.** One parameter point — 024 design row 868, the draw that fit 48/48
targets, the same point 026 used — at model-v1.5, 10 seeds, local. New analyzer
`CascadeByAge` recording per (3 age bins x {f,m} x year): `n_infected`,
`n_diagnosed`, `n_on_art`, `n_effective_art`, plus the three conditionals and
the two unconditionals. Reported at 2011, 2016, 2021 (data years) and 2030
(026's scenario horizon).

Additionally, instrument the ART stratum loop to record **target vs. achieved**
per stratum per year, so an unmet coverage target becomes visible rather than
silent.

**Two caveats to carry into the reading.**

1. *ART coverage is an input, not an emergent quantity.* The model is driven to
   observed coverage via `art_coverage`, so `on ART | PLHIV` is pinned to data
   by construction wherever the target is attainable. The informative quantities
   are therefore the **first step** (diagnosed | infected, which is free) and the
   **shortfall** (strata where the target could not be met).
2. *The 35-49 bin is not clean.* The ART coverage input is keyed on `[35,45)` and
   `[45,100)`, so the model's 45-49 slice is driven by the 45+ band — which
   `art_coverage_construction.py` documents as an extrapolation for 2011 and 2016
   and as 50+-weighted for 2021. The 35-49 row is reported as asked, with the
   `[35,45)` sub-band alongside it.

**Success criteria.** A good result: the first 95 is plausible (roughly 80-95%
by 2021, lowest in young men), no stratum shows a material unmet ART target, and
the two unconditional levels track the data within a few points. A bad result —
and the one that would matter most — is a large silent shortfall in men 25-34,
which would mean 026's `art_95` arm was partly unachievable and its headline
contrast needs re-reading.
