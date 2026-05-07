# Experiment 008 — Calibration targets

*First calibration experiment. Experiments 001–007 were model-development;
008 marks the transition to calibration proper.*

## Question

What exactly are we calibrating to? Produce a single, authoritative
"calibration targets" reference figure showing PHIA prevalence (by age ×
sex × year), UNAIDS deaths, and UNAIDS age distribution — each with its
uncertainty — and audit whether the current data files contain the
uncertainty information we need.

This is the "data + visualisation" step (steps 1–2 of `calibration-workflow`).
No calibration runs yet; no model runs at all.

## Why this first

A calibration is only as defensible as its targets. Before any prior
predictive check, method choice, or sweep, we need:

1. A frozen, named set of targets (file → variable → year/age/sex)
2. Plotted uncertainty for every target (95% CI or equivalent)
3. An explicit audit: where uncertainty info is missing, where files
   disagree, where definitions are ambiguous

If a target lacks defensible uncertainty, that becomes a follow-up task
(e.g. derive design effects from PHIA microdata) before fitting.

## Inputs

Existing files in `hivsim_eswatini/`:

| Target | Likely source file(s) |
|---|---|
| Prevalence by age × sex × year | `calibration_data/prevalence_by_age_sex.csv`; cross-check with `external_data/SwazilandPrevalencePHIA2017_binned.csv`, `external_data/SWAZILAND_nationalprevalence_all_updatedPHIA3.csv` |
| Deaths (UNAIDS) | `data/eswatini_deaths.csv` |
| Age distribution (UNAIDS) | `data/eswatini_age_1985.csv` (1985 baseline?) — confirm whether multi-year age structure exists elsewhere |
| **Validation hold-out (do not load into fitting workflow):** 2021 incidence | `calibration_data/incidence_2021_VALIDATION_ONLY.csv` |

## Plan

1. **Load.** Read each target into a single tidy `pandas` dataframe with
   columns: `quantity`, `year`, `sex`, `age_low`, `age_high`,
   `value`, `lower`, `upper`, `source`, `notes`.
2. **Audit.** For each row, check whether `lower`/`upper` is populated.
   Where missing, flag in the audit table and note what would be needed
   (e.g. "PHIA design effect not in CSV — needs microdata or DHS-style
   design-effect assumption").
3. **Plot.** One multi-panel figure: prevalence panels (faceted by sex,
   age band on x, one line per survey year with shaded CI), deaths over
   time (line + CI band), age distribution (bar + CI). Save to
   `outputs/calibration_targets.png`.
4. **Freeze.** Write `outputs/calibration_targets.csv` (the tidy
   dataframe). This becomes the canonical fitting target file referenced
   by all downstream experiments.
5. **Summarise.** Write `SUMMARY.md` covering: what's in (and what's
   not), uncertainty gaps to fix, any cross-file disagreements, and the
   tidy target file's path.

## Success criteria

- [ ] `outputs/calibration_targets.csv` exists, frozen, with all listed
      quantities present
- [ ] `outputs/calibration_targets.png` shows every target with
      uncertainty rendered (or explicit "uncertainty TBD" annotation)
- [ ] `SUMMARY.md` lists every target row missing uncertainty and the
      remediation plan
- [ ] 2021 incidence is **not** loaded by `run.py` — it's segregated for
      validation only

## Out of scope

- No model runs.
- No likelihood design (deferred to `likelihood-design` skill in
  experiment 03 or so).
- No method choice (deferred to `method-selection`).
- ART coverage, VMMC, condom use are not loaded — even though the files
  exist — because they are not in the fitting target set.

## Status

`README.md` written. `run.py` to be implemented next.
