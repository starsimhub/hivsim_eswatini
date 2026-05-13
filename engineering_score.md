# Project Engineering Score

- **Project**: `/home/cliffk/idm/hivsim_eswatini`
- **Tier**: 2 (Small-scale project used by multiple people or projects)
- **Overall Score**: 66/100
- **Status**: PASS
- **Date**: 2026-05-13
- **Version**: idm-eng-plugin:eng-quality-checker v1.3_2026.04.13
- **Time spent**: 123s

## Summary

| Category | Score | Weight |
| -- | -- | -- |
| Quality | 56/100 | 40% |
| Usability | 77/100 | 40% |
| Safety | 66/100 | 20% |
| **Total** | **66/100** | 100% |

| Metric | Score | Notes |
| -- | -- | -- |
| correct | 5/10 | No tests anywhere; correctness depends on visual inspection of calibration figures. |
| clear | 7/10 | Good modular split with docstrings, but a few public entry points lack docstrings and `plot_calibration.py` runs sims at import time. |
| concise | 7/10 | Mostly lean and library-driven, but `plot_calibration.py` and `plot_incidence.py` duplicate `run_one`/`summarize` helpers; some dead commented code. |
| simple | 8/10 | README provides a clear numbered workflow and `make_sim()` has sensible defaults; scripts depend on hardcoded relative paths with no validation. |
| powerful | 7/10 | Calibration parameters and key sim args are exposed, but many model assumptions (network proportions, beta, condom efficacy) remain hardcoded inside `make_sim()`. |
| performant | 8/10 | Vectorized agent step logic and parallelized `to_df()`/network data runs; minor non-vectorized loop in `NetworkSnapshot._capture_snapshot` over 36 age bins. |
| documented | 8/10 | Strong README with workflow, file table, key concepts; major classes documented; no tutorials/notebooks and several helpers lack docstrings. |
| accessible | 7/10 | Public on GitHub with MIT LICENSE, but install requires cloning three separate repos (>3 commands) and `raw_data/` is not redistributable. |
| compliant | 9/10 | MIT license, no exposed secrets, all dependencies permissively licensed; raw sensitive data correctly excluded via `.gitignore`. |
| reproducible | 3/10 | No `requirements.txt`/`pyproject.toml` and no version pins for starsim/stisim/sciris; seeds set consistently but environment cannot be reliably reconstructed. |

The project is a well-organized small-scale HIV transmission model for Eswatini with a strong README, sensible modular layout, and clean MIT-licensed code. The biggest gaps are the complete absence of automated tests (driving the `correct` score down) and the lack of any dependency-specification file with pinned versions (dragging `reproducible` to 3/10). Fixing reproducibility and adding a minimal test suite would move this project comfortably into the 80s.

## Recommendations

1. **[reproducible] — Add a `pyproject.toml` with pinned dependencies** *(effort: quick; automated: yes)*
   Create a `pyproject.toml` (or at minimum `requirements.txt`) that lists `numpy`, `pandas`, `sciris`, `matplotlib`, `starsim`, and `stisim` with `>=` version pins matching what was used to produce the published calibration figures. Document the tested Python version. This is the single highest-impact fix (`weight=4`, score gap of 7).

2. **[correct] — Add a minimal `tests/` directory with pytest tests for the main workflows** *(effort: medium; automated: partial)*
   Create `tests/test_smoke.py` that runs `make_sim()` with a tiny population (e.g., `n_agents=500`, short horizon) and asserts results structure; add `tests/test_calibration.py` exercising `make_calibration()` with `n_trials=1`; add a test verifying that `hiv_epi` analyzer denominators (susceptible counts) are non-negative and that prevalence stays in [0, 1]. Wire it into a GitHub Action so PRs run tests. Highest impact on the `correct` metric (`weight=7`).

3. **[concise + clear] — Factor `run_one` and `summarize` helpers into `utils.py`** *(effort: quick; automated: yes)*
   Move the duplicated `run_one(seed, ...)` and `summarize(stack)` functions out of `plot_calibration.py` and `plot_incidence.py` into a single shared location in `utils.py`. Import from both plot scripts. Also delete the commented-out `vmmc`/`n_vmmc` lines in `interventions.py`.

4. **[clear] — Move module-level simulation calls inside `if __name__ == "__main__":` guards** *(effort: quick; automated: yes)*
   `plot_calibration.py` lines 31-33 execute 10 simulations on import. Wrap the script body in a `main()` function called only when run as `__main__` so the module can be imported without side effects.

5. **[clear] — Add docstrings to remaining public functions** *(effort: quick; automated: partial)*
   Add Google-style docstrings to `make_interventions()` (interventions.py), `run_msim()` (run_msim.py), `check_hiv_alive`, `prune_columns`, and `plot_panel`. Each should document arguments, returns, and a one-line purpose.

6. **[correct] — Replace magic numbers in `plot_incidence.py` with cited constants** *(effort: quick; automated: no)*
   The PHIA reference values `0.0314`, `0.0165`, `0.0173`, `0.0085` on lines 19-22 should be moved to a named constant block (e.g., `PHIA_INCIDENCE_2021 = {...}`) with an inline comment citing the PHIA report version and table number.

7. **[simple] — Make scripts robust to working directory** *(effort: quick; automated: yes)*
   Replace hardcoded relative paths (`raw_data/...`, `calibration_data/...`, `results/...`) with paths anchored to `Path(__file__).parent`. Add a friendly `FileNotFoundError` message pointing users to the section of the README that explains how to get the raw data.

8. **[powerful] — Expose network and transmission assumptions as `make_sim()` kwargs** *(effort: medium; automated: no)*
   Promote hardcoded values (network proportions, `beta`, condom efficacy, age structure) to keyword arguments of `make_sim()` with sensible defaults so users can run scenarios without editing the function body.

9. **[documented] — Add a runnable tutorial notebook or `examples/` script** *(effort: medium; automated: no)*
   Create `examples/quickstart.py` (or a Jupyter notebook) that walks through: build a sim, run it, plot prevalence by age, and run a tiny calibration. Reference it from the README.

10. **[accessible] — Add a project-root `CHANGELOG.md` and `CONTRIBUTING.md`** *(effort: quick; automated: no)*
    The existing `figures/CHANGELOG.md` is figure-specific. Add a top-level `CHANGELOG.md` summarizing model version changes, and a short `CONTRIBUTING.md` describing how to set up the dev environment and run tests.

11. **[performant] — Vectorize the per-age-bin loop in `NetworkSnapshot._capture_snapshot`** *(effort: quick; automated: yes)*
    Replace the Python loop over 36 age bins with a `numpy.bincount`/`groupby` approach. Low priority since this only runs at analyzer checkpoints.

12. **[compliant] — Document data provenance in a `data/README.md`** *(effort: quick; automated: no)*
    Briefly note source (UNAIDS, PHIA), version/year, and link to the public dataset for each file in `calibration_data/` and `external_data/`. Resolves the minor "data provenance" uncertainty noted by the safety scorer.

## Full Results

```yaml
project: /home/cliffk/idm/hivsim_eswatini
tier: 2
overall_score: 66
failed: false
quality:
  correct:
    score: 5
    weight: 7
    reason: "No test files exist anywhere in the project (no test_*.py, tests/ directory, or CI/CD); the scientific logic in analyzers.py and utils.py appears correct (incidence denominator uses susceptible counts, prevalence uses np.mean over boolean array), but correctness relies entirely on visual inspection of calibration figures rather than automated validation. Magic numbers in plot_incidence.py lines 19-22 (e.g., REF_DATA with hardcoded values 0.0314, 0.0165) have no inline citation to their source."
  clear:
    score: 7
    weight: 2
    reason: "Good modular structure with separate files for interventions, analyzers, utilities, run scripts, and plot scripts, and most functions have docstrings; however, make_interventions() in interventions.py has no docstring, run_msim() in run_msim.py has no docstring, and plot_calibration.py has module-level simulation code (lines 31-33) that executes on import, violating the module boundary pattern used consistently elsewhere."
  concise:
    score: 7
    weight: 1
    reason: "plot_calibration.py and plot_incidence.py each independently define near-identical run_one() and summarize() helper functions that could be shared; commented-out dead code remains in interventions.py (n_vmmc, vmmc lines). Otherwise the codebase makes good use of numpy, pandas, and sciris without hand-rolling common operations."
usability:
  simple:
    score: 8
    weight: 3
    reason: "The README provides a clear step-by-step workflow (5 numbered steps) with a file descriptions table and the central make_sim() function has sensible defaults (seed, start, stop, verbose). However, plot_calibration.py runs 10 sims at module-level import time (lines 32-33) rather than inside a guard, and many scripts use hardcoded relative paths that break if not run from the repo root, with no validation or error messages to guide the user."
  powerful:
    score: 7
    weight: 2
    reason: "make_sim() exposes key parameters (seed, start, stop, analyzers) and make_calibration() exposes calibration ranges via calib_pars. However, many model assumptions (network proportions, beta, condom efficacy) are hardcoded inside make_sim() rather than exposed as arguments, requiring code edits to explore alternative parameterizations or scenarios."
  performant:
    score: 8
    weight: 2
    reason: "The agent-level step logic in analyzers.py uses numpy boolean masking rather than Python loops, run_msim.py parallelizes to_df() via sc.parallelize, and run_network_data.py uses ss.parallel. The NetworkSnapshot._capture_snapshot iterates over 36 individual age bins in a Python loop rather than vectorizing, a minor but observable inefficiency; no performance tests or profiling infrastructure is present."
  documented:
    score: 8
    weight: 2
    reason: "The README is comprehensive — it covers purpose, installation, a full numbered workflow, a file-description table, key concepts, and git hygiene. Major classes (hiv_epi, NetworkSnapshot) and most public functions have docstrings. However, there are no tutorials or notebooks, the hivsim_swz_docs.md is sparse (mostly meeting notes), and several helper functions (check_hiv_alive, prune_columns, plot_panel) lack docstrings."
  accessible:
    score: 7
    weight: 1
    reason: "The repo is public on GitHub under starsimhub/hivsim_eswatini with an MIT LICENSE. Installation requires cloning two external dependencies (starsim, stisim) plus this repo — four non-trivial steps, more than the ideal 1-3 commands — and critical raw_data/ is not tracked in git (users must ask the project lead), making it impossible to fully reproduce results from the public repo alone. No CHANGELOG at the repo root, no contributing guidelines, and no pyproject.toml/setup.py."
safety:
  compliant:
    score: 9
    weight: 6
    reason: "MIT license is present, no exposed secrets or PII found, and dependencies (numpy, pandas, sciris, starsim, stisim, matplotlib) are all permissive (MIT/BSD). Raw sensitive data is excluded via .gitignore; the only minor uncertainty is that external data provenance (UNAIDS, PHIA) is described informally in the README without formal data-use agreements documented in the repo."
  reproducible:
    score: 3
    weight: 4
    reason: "No dependency specification file (requirements.txt, pyproject.toml, etc.) exists — the five key dependencies (starsim, stisim, sciris, numpy, pandas) are installed via manual git clone at unspecified versions, making environment reconstruction unreliable. Random seeds are consistently set (seed=1 throughout), and git version history exists, but there are no semantic version tags and no lock file or pinned versions."
```
