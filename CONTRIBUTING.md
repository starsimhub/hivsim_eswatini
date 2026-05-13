# Contributing

Thanks for contributing to `hivsim_eswatini`. This project is a small research codebase; the workflow is light.

## Dev setup

1. Install the dependencies. The simplest path:

   ```bash
   # In a fresh conda/virtualenv
   git clone https://github.com/starsimhub/starsim.git && pip install -e ./starsim
   git clone https://github.com/starsimhub/stisim.git  && pip install -e ./stisim
   git clone https://github.com/starsimhub/hivsim_eswatini.git
   cd hivsim_eswatini
   pip install -e ".[test]"
   ```

2. Get the raw input data: place a `raw_data/` directory in the repo root (ask the
   project lead — it is not redistributable). The processed inputs in `data/` are
   sufficient for running sims and tests; `raw_data/` is only needed when you
   regenerate `data/eswatini_hiv_calib.csv` via `python utils.py`.

## Running the tests

```bash
pytest
```

The smoke test (`tests/test_smoke.py`) runs a tiny sim to check the build/run path.
Analyzer-invariant tests check that prevalence stays in `[0, 1]` and that
incidence denominators are non-negative.

## Style

- Google-style docstrings on public functions/classes.
- Prefer `pathlib.Path` and `Path(__file__).parent` for file paths so scripts work
  regardless of the working directory.
- Keep the public surface of `make_sim()` minimal — push model-knob exposure
  through new keyword arguments rather than monkey-patching.

## Branches and PRs

- Work on a feature branch off `main`.
- Open a PR with a short description of *why* — the diff already shows *what*.
- The CI suite is minimal; please run `pytest` locally before requesting review.

## Releases

Tag with `v<major>.<minor>.<patch>` (semver). Update `CHANGELOG.md` in the same commit.
