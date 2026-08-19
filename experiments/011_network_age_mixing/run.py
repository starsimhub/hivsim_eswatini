"""Experiment 011 — network age-mixing fix: entry point.

This is a model-development / debugging experiment. The real work is split
across three scripts, run in this order:

  1. test_rank_matching.py        — the diagnostic ensemble (5 seeds x 3 configs)
                                     used to confirm the bug and produce the
                                     per-variant realized-gap CSVs/plots.
  2. run_pr477_comparison.py       — orchestrates the algorithm comparison:
                                     switches the local stisim checkout between
                                     the released baseline, the local Gaussian
                                     patch, and PR #477 (origin/459) across a
                                     max_deviation sweep, running the diagnostic
                                     for each and logging to
                                     outputs/comparison_run.log.
  3. plot_algorithm_comparison.py  — builds the consolidated comparison figures
                                     in figures/.

See README.md for the plan and SUMMARY.md for the result. Configuration
(variants, seeds, n_agents, stisim refs) is recorded in config.yaml.

Run the comparison directly:

    python run_pr477_comparison.py
"""

import runpy
import pathlib

if __name__ == "__main__":
    here = pathlib.Path(__file__).parent
    runpy.run_path(str(here / "run_pr477_comparison.py"), run_name="__main__")
