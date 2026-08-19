"""
Orchestrate unattended algorithm comparison between:

    Variant E   — your Gaussian-acceptance patch on stisim main (v1.5.5)
    Variant A   — PR #477 stock (max_deviation=1)
    Variant B   — PR #477 with max_deviation=3
    Variant C   — PR #477 with max_deviation=5
    Variant D   — PR #477 with DHS Eswatini 2006-07 age_diff_pars

Each variant runs test_rank_matching.py in a fresh subprocess (mandatory:
Python won't re-import stisim across git checkouts of the editable install).
All output goes to outputs/ with a variant-label suffix. A consolidated
comparison plot is built at the end. Status is logged continuously to
outputs/comparison_run.log.

The script always restores stisim to main + Gaussian patch in its finally
clause, so a crash mid-run still leaves the working tree in a sane state.
"""

from __future__ import annotations

import datetime
import json
import os
import re
import subprocess
import sys
import shutil
import traceback
from pathlib import Path


EXP_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXP_DIR.parents[1]
OUT_DIR = EXP_DIR / "outputs"
OUT_DIR.mkdir(exist_ok=True)

STISIM_REPO = Path("c:/Users/adamak/OneDrive - Gates Foundation/Dropbox/star_sim/stisim")
NETWORKS_PY = STISIM_REPO / "stisim" / "networks.py"
PYTHON = r"C:\Python314\python.exe"

LOG_PATH = OUT_DIR / "comparison_run.log"


def log(msg: str) -> None:
    ts = datetime.datetime.now().isoformat(timespec="seconds")
    line = f"[{ts}] {msg}\n"
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line)
    print(line, end="", flush=True)


def run_git(args: list[str], cwd: Path = STISIM_REPO, check: bool = True) -> subprocess.CompletedProcess:
    cmd = ["git"] + args
    log(f"(in {cwd.name}) git {' '.join(args)}")
    r = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    if r.stdout.strip():
        log(f"  stdout: {r.stdout.strip()[:400]}")
    if r.returncode and r.stderr.strip():
        log(f"  stderr: {r.stderr.strip()[:400]}")
    if check and r.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed exit {r.returncode}")
    return r


def patch_max_deviation(value: int) -> None:
    """Edit the literal `max_deviation=N` in the PR-#477 call site of networks.py."""
    txt = NETWORKS_PY.read_text(encoding="utf-8")
    new = re.sub(r"max_deviation=\d+", f"max_deviation={value}", txt)
    if new == txt:
        log(f"WARNING: max_deviation literal not found in networks.py — value {value} not applied")
    else:
        NETWORKS_PY.write_text(new, encoding="utf-8")
        log(f"Patched max_deviation = {value} in stisim/networks.py")


def run_diagnostic(label: str, configs_override: dict | None = None) -> bool:
    """Run test_rank_matching.py as a subprocess; outputs land with `_{label}` suffix."""
    log(f"=== Variant {label}: starting diagnostic ===")
    env = os.environ.copy()
    env["EXP011_LABEL"] = label
    if configs_override is not None:
        env["EXP011_CONFIGS_JSON"] = json.dumps(configs_override)

    cmd = [PYTHON, str(EXP_DIR / "test_rank_matching.py")]
    t0 = datetime.datetime.now()
    proc = subprocess.run(cmd, env=env, cwd=str(REPO_ROOT),
                          capture_output=True, text=True)
    dt = (datetime.datetime.now() - t0).total_seconds()
    ok = proc.returncode == 0
    log(f"Variant {label} exited {proc.returncode} after {dt:.0f}s")
    if proc.stdout:
        for line in proc.stdout.splitlines()[-15:]:
            log(f"  out: {line}")
    if not ok and proc.stderr:
        for line in proc.stderr.splitlines()[-25:]:
            log(f"  err: {line}")
    return ok


# DHS-derived age_diff_pars for variant D, from experiments/011/README.md
# (level 0 = marital → larger gaps; level 2 = casual → smaller).
DHS_PARS = {
    "DHS_eswatini": {
        "teens": [(9, 5), (8, 5), (7, 5)],
        "young": [(8, 5), (7, 5), (6, 5)],
        "adult": [(8, 6), (7, 6), (6, 5)],
    }
}


def main() -> int:
    LOG_PATH.write_text("")  # truncate prior log
    log("Starting PR #477 vs Gaussian-patch comparison")
    log(f"Stisim repo: {STISIM_REPO}")
    log(f"Eswatini repo: {REPO_ROOT}")

    # Snapshot starting state so we can always get back
    stash_msg = "pr477-test: park Gaussian patch during orchestration"
    on_pr477_branch = False
    stashed = False

    # Save the existing pre-fix baseline CSV (only if not already saved)
    baseline = OUT_DIR / "rank_test_realized_gaps.csv"
    baseline_backup = OUT_DIR / "rank_test_realized_gaps_baseline_no_fix.csv"
    if baseline.exists() and not baseline_backup.exists():
        shutil.copy(baseline, baseline_backup)
        log(f"Backed up pre-fix baseline to {baseline_backup.name}")

    try:
        # === Variant E: Gaussian patch (current working tree) ===
        # Verify we're on stisim main with the Gaussian patch applied
        run_git(["status", "-s"])
        run_diagnostic("E_gaussian_patch")

        # === Switch to PR #477 ===
        log("--- Switching stisim to PR #477 branch ---")
        run_git(["stash", "push", "-m", stash_msg, "--", "stisim/networks.py"], check=False)
        stashed = True
        # Confirm fetched
        run_git(["fetch", "origin", "459"], check=False)
        # Create or reset our local test branch from origin/459
        # Delete first if it exists (idempotent re-runs)
        existing = run_git(["branch", "--list", "test/pr-477"], check=False).stdout.strip()
        if existing:
            run_git(["checkout", "main"], check=False)
            run_git(["branch", "-D", "test/pr-477"], check=False)
        run_git(["checkout", "-b", "test/pr-477", "origin/459"])
        on_pr477_branch = True
        log(f"On PR #477 branch (HEAD: {run_git(['rev-parse', '--short', 'HEAD']).stdout.strip()})")

        # === Variant A: PR #477 stock (max_deviation=1) ===
        run_diagnostic("A_pr477_md1")

        # === Variant B: max_deviation=3 ===
        patch_max_deviation(3)
        run_diagnostic("B_pr477_md3")

        # === Variant C: max_deviation=5 ===
        patch_max_deviation(5)
        run_diagnostic("C_pr477_md5")

        # Restore PR #477 source before next variant
        run_git(["checkout", "--", "stisim/networks.py"])

        # === Variant D: PR #477 stock with DHS-derived pars ===
        run_diagnostic("D_pr477_dhs", configs_override=DHS_PARS)

    except Exception:
        log(f"FATAL during variant runs:\n{traceback.format_exc()}")

    finally:
        # === Restore stisim to main + Gaussian patch ===
        log("--- Restoring stisim main + Gaussian patch ---")
        try:
            if on_pr477_branch:
                run_git(["checkout", "--", "stisim/networks.py"], check=False)
                run_git(["checkout", "main"], check=False)
                run_git(["branch", "-D", "test/pr-477"], check=False)
            if stashed:
                run_git(["stash", "pop"], check=False)
            log("Stisim restoration complete")
        except Exception:
            log(f"WARNING during stisim restoration:\n{traceback.format_exc()}")

        # === Comparison plot ===
        log("--- Building consolidated comparison plot ---")
        try:
            plot_cmd = [PYTHON, str(EXP_DIR / "plot_algorithm_comparison.py")]
            r = subprocess.run(plot_cmd, capture_output=True, text=True)
            for line in (r.stdout + r.stderr).splitlines()[-15:]:
                log(f"  plot: {line}")
            log(f"Plot subprocess exit {r.returncode}")
        except Exception:
            log(f"WARNING during plot:\n{traceback.format_exc()}")

        log("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
