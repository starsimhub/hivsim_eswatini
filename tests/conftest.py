"""
Shared pytest fixtures and path setup for the hivsim_eswatini test suite.

Adds the repo root to ``sys.path`` so test modules can import the project's
top-level modules (``run_sims``, ``analyzers``, etc.) without an editable install.
"""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
