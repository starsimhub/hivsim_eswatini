"""Import `historymatching` on Python 3.14, where TensorFlow cannot be installed.

Use instead of importing the package directly:

    from hm_shim import hm          # equivalent to `import historymatching as hm`

The problem
-----------
`historymatching` 2.0.1 declares `tensorflow>=2.18`, `tf-keras` and
`gpflow>=2.5` as hard dependencies, and **TensorFlow publishes no wheels for
Python 3.14** -- the interpreter this project runs on (`C:\\Python314`). So
`pip install` fails outright.

Installing with `--no-deps` gets the package in, but importing it still fails:
`historymatching/emulators/__init__.py` imports every emulator eagerly,
including `GPR`, and `emulators/gpr.py` line 7 is a bare `import gpflow`.

The fix, and why it is safe
---------------------------
`gpflow` is referenced **only inside GPR's methods** -- `gpflow.models.GPR`,
`gpflow.kernels.SquaredExponential`, `gpflow.Parameter`,
`gpflow.optimizers.Scipy` and `gpflow.utilities` all appear in method bodies
(gpr.py lines 81-99 and 194), never at class-definition time. So the top-level
import is the *only* thing that needs to exist for the module to load.

This installs a stub `gpflow` (and `tensorflow`) into `sys.modules` before
importing the package. `BayesLinear` -- the emulator this skill recommends as
first choice, and pure NumPy/SciPy -- then works completely. Any attempt to
actually *use* `GPR` raises a clear error naming this file rather than failing
obscurely deep inside gpflow.

Why a shim and not a patch to site-packages
-------------------------------------------
Editing the installed package would work and would be silently destroyed by the
next `pip install --upgrade`. That is exactly the failure this project already
recorded once: the exp-005 VMMC fix was patched into the editable stisim
checkout and wiped by a `git pull` (see exp 015, and the note in CLAUDE.md).
In-repo shims survive reinstalls; patches to site-packages do not.

Worth reporting upstream
------------------------
Moving `import gpflow` inside GPR's methods (or wrapping the `from .gpr import
GPR` line in a try/except) would make the whole package importable without
TensorFlow, which matters for any environment on a Python version TF has not
caught up to yet. Filed as a note in exp 024's SUMMARY.

If GPR is ever genuinely needed, the options are a Python 3.12 environment or
raccoon (which runs uv-managed 3.12 per CLAUDE.md).
"""

import sys
import types


class _MissingGpflow(types.ModuleType):
    """Stands in for gpflow so `historymatching` can import. Not usable."""

    _MSG = (
        "gpflow/tensorflow are not installed -- TensorFlow has no wheels for "
        "Python 3.14. The 'gpr' emulator is therefore unavailable. Use "
        "'bayes_linear' (pure NumPy/SciPy, and the recommended first choice) "
        "or run on a Python 3.12 environment. See hm_shim.py."
    )

    def __getattr__(self, name):
        raise ImportError(f"{self._MSG}\n(tried to access gpflow.{name})")


def _install_stubs():
    """Stub ONLY gpflow -- never tensorflow.

    An earlier version stubbed both, which broke every diagnostics plot.
    matplotlib's `cbook._is_tensorflow_array` does
    `sys.modules.get("tensorflow").is_tensor` on arrays it is asked to plot;
    with a stub present that lookup succeeds and then raises, so pairplot,
    convergence, z-score and constrained-dims plots all failed. With no
    tensorflow entry the lookup returns None and matplotlib skips the check.

    `historymatching` never imports tensorflow directly -- only gpflow does,
    and only inside GPR's methods.
    """
    if "gpflow" not in sys.modules:
        try:
            __import__("gpflow")          # use the real thing if it exists
        except ImportError:
            sys.modules["gpflow"] = _MissingGpflow("gpflow")


_install_stubs()

import historymatching as hm  # noqa: E402

# The emulator that works in this environment. Referenced by experiment scripts
# so the constraint is visible at the call site rather than buried here.
USABLE_EMULATORS = ("bayes_linear", "linear", "glm")
DEFAULT_EMULATOR = "bayes_linear"

__all__ = ["hm", "USABLE_EMULATORS", "DEFAULT_EMULATOR"]
