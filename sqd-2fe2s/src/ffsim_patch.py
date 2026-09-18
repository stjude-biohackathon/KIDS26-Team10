"""Work around an upstream ffsim bug that blocks the tutorial's ``optimize=True``.

``ffsim.linalg.util.unitaries_to_parameters`` calls ``scipy.linalg.logm(mats)`` on a
batch of shape ``(n_mats, dim, dim)``, but ``scipy.linalg.logm`` only accepts a single
square matrix, so the call raises

    ValueError: expected square array_like input

Every path that round-trips an operator through its parameter vector hits this:
``UCJOp*.to_parameters()`` and, through it, ``from_t_amplitudes(..., optimize=True)``.
Reproduced on ffsim 0.0.83 and 0.0.84.

The fix is to take the matrix logarithm of each matrix in the batch and stack the
results, which is what the docstring and every caller already assume.  Importing this
module patches the function in place; call :func:`apply` explicitly if you prefer.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg

import ffsim.linalg.util as _util

_PATCHED_FLAG = "_sqd_hackathon_batched_logm_patch"


def _unitaries_to_parameters(mats: np.ndarray, real: bool = False) -> np.ndarray:
    mats = np.asarray(mats)
    logs = np.stack([scipy.linalg.logm(m) for m in mats])
    return _util.antihermitians_to_parameters(logs, real=real)


def apply() -> bool:
    """Patch ffsim in place.  Returns True if a patch was installed."""
    if getattr(_util.unitaries_to_parameters, _PATCHED_FLAG, False):
        return False
    try:
        _util.unitaries_to_parameters(np.eye(2)[None, :])
        return False  # upstream is fixed, nothing to do
    except ValueError:
        pass
    setattr(_unitaries_to_parameters, _PATCHED_FLAG, True)
    _util.unitaries_to_parameters = _unitaries_to_parameters
    # the variational modules bind the symbol at import time
    for mod in (
        "ffsim.variational.ucj_spin_balanced",
        "ffsim.variational.ucj_spin_unbalanced",
        "ffsim.variational.ucj_spinless",
        "ffsim.variational.util",
    ):
        try:
            module = __import__(mod, fromlist=["_"])
        except ImportError:
            continue
        if hasattr(module, "unitaries_to_parameters"):
            module.unitaries_to_parameters = _unitaries_to_parameters
        if hasattr(module, "unitary_to_parameters"):
            module.unitary_to_parameters = lambda m, real=False: _unitaries_to_parameters(
                np.asarray(m)[None, :], real=real
            )
    _util.unitary_to_parameters = lambda m, real=False: _unitaries_to_parameters(
        np.asarray(m)[None, :], real=real
    )
    return True


PATCHED = apply()
