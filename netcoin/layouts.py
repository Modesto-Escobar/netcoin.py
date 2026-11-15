"""Node coordinates computed from coincidences (ports of netCoin's layouts)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .coincidence import Coin


def _fix_signs(coords: np.ndarray) -> np.ndarray:
    """Make the sign of every axis deterministic (the largest coordinate is positive).

    The sign of an eigenvector or singular vector is arbitrary, and R and numpy
    pick different ones; a fixed convention makes the layout reproducible.
    """
    coords = coords.copy()
    for j in range(coords.shape[1]):
        if coords[np.argmax(np.abs(coords[:, j])), j] < 0:
            coords[:, j] = -coords[:, j]
    return coords


def layout_mca(incidences: pd.DataFrame) -> pd.DataFrame:
    """Coordinates of the variables from a simple correspondence analysis.

    Port of ``layoutMCA``: classical correspondence analysis of the table of
    incidences completed with its complement. Returns principal coordinates
    ``F1`` and ``F2``, one row per variable.
    """
    x = incidences.to_numpy(dtype=float)
    full = np.column_stack([x, 1 - x])
    p = full / full.sum()
    col_mass = p.sum(axis=0)
    if (col_mass == 0).any():
        raise ValueError("a variable is present in every case or in none: "
                         "the correspondence analysis is not defined")
    row_mass = p.sum(axis=1)
    expected = np.outer(row_mass, col_mass)
    z = (p - expected) / np.sqrt(expected)
    _, d, vt = np.linalg.svd(z, full_matrices=False)
    k = x.shape[1]
    standard = vt.T[:k, :2] / np.sqrt(col_mass[:k])[:, None]
    principal = _fix_signs(standard * d[:2])
    return pd.DataFrame(principal, index=list(incidences.columns), columns=["F1", "F2"])


def layout_pca(coincidences: Coin | pd.DataFrame) -> pd.DataFrame:
    """Coordinates of the variables from the eigenvectors of their Pearson correlations.

    Port of ``layoutPCA``. Takes a :class:`Coin` (or a correlation matrix) and
    returns ``F1`` and ``F2``, one row per variable.
    """
    if isinstance(coincidences, Coin):
        from .similarity import sim

        matrix = sim(coincidences, "Pearson")
    else:
        matrix = pd.DataFrame(coincidences)
    values, vectors = np.linalg.eigh(matrix.to_numpy(dtype=float))
    order = np.argsort(values)[::-1][:2]
    coords = vectors[:, order] * np.sqrt(np.maximum(values[order], 0))
    return pd.DataFrame(_fix_signs(coords), index=matrix.index, columns=["F1", "F2"])
