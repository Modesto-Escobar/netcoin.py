"""Rotations of the components of a principal components analysis (``varimax``, ``quartimax``, ``equamax``)."""
from __future__ import annotations

import warnings
from itertools import combinations

import numpy as np
import pandas as pd

ROTATIONS = ("varimax", "quartimax", "equamax")


def rotation_name(x):
    """Lower case, and ``equimax`` read as ``equamax``."""
    x = [str(v).lower() for v in (x if isinstance(x, (list, tuple, np.ndarray)) else [x])]
    return ["equamax" if v == "equimax" else v for v in x]


def varimax(x: np.ndarray, normalize: bool = True, eps: float = 1e-5):
    """Varimax rotation: ``(rotated loadings, rotation matrix)`` (as ``stats::varimax``)."""
    p, k = x.shape
    if k < 2:
        return x, np.eye(k)
    scale = np.ones(p)
    if normalize:
        scale = np.sqrt((x ** 2).sum(axis=1))
        x = x / scale[:, None]
    t, d = np.eye(k), 0.0
    for _ in range(1000):
        z = x @ t
        b = x.T @ (z ** 3 - z @ np.diag((z ** 2).sum(axis=0)) / p)
        u, s, vt = np.linalg.svd(b)
        t = u @ vt
        previous, d = d, s.sum()
        if d < previous * (1 + eps):
            break
    return (x @ t) * scale[:, None], t


def _gpforth(a: np.ndarray, criterion, normalize: bool = True, eps: float = 1e-5, maxit: int = 5000):
    """Orthogonal rotation by gradient projection (GPArotation's ``GPForth``) of the loadings ``a``."""
    p, k = a.shape
    weight = np.ones(p)
    if normalize:
        weight = np.sqrt((a ** 2).sum(axis=1))
        a = a / weight[:, None]
    t = np.eye(k)
    f, gq = criterion(a @ t)
    grad = a.T @ gq
    alpha = 1.0
    for _ in range(maxit):
        m = t.T @ grad
        s_mat = (m + m.T) / 2
        gp = grad - t @ s_mat
        s = np.sqrt((gp ** 2).sum())
        if s < eps:
            break
        alpha *= 2
        for _ in range(11):
            x = t - alpha * gp
            u, _, vt = np.linalg.svd(x)
            candidate = u @ vt
            f_new, gq_new = criterion(a @ candidate)
            if f_new < f - 0.5 * s ** 2 * alpha:
                break
            alpha /= 2
        t, f, grad = candidate, f_new, a.T @ gq_new
    return (a @ t) * weight[:, None], t


def quartimax(a: np.ndarray, normalize: bool = True):
    def criterion(l):
        return -(l ** 4).sum() / 4, -(l ** 3)
    return _gpforth(a, criterion, normalize)


def equamax(a: np.ndarray, kappa: float, normalize: bool = True):
    """Crawford-Ferguson family with the given ``kappa`` (equamax: components / (2 * variables))."""
    p, k = a.shape
    n_mat = np.ones((k, k)) - np.eye(k)
    m_mat = np.ones((p, p)) - np.eye(p)

    def criterion(l):
        l2 = l ** 2
        f1 = (1 - kappa) * np.trace(l2.T @ (l2 @ n_mat)) / 4
        f2 = kappa * np.trace(l2.T @ (m_mat @ l2)) / 4
        gq = (1 - kappa) * l * (l2 @ n_mat) + kappa * l * (m_mat @ l2)
        return f1 + f2, gq
    return _gpforth(a, criterion, normalize)


def rotation_axes(ff: np.ndarray, vm: np.ndarray, eig: np.ndarray, weight, nfactors: int, rotations: list[str],
                  kind: str, quiet: bool = False, asked: list[str] | None = None, index=None) -> dict:
    """The plane of the first ``nfactors`` rotated components, one row per case, as ``{plane name: DataFrame}``.

    Port of ``rotationAxes``: the components are brought to unit variance, rotated, ordered by the variance
    they account for and signed so that their loadings add up to a positive number. With two components
    every rotation turns the plane within itself; from three on there is one plane for each pair
    (``"varimax 1-2"``...).
    """
    from .surcoin import _corr

    asked = asked or []
    if not rotations:
        return {}
    if kind != "pca":
        if not quiet:
            warnings.warn("rotations are only made of principal components: " + ", ".join(rotations)
                          + f' ignored with type="{kind}"')
        return {}
    k = min(nfactors, ff.shape[1])
    if k < 2:
        if not quiet:
            warnings.warn("a rotation needs at least two components: " + ", ".join(rotations)
                          + f" ignored with nfactors={nfactors}")
        return {}
    if quiet and k == 2:
        rotations = [r for r in rotations if r in ["varimax"] + list(asked)]
    z = ff[:, :k] / np.sqrt(eig[:k])
    loadings = _corr(vm, z, weight)
    out: dict = {}
    for rot in rotations:
        try:
            if rot == "varimax":
                load, theta = varimax(loadings, normalize=True)
            elif rot == "quartimax":
                load, theta = quartimax(loadings, normalize=True)
            else:
                load, theta = equamax(loadings, k / (2 * loadings.shape[0]), normalize=True)
        except Exception as exc:                                              # pragma: no cover
            warnings.warn(f"the {rot} rotation could not be made: {exc}")
            continue
        scores = z @ theta
        order = np.argsort(-(load ** 2).sum(axis=0), kind="stable")
        scores, load = scores[:, order], load[:, order]
        flip = np.where(load.sum(axis=0) < 0, -1.0, 1.0)
        scores = scores * flip
        pct = (load ** 2).sum(axis=0) / vm.shape[1] * 100
        columns = [f"RC{i + 1} ({p:.1f}%)" for i, p in enumerate(pct)]
        frame = pd.DataFrame(scores, columns=columns)
        if index is not None:
            frame.index = index
        if k == 2:
            out[rot] = frame
        else:
            for i, j in combinations(range(k), 2):
                out[f"{rot} {i + 1}-{j + 1}"] = frame.iloc[:, [i, j]]
    return out
