"""Latent models for the groups of a scattergram: profiles (LPA) and classes (LCA).

``run_latent`` estimates, for several numbers of groups, a latent profile analysis (quantitative variables:
a mixture of normals with the same diagonal covariance in every group, the ``EEI`` of mclust) or a latent
class analysis (qualitative variables: independent categorical variables inside each class, as poLCA).
Both are fitted by the EM algorithm from several starts, keeping the best.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

LCA_NREP = 5


def latent_model(kind: str) -> str:
    """``"LPA"`` for a principal components analysis, ``"LCA"`` for a correspondence analysis."""
    return "LPA" if kind == "pca" else "LCA"


# --------------------------------------------------------------------------- LPA ----
def _logsumexp(a: np.ndarray, axis: int = 1) -> np.ndarray:
    m = a.max(axis=axis, keepdims=True)
    return (m + np.log(np.exp(a - m).sum(axis=axis, keepdims=True))).ravel()


def _lpa_once(x: np.ndarray, k: int, rng: np.random.Generator, iterations: int = 500, tol: float = 1e-10):
    n, d = x.shape
    from sklearn.cluster import KMeans

    labels = KMeans(n_clusters=k, n_init=1, random_state=int(rng.integers(0, 2 ** 31 - 1))).fit(x).labels_ if k > 1 \
        else np.zeros(n, int)
    resp = np.eye(k)[labels]
    ll_old = -np.inf
    for _ in range(iterations):
        weights = resp.sum(axis=0)
        pi = weights / n
        means = (resp.T @ x) / np.maximum(weights, 1e-12)[:, None]
        var = np.maximum(sum((resp[:, [j]] * (x - means[j]) ** 2).sum(axis=0) for j in range(k)) / n, 1e-8)
        log_prob = np.column_stack([
            np.log(pi[j] + 1e-300) - 0.5 * (np.log(2 * np.pi * var) + (x - means[j]) ** 2 / var).sum(axis=1)
            for j in range(k)])
        total = _logsumexp(log_prob)
        ll = float(total.sum())
        resp = np.exp(log_prob - total[:, None])
        if abs(ll - ll_old) < tol * (abs(ll) + 1):
            break
        ll_old = ll
    return ll, resp.argmax(axis=1), k - 1 + k * d + d


def fit_lpa(x: np.ndarray, k: int, rng: np.random.Generator, starts: int = 10):
    """Latent profile model with equal diagonal variances: ``(log-likelihood, classes, BIC)`` (lower BIC is better)."""
    best = max((_lpa_once(x, k, rng) for _ in range(starts if k > 1 else 1)), key=lambda r: r[0])
    ll, classes, npar = best
    return ll, classes, -2 * ll + npar * np.log(len(x))


# --------------------------------------------------------------------------- LCA ----
def _lca_once(codes: np.ndarray, levels: list[int], k: int, rng: np.random.Generator, iterations: int = 1000,
              tol: float = 1e-10):
    n, d = codes.shape
    pi = np.full(k, 1.0 / k)
    probs = [rng.dirichlet(np.ones(m), size=k) for m in levels]            # per variable: k x categories
    ll_old = -np.inf
    for _ in range(iterations):
        log_resp = np.tile(np.log(pi), (n, 1))
        for v in range(d):
            log_resp += np.log(probs[v][:, codes[:, v]].T + 1e-300)
        total = _logsumexp(log_resp)
        ll = float(total.sum())
        resp = np.exp(log_resp - total[:, None])
        pi = resp.sum(axis=0) / n
        for v in range(d):
            counts = np.zeros((k, levels[v]))
            for c in range(levels[v]):
                counts[:, c] = resp[codes[:, v] == c].sum(axis=0)
            probs[v] = counts / np.maximum(counts.sum(axis=1, keepdims=True), 1e-300)
        if abs(ll - ll_old) < tol * (abs(ll) + 1):
            break
        ll_old = ll
    npar = (k - 1) + k * sum(m - 1 for m in levels)
    return ll, resp.argmax(axis=1), npar


def fit_lca(data: pd.DataFrame, k: int, rng: np.random.Generator, nrep: int = LCA_NREP):
    """Latent class model: ``(log-likelihood, classes, BIC)``; refused when it has more parameters than patterns."""
    codes = np.column_stack([pd.Categorical(data[c].astype(str)).codes for c in data.columns]).astype(int)
    levels = [int(codes[:, j].max()) + 1 for j in range(codes.shape[1])]
    best = max((_lca_once(codes, levels, k, rng) for _ in range(nrep)), key=lambda r: r[0])
    ll, classes, npar = best
    if k > 1 and npar >= int(np.prod(levels)):
        raise ValueError("there are more parameters than response patterns")
    return ll, classes, -2 * ll + npar * np.log(len(codes))


def run_latent(data, model: str, ks, seed: int | None = None) -> dict:
    """Estimate the latent ``model`` ("LPA" or "LCA") for every number of groups in ``ks``.

    Returns ``{"clusters": {k: classes (1-based)}, "bic": {k: BIC}}``; a number of groups that cannot be
    estimated is reported with a warning and left out (its BIC is ``NaN``).
    """
    clusters, bic = {}, {}
    for k in ks:
        rng = np.random.default_rng(seed)
        try:
            if model == "LPA":
                x = np.asarray(data.to_numpy() if isinstance(data, pd.DataFrame) else data, dtype=float)
                _, classes, value = fit_lpa(x, int(k), rng)
            else:
                _, classes, value = fit_lca(pd.DataFrame(data), int(k), rng)
        except Exception as exc:
            warnings.warn(f"{model} could not be computed for {k} groups: {exc}")
            bic[int(k)] = np.nan
            continue
        clusters[int(k)] = np.asarray(classes) + 1
        bic[int(k)] = float(value)
    return {"clusters": clusters, "bic": bic}
