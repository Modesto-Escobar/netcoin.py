"""Clusters and dimensional reductions of a table of cases: ``looking4clusters``."""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform

from .dichotomize import _as_character, _collation_key

DISTANCES = ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
AGGLOMERATIONS = ("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")


# ---------------------------------------------------------------- helpers ----
def as_factor(values) -> pd.Categorical:
    """R's ``as.factor``: numbers give their sorted values as levels, text its sorted words."""
    series = pd.Series(values)
    if isinstance(series.dtype, pd.CategoricalDtype):
        return pd.Categorical(series.astype(object), categories=list(series.cat.categories))
    if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
        levels = sorted(series.dropna().unique())
        text = dict(zip(levels, _as_character(pd.Series(levels))))
        return pd.Categorical(series.map(text), categories=[text[v] for v in levels])
    text = pd.Series(_as_character(series), index=series.index).where(series.notna())
    return pd.Categorical(text, categories=sorted(set(text.dropna()), key=_collation_key))


def clean_names(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "_", name)


def dissimilarity(data: np.ndarray, distance: str = "euclidean") -> np.ndarray:
    """Condensed dissimilarities between the rows, as R's ``dist`` (unknown names mean euclidean)."""
    if distance not in DISTANCES:
        distance = "euclidean"
    if distance == "canberra":
        return _canberra(data)
    metric = {"euclidean": "euclidean", "maximum": "chebyshev", "manhattan": "cityblock",
              "binary": "jaccard", "minkowski": "euclidean"}[distance]
    x = data != 0 if distance == "binary" else data
    return pdist(x, metric=metric)


def _canberra(data: np.ndarray) -> np.ndarray:
    """R's Canberra distance: the terms 0/0 are left out and the sum is scaled to all the variables."""
    n, m = data.shape
    out = []
    for i in range(n - 1):
        a, b = data[i], data[i + 1:]
        num, den = np.abs(a - b), np.abs(a + b)
        with np.errstate(divide="ignore", invalid="ignore"):
            terms = np.where((num == 0) & (den == 0), np.nan, num / den)
        count = (~np.isnan(terms)).sum(axis=1)
        total = np.nansum(terms, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            out.append(np.where(count > 0, total * m / count, np.nan))
    return np.concatenate(out)


def calinski_harabasz(condensed: np.ndarray, labels: np.ndarray) -> float:
    """The Calinski-Harabasz index of a clustering from dissimilarities (as ``fpc::cluster.stats``)."""
    dmat = squareform(condensed)
    n = len(labels)
    clusters = sorted(set(labels))
    k = max(clusters)
    overall = (condensed ** 2).sum() / n
    within = 0.0
    for c in clusters:
        idx = np.where(labels == c)[0]
        if len(idx) > 1:
            within += (squareform(dmat[np.ix_(idx, idx)], checks=False) ** 2).sum() / len(idx)
    return float((overall - within) * (n - k) / (within * (k - 1)))


def cutree(tree: np.ndarray, k: int) -> np.ndarray:
    """The ``k`` groups after the first ``n - k`` merges of a linkage matrix (R's ``cutree``).

    Unlike a cut at a height it follows the order of the merges, so it also works when the tree
    has inversions (``centroid`` and ``median``). Groups are numbered by first appearance.
    """
    n = len(tree) + 1
    parent = list(range(2 * n - 1))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for row, (a, b, _, _) in enumerate(tree[:n - k]):
        parent[find(int(a))] = parent[find(int(b))] = n + row
    return _first_appearance([find(i) for i in range(n)])


def _first_appearance(labels) -> np.ndarray:
    """Number the groups in order of appearance of the cases (as ``cutree`` does)."""
    seen: dict = {}
    return np.array([seen.setdefault(l, len(seen) + 1) for l in labels])


# -------------------------------------------------------------------- pam ----
def pam(dmat: np.ndarray, k: int) -> tuple[list[int], np.ndarray]:
    """Partitioning around medoids (Kaufman and Rousseeuw: BUILD and then SWAP).

    ``dmat`` is the square matrix of dissimilarities. Returns the (sorted) positions of the
    medoids and the cluster (1..k, numbered by medoid) of every case.
    """
    n = len(dmat)
    medoids: list[int] = [int(np.argmin(dmat.sum(axis=1)))]
    nearest = dmat[medoids[0]].copy()
    while len(medoids) < k:
        chosen = np.zeros(n, bool)
        chosen[medoids] = True
        gains = np.maximum(nearest[None, :] - dmat, 0)[:, ~chosen].sum(axis=1)
        gains[chosen] = -np.inf                                  # candidates are the cases not chosen yet
        pick = int(np.argmax(gains))
        medoids.append(pick)
        nearest = np.minimum(nearest, dmat[pick])
    while True:
        sub = dmat[medoids]                                    # k x n
        order = np.argsort(sub, axis=0, kind="stable")
        first = sub[order[0], np.arange(n)]
        second = sub[order[1], np.arange(n)] if k > 1 else np.full(n, np.inf)
        which = order[0]
        best, move = -1e-10, None
        for h in range(n):
            if h in medoids:
                continue
            dh = dmat[h]
            for slot in range(k):
                owned = which == slot
                change = np.where(owned, np.minimum(dh, second) - first, np.minimum(dh - first, 0.0))
                total = change.sum()
                if total < best:
                    best, move = total, (slot, h)
        if move is None:
            break
        medoids[move[0]] = move[1]
    medoids = sorted(medoids)
    labels = _first_appearance(np.argmin(dmat[medoids], axis=0))
    return medoids, labels


# -------------------------------------------------------------- the object ----
@dataclass
class Looking4Clusters:
    """Result of :func:`looking4clusters`: the data, its clusterings and its reductions.

    ``clusters`` maps the name of a method (``kmeans``, ``pam``, ``hclust``...) to a categorical
    Series, or to a DataFrame with one column per number of groups (``levels_2``, ``levels_3``...);
    ``optim`` holds the number of groups the method found best. ``reductions`` maps a name
    (``pca``, ``tsne``, ``cmds``...) to a DataFrame of two coordinates, and ``pctvar`` the
    percentage of the variance of those that have it.
    """

    data: np.ndarray
    samples: list[str]
    variables: list[str]
    options: dict = field(default_factory=dict)
    clusters: dict = field(default_factory=dict)
    reductions: dict = field(default_factory=dict)
    optim: dict = field(default_factory=dict)
    pctvar: dict = field(default_factory=dict)

    def add_cluster(self, data, name: str | None = None, group_stats_by: bool = False,
                    my_groups: bool = False, optim_cluster: bool = False) -> "Looking4Clusters":
        """Add a clustering (one group per case); the same name again adds a column."""
        if len(data) != len(self.samples):
            raise ValueError("data: there must be one per sample")
        name = f"cluster_{len(self.clusters) + 1}" if name is None else clean_names(name)
        factor = pd.Series(as_factor(data), index=self.samples)
        if name in self.clusters:
            current = self.clusters[name]
            if isinstance(current, pd.Series):
                frame = pd.DataFrame({f"levels_{len(current.cat.categories)}": current,
                                      f"levels_{len(factor.cat.categories)}": factor})
                self.clusters[name] = frame
                self.optim[name] = len(current.cat.categories)
            else:
                current[f"levels_{len(factor.cat.categories)}"] = factor
            if optim_cluster:
                self.optim[name] = len(factor.cat.categories)
        else:
            self.clusters[name] = factor
        if group_stats_by:
            self.options["groupStatsBy"] = self.options.get("groupStatsBy", []) + [name]
        if my_groups:
            self.options["myGroups"] = name
        return self

    def add_reduction(self, data, name: str | None = None) -> "Looking4Clusters":
        """Add the first two columns of ``data`` as a dimensional reduction."""
        frame = pd.DataFrame(data)
        if len(frame) != len(self.samples):
            raise ValueError("data: there must be one row per sample")
        name = f"reduction_{len(self.reductions) + 1}" if name is None else clean_names(name)
        out = frame.iloc[:, :2].copy()
        out.index = self.samples
        self.reductions[name] = out
        return self

    def __repr__(self) -> str:
        lines = ["An object of class looking4clusters",
                 f"{len(self.variables)} variables across {len(self.samples)} samples"]
        if self.clusters:
            lines.append(f"{len(self.clusters)} clusters added: {', '.join(self.clusters)}")
        if self.reductions:
            lines.append(f"{len(self.reductions)} dimensional reductions added: {', '.join(self.reductions)}")
        return "\n".join(lines)


# ---------------------------------------------------------------- the runs ----
def _k_range(selectedk, n: int) -> tuple[int, int]:
    if selectedk is None:
        kmin, kmax = 2, 10
    else:
        kmin, kmax = selectedk - 5, selectedk + 5
    if kmin <= 1:
        kmin = 2
    if kmax >= n:
        kmax = n - 1
    if kmin >= kmax:
        kmin = kmax - 1
    return kmin, kmax


def run_kmeans(obj: Looking4Clusters, selectedk=None, seed=None) -> Looking4Clusters:
    """k-means for a range of numbers of groups; the best by Calinski-Harabasz (needs scikit-learn)."""
    try:
        from sklearn.cluster import KMeans
    except ImportError:                                                   # pragma: no cover
        warnings.warn("Install 'scikit-learn' to get kmeans clusters.")
        return obj
    data = obj.data
    n = len(data)
    kmin, kmax = _k_range(selectedk, n)
    total = ((data - data.mean(axis=0)) ** 2).sum()
    fits = {}
    for k in range(kmin, kmax + 1):
        model = KMeans(n_clusters=k, n_init=1, max_iter=15, init="random", random_state=seed).fit(data)
        fits[k] = (model.labels_, float(model.inertia_))
    if selectedk is None:
        index = {k: ((total - w) * (n - k)) / (w * (k - 1)) for k, (_, w) in fits.items()}
        best = max(index, key=lambda k: (index[k], -k))
    else:
        best = selectedk
    for k, (labels, _) in fits.items():
        obj.add_cluster([f"kmeans_{l + 1}" for l in labels], "kmeans", optim_cluster=len(set(labels)) == best)
    return obj


def run_pam_hclust(obj: Looking4Clusters, distance="euclidean", agglomeration="complete",
                   selectedk=None) -> Looking4Clusters:
    """PAM and hierarchical clusters for a range of numbers of groups."""
    data = obj.data
    n = len(data)
    kmin, kmax = _k_range(selectedk, n)
    condensed = dissimilarity(data, distance)
    dmat = squareform(condensed)

    pams = {k: pam(dmat, k)[1] for k in range(kmin, kmax + 1)}
    if selectedk is None:
        scores = {k: calinski_harabasz(condensed, pams[k]) for k in pams}
        best_pam = max(scores, key=lambda k: (scores[k], -k))
    else:
        best_pam = selectedk

    method = agglomeration if agglomeration in AGGLOMERATIONS else "complete"
    scipy_method, y = {"ward.D2": ("ward", condensed), "ward.D": ("ward", np.sqrt(condensed)),
                       "single": ("single", condensed), "complete": ("complete", condensed),
                       "average": ("average", condensed), "mcquitty": ("weighted", condensed),
                       "median": ("median", np.sqrt(condensed)),
                       "centroid": ("centroid", np.sqrt(condensed))}[method]
    tree = linkage(y, method=scipy_method)
    cuts = {k: cutree(tree, k) for k in range(kmin, kmax + 1)}
    if selectedk is None:
        scores = {k: calinski_harabasz(condensed, cuts[k]) for k in cuts}
        best_hclust = max(scores, key=lambda k: (scores[k], -k))
    else:
        best_hclust = selectedk

    for k in range(kmin, kmax + 1):
        obj.add_cluster([f"pam_{l}" for l in pams[k]], "pam", optim_cluster=len(set(pams[k])) == best_pam)
        obj.add_cluster([f"hclust_{l}" for l in cuts[k]], "hclust", optim_cluster=len(set(cuts[k])) == best_hclust)
    return obj


def run_pca(obj: Looking4Clusters) -> Looking4Clusters:
    """The first two principal components (not scaled), with the percentage of variance of each."""
    centered = obj.data - obj.data.mean(axis=0)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    scores = centered @ vt.T[:, :2]
    obj.add_reduction(pd.DataFrame(scores, columns=["PC1", "PC2"]), "pca")
    variance = s ** 2
    obj.pctvar["pca"] = (variance[:2] / variance.sum() * 100).tolist()
    return obj


def run_tsne(obj: Looking4Clusters, perplex: float = 30, max_iter: int = 1000, seed=None) -> Looking4Clusters:
    """t-SNE (needs scikit-learn); the result changes with ``seed``."""
    try:
        from sklearn.manifold import TSNE
    except ImportError:                                                   # pragma: no cover
        warnings.warn("Install 'scikit-learn' to get tsne dimensionality reduction.")
        return obj
    n = len(obj.data)
    if 3 * perplex > n - 1:
        perplex = (n - 1) / 3
    points = TSNE(n_components=2, perplexity=perplex, max_iter=max_iter, random_state=seed,
                  init="pca").fit_transform(obj.data)
    return obj.add_reduction(pd.DataFrame(points, columns=["tSNE1", "tSNE2"]), "tsne")


def run_mds(obj: Looking4Clusters) -> Looking4Clusters:
    """Classical multidimensional scaling of the Canberra distances."""
    dmat = squareform(_canberra(obj.data))
    n = len(dmat)
    centering = np.eye(n) - 1.0 / n
    b = -0.5 * centering @ (dmat ** 2) @ centering
    values, vectors = np.linalg.eigh((b + b.T) / 2)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    points = vectors[:, :2] * np.sqrt(np.maximum(values[:2], 0))
    obj.add_reduction(pd.DataFrame(points, columns=["cMDS1", "cMDS2"]), "cmds")
    obj.pctvar["cmds"] = (values[:2] / np.maximum(values, 0).sum() * 100).tolist()
    return obj


def run_nmf(obj: Looking4Clusters, seed=None) -> Looking4Clusters:
    """Non-negative matrix factorisation with two components (needs scikit-learn)."""
    data = obj.data
    if (data.sum(axis=0) == 0).any():
        warnings.warn("Your data has columns that are all zeros, this is not supported for nmf.")
    elif (data.sum(axis=1) == 0).any():
        warnings.warn("Your data has rows that are all zeros, this is not supported for nmf.")
    elif (data < 0).any():
        warnings.warn("Your data contains some negative entries, this is not supported for nmf.")
    else:
        try:
            from sklearn.decomposition import NMF
        except ImportError:                                               # pragma: no cover
            warnings.warn("Install 'scikit-learn' to get nmf dimensionality reduction.")
            return obj
        basis = NMF(n_components=2, init="nndsvda", random_state=seed, max_iter=1000).fit_transform(data)
        obj.add_reduction(pd.DataFrame(basis, columns=["NMF1", "NMF2"]), "nmf")
    return obj


def run_umap(obj: Looking4Clusters, seed=None) -> Looking4Clusters:
    """UMAP (needs ``umap-learn``)."""
    try:
        import umap
    except ImportError:
        warnings.warn("Install 'umap-learn' to get umap dimensionality reduction.")
        return obj
    n = len(obj.data)
    neighbors = 15
    if neighbors > n / 3:
        neighbors = int(np.floor(n / 3))
    neighbors = max(neighbors, 2)
    points = umap.UMAP(n_neighbors=neighbors, n_components=2, random_state=seed).fit_transform(obj.data)
    return obj.add_reduction(pd.DataFrame(points, columns=["UMAP1", "UMAP2"]), "umap")


def n_components(variability: Sequence[float]) -> int:
    """How many principal components to keep, from the proportion of variance of each (R's ``Ncomponents``)."""
    keep, divisor = 4, 10
    variability = list(variability)
    if len(variability) * 0.9 > 1:
        variability = variability[:int(len(variability) * 0.9)]
    change = [a - b for a, b in zip(variability[:-1], variability[1:])]
    prop = np.mean(change) / np.mean(sorted(change, reverse=True)[:3]) * 100
    if 10 < prop <= 15:
        divisor = 5
    if prop > 15:
        return keep
    min_var = max(change) / divisor
    cutoff = [c < min_var for c in change]
    position, k = [], 0
    for i in range(len(variability) - 1):
        k = k + 1 if cutoff[i] else 0
        position.append(k)
        if k > 7 or (i + 1 > 10 and k > 3):
            ones = [j + 1 for j, p in enumerate(position) if p == 1]
            n = max(ones) - 1
            if n > 1:
                keep = n
    return keep


def _create(data, components: bool) -> Looking4Clusters:
    frame = pd.DataFrame(data)
    samples = [str(i) for i in frame.index] if not isinstance(frame.index, pd.RangeIndex) \
        else [f"sample_{i + 1}" for i in range(len(frame))]
    variables = [str(c) for c in frame.columns]
    matrix = frame.to_numpy(dtype=float)
    if components:
        if len(matrix) < 5:
            raise ValueError("You cannot apply to components for less than five samples.")
        if matrix.shape[1] < 5:
            raise ValueError("You cannot apply to components for less than five variables.")
        logged = np.log2(matrix + 1)
        centered = logged - logged.mean(axis=0)
        variance = np.linalg.svd(centered, compute_uv=False) ** 2
        nc = n_components(np.round(variance / variance.sum(), 5))      # summary.prcomp rounds to 5 digits
        raw = matrix - matrix.mean(axis=0)
        u, s, vt = np.linalg.svd(raw, full_matrices=False)
        matrix = raw @ vt.T[:, :nc]
        variables = [f"PC{i + 1}" for i in range(nc)]
    return Looking4Clusters(matrix, samples, variables)


def looking4clusters(
    data,
    groups=None,
    components: bool = False,
    running_all: bool = True,
    distance: str = "euclidean",
    agglomeration: str = "complete",
    selectedk: int | None = None,
    perplex: float = 30,
    max_iter: int = 1000,
    force_execution: bool = False,
    seed: int | None = None,
) -> Looking4Clusters:
    """Look for groups of cases and for two-dimensional views of them.

    Port of ``looking4clusters``. With ``running_all`` it runs k-means, PAM and hierarchical
    clustering for 2 to 10 groups (or ``selectedk`` plus or minus five) and the reductions PCA,
    t-SNE, classical MDS (on Canberra distances), NMF and UMAP; each needs its own library and
    is skipped with a warning if it is missing. The best number of groups of every method is
    stored in ``optim`` (Calinski-Harabasz).

    Parameters
    ----------
    data : DataFrame or 2-d array
        One row per case, one numeric column per variable.
    groups : sequence, optional
        Groups of your own, one per case.
    components : bool
        Work on the principal components of ``log2(data + 1)`` that explain the data instead of
        on the variables (at least five cases and five variables).
    running_all : bool
        Run the methods; otherwise only your ``groups`` are stored.
    distance : str
        ``"euclidean"``, ``"maximum"``, ``"manhattan"``, ``"canberra"``, ``"binary"`` or
        ``"minkowski"`` for PAM and hierarchical clustering.
    agglomeration : str
        ``"ward.D"``, ``"ward.D2"``, ``"single"``, ``"complete"``, ``"average"``,
        ``"mcquitty"``, ``"median"`` or ``"centroid"``.
    selectedk : int, optional
        The number of groups you prefer.
    perplex, max_iter
        Of t-SNE.
    force_execution : bool
        Run the slow methods also with more than 5000 cases (or variables).
    seed : int, optional
        Seed of the methods that start at random (k-means, t-SNE, NMF, UMAP).

    Returns
    -------
    Looking4Clusters
    """
    if np.ndim(data) != 2:
        raise ValueError("data: incorrect dimensions. A matrix like object is required.")
    frame = pd.DataFrame(data)
    if running_all and (frame.shape[1] > 5000 or frame.shape[0] > 5000):
        warnings.warn("Too large matrix, could cause performance problems with some methods"
                      + ("" if force_execution else ", they will be omitted"))
    obj = _create(data, components)
    if groups is not None:
        obj.add_cluster(groups, my_groups=True)
    if running_all:
        run_kmeans(obj, selectedk, seed)
        if force_execution or not frame.shape[0] > 5000:
            run_pam_hclust(obj, distance, agglomeration, selectedk)
        run_pca(obj)
        run_tsne(obj, perplex, max_iter, seed)
        if force_execution or not frame.shape[0] > 5000:
            run_mds(obj)
        if force_execution or not frame.shape[1] > 5000:
            run_nmf(obj, seed)
        run_umap(obj, seed)
    return obj
