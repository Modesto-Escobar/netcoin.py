"""Scattergram of the cases of a survey: ``sur_scat`` (port of ``surScat``)."""
from __future__ import annotations

import logging
import warnings
from typing import Sequence

import numpy as np
import pandas as pd

from .core import _normalize_legacy, name_by_language, netCoin
from .dichotomize import _as_character, _collation_key, dichotomize
from .layouts import _fix_signs
from .surcoin import _corr
from .utils import casesList, getByLanguage, groupList, groupsList

log = logging.getLogger("netcoin")
_INDEXES = ("ch", "silhouette", "bic")


def _as_list(value) -> list:
    if value is None:
        return []
    return [value] if isinstance(value, str) or not hasattr(value, "__iter__") else list(value)


def _is_numeric(series: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series)


# ------------------------------------------------------------------- analyses ----
def mca_rows(matrix: np.ndarray, nfactors: int, weight=None) -> tuple[np.ndarray, np.ndarray]:
    """Principal coordinates of the cases of a simple correspondence analysis of an indicator matrix.

    Port of ``layoutMca(rows = TRUE)``. Returns the coordinates on factors 2 to ``nfactors + 1`` of
    the analysis (the first dimension is the trivial one) and the non-trivial eigenvalues.
    """
    m = np.asarray(matrix, dtype=float)
    p = m / m.shape[0] if weight is None else m * (np.asarray(weight, dtype=float) / np.sum(weight))[:, None]
    col_mass, row_mass = p.sum(axis=0), p.sum(axis=1)
    expected = np.outer(row_mass, col_mass)
    z = (p - expected) / np.sqrt(expected)
    u, d, _ = np.linalg.svd(z, full_matrices=False)
    standard = u / np.sqrt(row_mass)[:, None]
    principal = standard * d
    return _fix_signs(principal[:, 1:nfactors + 1]), d[1:] ** 2


def pca_scores(b: np.ndarray, weight=None, n_extract: int = 2) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scores on the principal components of the correlations, eigenvalues and the standardized data."""
    if weight is None:
        std = (b - b.mean(axis=0)) / b.std(axis=0, ddof=1)
        u, s, vt = np.linalg.svd(std, full_matrices=False)
        scores = _fix_signs(std @ vt.T[:, :n_extract])
        return scores, s ** 2 / (len(b) - 1), std
    w = np.asarray(weight, dtype=float) / np.sum(weight)
    centered = b - (b * w[:, None]).sum(axis=0)
    std = centered / np.sqrt((centered ** 2 * w[:, None]).sum(axis=0))
    values, vectors = np.linalg.eigh(_corr(b, b, weight))
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    return _fix_signs(std @ vectors[:, :n_extract]), values, std


# ---------------------------------------------------------------------- axes ----
def _axis_spec(spec, default: int, argname: str, data: pd.DataFrame) -> dict:
    out = {"variable": None, "factor": default, "scale": None}
    if spec is None:
        return out
    items = _as_list(spec) if not isinstance(spec, tuple) else list(spec)
    if len(items) < 1 or len(items) > 2:
        raise ValueError(f"{argname} must be an axis, optionally followed by True or False")
    if len(items) == 2:
        if not isinstance(items[1], (bool, np.bool_)):
            raise ValueError(f"the second element of {argname} must be True or False")
        out["scale"] = bool(items[1])
    first = str(items[0])
    if first in data.columns:                                  # a column prevails over the "F1" reading
        out["variable"] = first
        return out
    text = first.upper()
    for prefix in ("DIM", "PC", "F"):
        if text.startswith(prefix):
            text = text[len(prefix):]
            break
    try:
        n = int(text)
    except ValueError:
        n = 0
    if n < 1:
        raise ValueError(f"{argname} must be a factorial axis (1, \"F1\", \"PC2\"...) or a variable of the data frame")
    out["factor"] = n
    return out


# --------------------------------------------------------------- group helpers ----
def pattern_groups(values, groups: list[str]) -> pd.Categorical:
    """Ordered factor of the group (or combination of groups) of every pattern; groups lead the levels."""
    values = [str(v) for v in values]
    extra = sorted(set(values) - set(groups), key=_collation_key)
    return pd.Categorical(values, categories=list(groups) + extra, ordered=True)


def cluster_shares(node: np.ndarray, n_nodes: int, group: pd.Categorical, w=None) -> pd.DataFrame:
    """The (weighted) cases of every group in every node: one row per node, one column per group."""
    w = np.ones(len(node)) if w is None else np.asarray(w, dtype=float)
    if isinstance(group, pd.Series):
        group = group.cat.as_unordered() if hasattr(group, "cat") else pd.Categorical(group)
        group = pd.Categorical(group.astype(object), categories=list(group.cat.categories))
    present = [g for g in group.categories if (group == g).any()]
    out = np.zeros((n_nodes, len(present)))
    codes = {g: j for j, g in enumerate(present)}
    for i, (n, g) in enumerate(zip(node, group)):
        out[int(n) - 1, codes[g]] += w[i]
    return pd.DataFrame(out, columns=present)


def _group_weights(net, column: str, w=None) -> pd.DataFrame:
    n = len(net.nodes)
    shares = net.scat.get("cluster_shares", {}).get(column)
    if shares is not None and len(shares) == n:
        return shares
    group = net.nodes[column].astype(str)
    group = group.map(lambda s: s.split("|")[0])
    levels = sorted(set(group), key=_collation_key)
    weights = np.ones(n) if w is None else np.asarray(w, dtype=float)
    return pd.DataFrame({lv: np.where(group == lv, weights, 0.0) for lv in levels})


def _group_means(x: pd.DataFrame, w: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(np.nan, index=w.columns, columns=x.columns)
    for col in x.columns:
        ok = x[col].notna().to_numpy()
        weights = w.to_numpy()[ok]
        total = weights.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[col] = np.where(total > 0, (weights * x[col].to_numpy()[ok][:, None]).sum(axis=0) / total, np.nan)
    return out


def current_layouts(net) -> dict:
    """The planes of a scattergram: the one in ``fx``/``fy`` (named after the analysis) and any added."""
    if net.layouts:
        return net.layouts
    nodes = net.nodes
    if not {"fx", "fy"} <= set(nodes.columns):
        raise ValueError("the object holds no coordinates")
    labels = net.options.get("axesLabels")
    first = pd.DataFrame({"x": nodes["fx"].to_numpy(), "y": nodes["fy"].to_numpy()})
    if labels is not None and len(labels) >= 2:
        first.columns = [str(l) for l in list(labels)[:2]]
        import re
        key = "-".join(re.sub(r" *\([^()]*\) *$", "", c) for c in first.columns)
        if net.scat.get("plane_name"):
            key = net.scat["plane_name"]
    else:
        first.columns, key = ["fx", "fy"], "layout1"
    return {key: first}


def current_clusters(net) -> dict:
    """The centroid of every group of every clusterization on every plane (``name.plane``)."""
    columns = [c for c in net.scat.get("cluster_columns", []) if c in net.nodes.columns]
    out = {}
    if not columns:
        return out
    layouts = current_layouts(net)
    idx = net.scat.get("case_to_pattern")
    weights = None
    if idx is not None:
        case_w = net.scat.get("case_weight")
        cw = np.ones(len(idx)) if case_w is None or len(case_w) != len(idx) else np.asarray(case_w, dtype=float)
        weights = np.bincount(np.asarray(idx) - 1, weights=cw, minlength=len(net.nodes))
    for column in columns:
        w = _group_weights(net, column, weights)
        for plane, frame in layouts.items():
            out[f"{column}.{plane}"] = _group_means(frame.reset_index(drop=True), w)
    return out


def current_cluster_means(net) -> dict | None:
    """The mean of every group on every numeric column of the node table (``fx`` and ``fy`` left out)."""
    columns = [c for c in net.scat.get("cluster_columns", []) if c in net.nodes.columns]
    nodes = net.nodes
    variables = [c for c in nodes.columns if _is_numeric(nodes[c]) and c not in ("fx", "fy")]
    if not columns or not variables:
        return None
    idx = net.scat.get("case_to_pattern")
    weights = None
    if idx is not None:
        case_w = net.scat.get("case_weight")
        cw = np.ones(len(idx)) if case_w is None or len(case_w) != len(idx) else np.asarray(case_w, dtype=float)
        weights = np.bincount(np.asarray(idx) - 1, weights=cw, minlength=len(nodes))
    return {c: _group_means(nodes[variables].astype(float), _group_weights(net, c, weights)) for c in columns}


# -------------------------------------------------------------------- k-means ----
def hartigan_refine(x: np.ndarray, labels: np.ndarray, k: int, max_passes: int = 100) -> np.ndarray:
    """Move cases between groups while that lowers the within-group sum of squares.

    A case in a group of ``n1`` cases moves to another of ``n2`` cases if
    ``n2 / (n2 + 1) * d(x, c2)^2 < n1 / (n1 - 1) * d(x, c1)^2``: the optimal-transfer step of the
    Hartigan-Wong algorithm that R's ``kmeans`` uses, which leaves k-means at a better optimum
    than Lloyd's alone.
    """
    labels = labels.copy()
    centres = np.array([x[labels == g].mean(axis=0) for g in range(k)])
    sizes = np.bincount(labels, minlength=k).astype(float)
    for _ in range(max_passes):
        moved = False
        for i in range(len(x)):
            g = labels[i]
            if sizes[g] <= 1:
                continue
            d2 = ((x[i] - centres) ** 2).sum(axis=1)
            cost = sizes / (sizes + 1) * d2
            cost[g] = sizes[g] / (sizes[g] - 1) * d2[g]
            h = int(np.argmin(cost))
            if h != g and cost[h] < cost[g] - 1e-12:
                centres[g] = (centres[g] * sizes[g] - x[i]) / (sizes[g] - 1)
                centres[h] = (centres[h] * sizes[h] + x[i]) / (sizes[h] + 1)
                sizes[g] -= 1
                sizes[h] += 1
                labels[i] = h
                moved = True
        if not moved:
            break
    return labels


# -------------------------------------------------------------------- sur_scat ----
def sur_scat(
    data: pd.DataFrame,
    variables: Sequence[str] | None = None,
    active: Sequence[str] | None = None,
    columns: Sequence[str] | None = None,
    patterns: bool = False,
    v_patterns: Sequence[str] | None = None,
    weight: str | Sequence[float] | None = None,
    seed: int | None = 2020,
    max_n: int = 2000,
    jitter: float | bool = 0,
    type: str | None = None,
    xaxis=None,
    yaxis=None,
    scale_axes: bool = True,
    nfactors: int = 2,
    crit_factors: float = 0,
    layouts=None,
    suffix: str | None = None,
    align_to=None,
    nclusters=None,
    cluster_on: str = "factors",
    scale_clusters: bool | None = None,
    nstart: int = 25,
    cluster_index: str = "ch",
    sort_clusters: bool = True,
    clusters=None,
    force_execution: bool = False,
    **netcoin_args,
):
    """A scattergram of the cases of a survey, coloured by groups (port of ``surScat``).

    The cases (or the patterns of answers, with ``patterns``) are drawn on the plane of a
    multiple correspondence analysis (qualitative active variables) or of a principal components
    analysis (quantitative ones), or on observed variables, and grouped with k-means. The result
    is a :class:`netCoin` object without links, whose nodes hold the variables, the groups and
    the coordinates; its attribute ``scat`` keeps what the later steps use (the case-to-pattern
    map, the shares of every group in every pattern...), and ``clusters`` and ``cluster_means`` the
    centroids of the groups.

    Parameters
    ----------
    data : DataFrame
        One row per case.
    variables, active : list of str, optional
        The variables shown in the nodes and those the analysis is made of (``active`` defaults to
        ``variables``). A variable used on an axis joins ``variables``.
    columns : list of str, optional
        Columns of the node table to put first (``"_vpattern"``, ``"_actives"``, ``"_others"`` and
        ``"_nclusters"`` stand for those blocks).
    patterns, v_patterns : bool, list of str
        Merge the cases that give the same answers to ``v_patterns`` into one node (a node of
        patterns is sized by its number of cases). Stating ``v_patterns`` implies ``patterns``.
    weight : str or sequence, optional
        Weight of each case.
    seed : int, optional
        Seed of the random steps (k-means, samples, jitter).
    max_n : int
        At most this many nodes are drawn, a random sample of them if there are more.
    jitter : float or bool
        A little random spread of the coordinates, as a fraction of their range (``True`` is 0.02).
    type : {"mca", "pca"}, optional
        Kind of analysis; by default ``"pca"`` if all the active variables are numbers.
    xaxis, yaxis : int, str or (axis, bool), optional
        What is drawn on each axis: a factorial axis (``1``, ``"F1"``, ``"PC2"``) or a numeric (or
        ordered) variable, optionally followed by whether to standardize it. With a variable on both
        axes there are no groups by default.
    scale_axes : bool
        Standardize the factorial axes.
    nfactors, crit_factors : int, float
        Factors used to cluster: those of the first ``nfactors`` whose eigenvalue exceeds
        ``crit_factors`` times the mean eigenvalue.
    nclusters : int, "auto" or list
        Number(s) of k-means groups, one column each (the last stated colours the graph). 0 draws
        no groups. ``"auto"`` chooses among the numbers stated next to it (2 to 10) with
        ``cluster_index``.
    cluster_on : {"factors", "variables"}
        What the k-means works on.
    scale_clusters : bool, optional
        Standardize the clustering space (by default yes, except indicators of categories).
    nstart : int
        Random starts of the k-means (the best is kept).
    cluster_index : {"ch", "silhouette", "bic"}
        Index that chooses the number of groups: Calinski-Harabasz, silhouette, or BIC of a latent
        model (LPA for ``type="pca"``, LCA for ``"mca"``, lower is better); with ``"auto"`` next to a 1 the BIC
        also weighs the absence of groups.
    sort_clusters : bool
        Number the groups from left to right on the horizontal axis.
    layouts : bool or list of str, optional
        Planes added to the main one: the reductions of ``looking4clusters`` (``"pca"``, ``"tsne"``,
        ``"cmds"``, ``"nmf"``, ``"umap"``) and the rotations of the components (``"varimax"``, ``"quartimax"``,
        ``"equamax"``, only with ``type="pca"``); ``True`` or ``"all"`` asks for all of them.
    clusters : bool or list of str, optional
        Groups added to the main ones, by number of groups in ``nclusters``: ``"kmeans"``, ``"pam"``, ``"hclust"``
        (of ``looking4clusters``, on the active variables as they come) and the latent models ``"lpa"`` /
        ``"lca"``; ``True`` asks for all of them (the latent model that suits ``type``).
    suffix, align_to : str, str or int, optional
        With ``suffix``, every plane added is also shown aligned by Procrustes to the main plane (or to
        ``align_to``, a plane's name or position), under its name followed by the suffix.
    force_execution : bool
        Compute ``pam``, ``hclust`` and ``cmds`` even beyond 5000 cases.
    **netcoin_args
        The rest goes to :class:`netCoin` (``language``, ``main``, ``color``...).

    Returns
    -------
    netCoin
    """
    args = _normalize_legacy(netcoin_args, strict=False)
    stated_v_patterns = v_patterns is not None
    if stated_v_patterns:
        patterns = True
    if cluster_on not in ("factors", "variables"):
        raise ValueError("cluster_on must be 'factors' or 'variables'")
    cluster_index = str(cluster_index).lower()
    if cluster_index not in _INDEXES:
        raise ValueError("cluster_index must be one of: " + ", ".join(_INDEXES))
    stated_labels = "axes_labels" in args
    data = data.reset_index(drop=True) if isinstance(data.index, pd.RangeIndex) else data.copy()
    variables = list(data.columns) if variables is None else _as_list(variables)
    active = list(variables) if active is None else _as_list(active)
    v_patterns = list(variables) if v_patterns is None else _as_list(v_patterns)

    # ---- the number of groups
    def is_var_axis(spec):
        items = _as_list(spec) if not isinstance(spec, tuple) else list(spec)
        return spec is not None and len(items) > 0 and str(items[0]) in data.columns

    if nclusters is None:
        nclusters = 0 if (is_var_axis(xaxis) and is_var_axis(yaxis)) else 2
    given = _as_list(nclusters)
    auto_k = "auto" in given
    numbers = [x for x in given if x != "auto"]
    try:
        ks = [float(x) for x in numbers]
    except (TypeError, ValueError):
        raise ValueError("nclusters must hold positive whole numbers, optionally next to \"auto\", or be 0")
    no_k = len(ks) > 0 and all(k == 0 for k in ks)
    if any(k == 0 for k in ks) and (not no_k or len(ks) > 1 or auto_k):
        raise ValueError("nclusters=0 asks for no groups at all, so it must stand alone")
    if no_k:
        ks = []
    elif any(k < 1 or k != round(k) for k in ks):
        raise ValueError("nclusters must hold positive whole numbers, optionally next to \"auto\", or be 0")
    ks = [int(k) for k in ks]
    keep_all_k = not auto_k or len(ks) > 0
    stated_ks = len(ks) > 0
    weigh_one = False
    if 1 in ks:
        ks = [k for k in ks if k != 1]
        weigh_one = auto_k and cluster_index == "bic"        # the BIC can weigh no groups at all
        if not auto_k and not ks:
            no_k = True
            warnings.warn("nclusters=1 makes no groups, as a single group is no clustering: the plane is drawn "
                          "alone, as with nclusters=0")
        elif not weigh_one:
            warnings.warn("nclusters: 1 group left out, as a single group is no clustering")
    if auto_k and not stated_ks:
        ks = list(range(2, 11))
    if auto_k and not any(k >= 2 for k in ks):
        raise ValueError("nclusters=\"auto\" needs at least one number of groups above 1")
    if not auto_k and not no_k and not ks:
        raise ValueError("nclusters must state at least one number of groups")
    ks = list(dict.fromkeys(reversed(ks)))[::-1]                 # a number stated twice keeps its last place
    # planes and groups to be added by looking4clusters, checked before any computation
    from .scatl4c import add_l4c, l4c_clusters, l4c_layout, l4c_run, latent_clusters, rotation_layout
    l4c_lay, rot_lay, l4c_meth = l4c_layout(layouts), rotation_layout(layouts), l4c_clusters(clusters)
    if suffix is not None and (not isinstance(suffix, str) or not suffix):
        raise ValueError("suffix must be a single non-empty string")
    if suffix is not None and not l4c_lay and not rot_lay:
        warnings.warn("suffix aligns the planes added by layouts, and none was asked for")
    if align_to is not None:
        if isinstance(align_to, bool) or not isinstance(align_to, (str, int, np.integer)):
            raise ValueError("align_to must be None, or the name or the index of one plane")
        if suffix is None:
            warnings.warn("align_to sets the plane suffix aligns to, and no suffix was stated")
    if nfactors < 1:
        raise ValueError("nfactors must be at least 1")
    if nstart < 1:
        raise ValueError("nstart must be at least 1")

    repeated = list(dict.fromkeys(v for v, c in pd.Series(variables + active + v_patterns).value_counts().items()
                                  if False))
    for lst, name in ((variables, "variables"), (active, "active"), (v_patterns, "v_patterns")):
        if len(set(lst)) != len(lst):
            warnings.warn("variables named more than once are kept once: "
                          + ", ".join(sorted({v for v in lst if lst.count(v) > 1})))
    variables, active, v_patterns = (list(dict.fromkeys(x)) for x in (variables, active, v_patterns))

    x_axis = _axis_spec(xaxis, 1, "xaxis", data)
    y_axis = _axis_spec(yaxis, 2, "yaxis", data)
    f_axes = [a["factor"] for a in (x_axis, y_axis) if a["variable"] is None]
    n_extract = max([nfactors, 2] + f_axes)
    rng = np.random.default_rng(seed)

    auto_type = type is None
    type = "mca" if type is None else str(type).lower()
    if type not in ("mca", "pca"):
        raise ValueError("type must be 'mca' or 'pca'")
    if auto_type and all(_is_numeric(data[a]) for a in active):
        type = "pca"

    # ---- weights and complete cases
    if isinstance(weight, str):
        if weight not in data.columns:
            raise ValueError("weight column not found in data")
        variables = [v for v in variables if v != weight]
        active = [v for v in active if v != weight]
        v_patterns = [v for v in v_patterns if v != weight]
        weight = data[weight].to_numpy(dtype=float)
    elif weight is not None:
        weight = np.asarray(weight, dtype=float)
    if weight is not None and len(weight) != len(data):
        raise ValueError("Weights have not the correct dimensions")
    if not set(v_patterns) <= set(variables):
        raise ValueError("v_patterns must be a subset of variables")
    variables = list(dict.fromkeys(variables + [a["variable"] for a in (x_axis, y_axis) if a["variable"]]))
    keep = data[variables].notna().all(axis=1).to_numpy()
    if weight is not None:
        keep &= ~np.isnan(weight)
    if keep.sum() < len(data) / 2:
        nas = data[variables].isna().sum()
        nas = nas[nas > 0].sort_values(ascending=False, kind="stable")
        worst = ", ".join(f"{k} ({v})" for k, v in list(nas.items())[:5])
        what = (f" over the {len(variables)} variables" + (" and the weight" if weight is not None and np.isnan(weight).any() else "")
                + f"; those with most missing values: {worst}"
                + (f" and {len(nas) - 5} more" if len(nas) > 5 else "")
                + ". State variables to keep only those the analysis needs.")
        if keep.sum() < 2:
            raise ValueError(("only one case is" if keep.sum() else "no case is") + " complete" + what)
        warnings.warn(f"only {int(keep.sum())} of {len(data)} cases are complete" + what)
    if weight is not None:
        weight = weight[keep]
    case_weight = None if weight is None else weight.copy()
    D = data.loc[keep, variables]
    row_names = [str(i) for i in D.index] if not isinstance(D.index, pd.RangeIndex) else [str(i + 1) for i in D.index]
    row_names = [str(i + 1) for i in np.where(keep)[0]] if isinstance(data.index, pd.RangeIndex) else [str(i) for i in D.index]
    D = D.reset_index(drop=True)

    # ---- the analysis
    B = D.copy()
    for v in variables:
        if not _is_numeric(B[v]) and not isinstance(B[v].dtype, pd.CategoricalDtype):
            B[v] = pd.Categorical(_as_character(B[v]))
    if type == "mca":
        b = B[active]
        m = dichotomize(b.astype(object).where(b.notna(), None), variables=list(b.columns), sort=False, add=False, nas=None)
        if n_extract > m.shape[1] - 1:
            raise ValueError(f"only {m.shape[1] - 1} factors can be extracted from {m.shape[1]} categories")
        ff, eig = mca_rows(m.to_numpy(dtype=float), n_extract, weight)
        ff_names = [f"F{i + 1}" for i in range(n_extract)]
        vm = m.to_numpy(dtype=float)
    else:
        if len(active) < 2:
            raise ValueError("Only one quantitative active variable was found, so type=\"pca\" was chosen "
                             "automatically. Include at least one more active variable, or set type=\"mca\" explicitly."
                             if auto_type else "type=\"pca\" requires at least two active variables")
        b = D[active].astype(float)
        for v in [v for v in variables if v not in active]:
            if not _is_numeric(B[v]) and not isinstance(B[v].dtype, pd.CategoricalDtype):
                B[v] = pd.Categorical(_as_character(B[v]))
        if n_extract > b.shape[1]:
            raise ValueError(f"only {b.shape[1]} factors can be extracted from {b.shape[1]} active variables")
        ff, eig, vm = pca_scores(b.to_numpy(), weight, n_extract)
        ff_names = [f"PC{i + 1}" for i in range(n_extract)]
    # the latent model suits the type of analysis, which is only known now: LPA on the quantitative active
    # variables, LCA on the qualitative ones, both as they come
    latent_m = latent_clusters(l4c_meth, type)
    l4c_meth = [m_ for m_ in l4c_meth if m_ in ("kmeans", "pam", "hclust")]
    if no_k and (l4c_meth or latent_m):
        warnings.warn("clusters left out: nclusters states no number of groups to make")
        l4c_meth, latent_m = [], []
    l4c_data = m if type == "mca" else b.to_numpy(dtype=float)       # the active variables as they come
    latent_data = b
    latent_found = None
    n_use = max(1, int((eig[:nfactors] > crit_factors * eig.mean()).sum()))
    fm = ff[:, :n_use]

    def axis_values(ax):
        if ax["variable"] is None:
            return ff[:, ax["factor"] - 1]
        v = D[ax["variable"]]
        if isinstance(v.dtype, pd.CategoricalDtype):
            if v.cat.ordered:
                return v.cat.codes.to_numpy(dtype=float) + 1
            raise ValueError(f"{ax['variable']} is a nominal variable, so placing it on an axis would impose an "
                             "arbitrary order. Turn it into an ordered factor or into a numeric variable.")
        if not _is_numeric(v):
            if v.map(lambda x: isinstance(x, str)).all():
                raise ValueError(f"{ax['variable']} is a nominal variable, so placing it on an axis would impose an "
                                 "arbitrary order. Turn it into an ordered factor or into a numeric variable.")
            raise ValueError(f"{ax['variable']} must be numeric or an ordered factor to be plotted on an axis")
        return v.to_numpy(dtype=float)

    def axis_name(ax):
        return ff_names[ax["factor"] - 1] if ax["variable"] is None else ax["variable"]

    def axis_title(ax):
        if ax["variable"] is None:
            return f"{'F' if type == 'mca' else 'PC'}{ax['factor']}"
        label = getattr(data[ax["variable"]], "attrs", {}).get("label")
        return ax["variable"] if label is None else str(label)

    cc = np.column_stack([axis_values(x_axis), axis_values(y_axis)])
    cc_names = [axis_name(x_axis), axis_name(y_axis)]
    if cc_names[0] == cc_names[1]:
        cc_names[1] += ".1"
    pctvar = eig / eig.sum() * 100
    if type == "mca" and len(active) > 1:
        q = len(active)
        adj = np.where(eig > 1 / q, (q / (q - 1)) ** 2 * (eig - 1 / q) ** 2, 0.0)
        if len(adj) > 1 and adj[1] > 0:
            pctvar = adj / adj.sum() * 100
    for j, ax in enumerate((x_axis, y_axis)):
        scale = ax["scale"] if ax["scale"] is not None else (ax["variable"] is None and scale_axes)
        sd = cc[:, j].std(ddof=1)
        if scale and sd > 0:
            cc[:, j] = (cc[:, j] - cc[:, j].mean()) / sd

    cm = vm if cluster_on == "variables" else fm
    if scale_clusters is None:
        scale_clusters = not (cluster_on == "variables" and type == "mca")
    if scale_clusters:
        sd = cm.std(axis=0, ddof=1)
        sd = np.where((sd == 0) | np.isnan(sd), 1.0, sd)
        cm = (cm - cm.mean(axis=0)) / sd

    language = args.setdefault("language", "en")
    language = language[0] if isinstance(language, (list, tuple)) else language
    args["language"] = language
    groups_word, group_word = getByLanguage(groupsList, language), getByLanguage(groupList, language)

    # ---- k-means
    if auto_k:
        n_distinct = len(np.unique(np.round(cm, 8), axis=0))
        if stated_ks and any(k >= n_distinct for k in ks):
            warnings.warn(f"nclusters: {', '.join(str(k) for k in ks if k >= n_distinct)} groups left out, as there "
                          f"are only {n_distinct} distinct cases to cluster")
        ks = [k for k in ks if k < n_distinct]
        if not any(k >= 2 for k in ks):
            raise ValueError("there are too few distinct cases to choose a number of groups")
    fits = {}
    if ks:
        try:
            from sklearn.cluster import KMeans
        except ImportError:                                               # pragma: no cover
            raise ImportError("sur_scat needs scikit-learn for the k-means groups")
        for k in ks:
            seeds = rng.integers(0, 2 ** 31 - 1, size=nstart)
            runs = []
            for sd in seeds:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="KMeans is known to have a memory leak")
                    warnings.filterwarnings("ignore", message="Number of distinct clusters")
                    model = KMeans(n_clusters=k, n_init=1, init="random", random_state=int(sd)).fit(cm)
                runs.append((model.inertia_, model.labels_))
            runs.sort(key=lambda r: r[0])
            best_labels, best_sse = None, np.inf          # the best few starts are refined, the best result kept
            for _, labels0 in runs[:3]:
                labels0 = hartigan_refine(cm, labels0, k)
                sse = sum(((cm[labels0 == g] - cm[labels0 == g].mean(axis=0)) ** 2).sum() for g in range(k))
                if sse < best_sse - 1e-12:
                    best_labels, best_sse = labels0, sse
            fits[k] = best_labels

    k_index = None
    best = None
    k_all = list(ks)                      # the numbers of groups the index weighs: those of the solutions, and 1 for the BIC
    bic_model = None
    if auto_k and cluster_index == "bic":
        from .latent import latent_model, run_latent
        bic_model = latent_model(type)
        if weigh_one:
            k_all = [1] + ks
        latent_found = run_latent(latent_data, bic_model, k_all, seed=seed)
        if all(np.isnan(v) for v in latent_found["bic"].values()):
            warnings.warn(f"the BIC of {bic_model} could not be computed for any number of groups: "
                          "Calinski-Harabasz is used instead")
            cluster_index = "ch"
            k_all = list(ks)
    if auto_k:
        em = np.where(np.isfinite(vm), vm, 0.0)
        n = len(em)
        total = ((em - em.mean(axis=0)) ** 2).sum()
        one_wins = False

        def within(labels):
            return sum(((em[labels == g] - em[labels == g].mean(axis=0)) ** 2).sum() for g in np.unique(labels))

        if cluster_index == "bic":        # the lower the better, and defined for a single group as well
            values = [latent_found["bic"][k] for k in k_all]
            best_all = int(np.nanargmin(values))
            one_wins = len(k_all) > len(ks) and best_all == 0
            best = None if one_wins else best_all - (len(k_all) - len(ks))     # no solution is the main one when 1 wins
        else:
            sample = np.sort(rng.choice(n, 2000, replace=False)) if n > 2000 else np.arange(n)
            values = []
            for k in ks:
                if cluster_index == "ch":
                    w = within(fits[k])
                    values.append(((total - w) / (k - 1)) / (w / (n - k)))
                else:
                    from sklearn.metrics import silhouette_score
                    labels = fits[k][sample]
                    values.append(float(silhouette_score(em[sample], labels)) if len(set(labels)) > 1 else np.nan)
            best = best_all = int(np.nanargmax(values))
        k_index = pd.DataFrame({"k": k_all, cluster_index: values, "optimal": np.arange(len(k_all)) == best_all})
        index_name = {"ch": "Calinski-Harabasz", "silhouette": "silhouette", "bic": f"BIC of {bic_model}"}[cluster_index]
        k_index.attrs["index"] = index_name
        if one_wins:
            log.info('nclusters="auto": 1 group (no clustering), by %s; every case takes the same color', index_name)
        else:
            log.info('nclusters="auto": %d groups, by %s', ks[best], index_name)
        j_order = list(range(len(ks))) if keep_all_k else [best]
    else:
        j_order = list(range(len(ks)))
    ks_used = [ks[j] for j in j_order]
    g_cols: list[str] = []
    g_main = None
    for j in j_order:
        k = ks[j]
        labels_k = fits[k] + 1
        if sort_clusters:
            means = []
            for c in range(1, k + 1):
                sel = labels_k == c
                means.append(np.average(cc[sel, 0], weights=None if weight is None else weight[sel]))
            rank = np.argsort(means, kind="stable")
            labels_k = np.array([int(np.where(rank == c - 1)[0][0]) + 1 for c in labels_k])
        name_g = f"{groups_word}({k:0{len(str(max(ks_used)))}d})"
        cats = [f"{group_word}: {i:0{len(str(k))}d}" for i in range(1, k + 1)]
        B[name_g] = pd.Categorical([cats[c - 1] for c in labels_k], categories=cats, ordered=True)
        g_cols.append(name_g)
        if not auto_k or j == best:
            g_main = name_g

    # ---- patterns
    idx = None
    shares = None
    if patterns:
        key = D[v_patterns].astype(str).agg("\r".join, axis=1)
        _, first_pos = np.unique(key, return_index=True)
        seen: dict = {}
        idx = np.array([seen.setdefault(k, len(seen) + 1) for k in key])
        groups = [np.where(idx == g)[0] for g in range(1, idx.max() + 1)]

        def wmean(x, w):
            return float(np.mean(x)) if w is None else float((x * w).sum() / w.sum())

        shares = {v: cluster_shares(idx, len(groups), B[v], weight) for v in g_cols}
        firsts = [g[0] for g in groups]
        new_b = B.iloc[firsts].reset_index(drop=True)
        for v in [c for c in B.columns if c not in v_patterns]:
            if v in g_cols:
                ordered = []
                for g in groups:
                    counts = B[v].iloc[g].astype(str).value_counts()
                    names = sorted(counts.index, key=lambda s: (-counts[s], _collation_key(s)))
                    ordered.append("|".join(names))
                new_b[v] = pattern_groups(ordered, list(B[v].cat.categories))
            elif _is_numeric(data[v]):
                new_b[v] = [wmean(B[v].to_numpy(dtype=float)[g], None if weight is None else weight[g]) for g in groups]
            else:
                cells = []
                for g in groups:
                    values = B[v].iloc[g].astype(str).to_numpy()
                    w_g = np.ones(len(g)) if weight is None else weight[g]
                    totals: dict = {}
                    for val, w in zip(values, w_g):
                        totals[val] = totals.get(val, 0.0) + w
                    names = sorted(totals, key=lambda s: (-totals[s], _collation_key(s)))
                    cells.append("|".join(names))
                new_b[v] = cells
        B = new_b
        cc = np.column_stack([[wmean(cc[g, j], None if weight is None else weight[g]) for g in groups]
                              for j in range(2)])
        weight = np.array([float(len(g)) if weight is None else float(weight[g].sum()) for g in groups])
        row_names = [f"{i + 1:0{len(str(len(B)))}d}" for i in range(len(B))]
        n_name = "n"
        while n_name in B.columns:
            n_name += ".1"
        B[n_name] = weight
        args.setdefault("size", n_name)
        text = [""] * len(B)
        for v in v_patterns:
            text = [t + f"<b>{v}</b>: {s}<br/>" for t, s in zip(text, B[v].astype(str))]
        text = [t + f"<b>{getByLanguage(casesList, language)}</b>: {int(round(n))}" for t, n in zip(text, B[n_name])]
        t_name = "ntext"
        while t_name in B.columns:
            t_name += ".1"
        B[t_name] = text
        args.setdefault("ntext", t_name)

    name = args["name"] = name_by_language(None, language, None)
    B[name] = row_names
    stated = v_patterns if stated_v_patterns else []
    blocks = {"_vpattern": stated, "_actives": active,
              "_others": [v for v in variables if v not in stated and v not in active], "_nclusters": g_cols}
    n_col = [n_name] if patterns else []
    order = (stated + (n_col if stated_v_patterns else []) + active + ([] if stated_v_patterns else n_col)
             + blocks["_others"] + g_cols + ([t_name] if patterns else []) + [name])
    if columns is not None:
        expanded = []
        for c in _as_list(columns):
            expanded += blocks[c.lower()] if c.lower() in blocks else [c]
        unknown = [c for c in expanded if c not in B.columns]
        if unknown:
            warnings.warn(f"{', '.join(unknown)} is/are not a column of the node data frame.")
        order = expanded + order
    order = list(dict.fromkeys(order))
    B = B[[c for c in order if c in B.columns] + [c for c in B.columns if c not in order]]

    # ---- a sample of the nodes, jitter
    n_all = len(B)
    kept = None
    if n_all > max_n:
        kept = np.sort(rng.choice(n_all, max_n, replace=False))
        B = B.iloc[kept].reset_index(drop=True)
        cc = cc[kept]
        if shares is not None:
            shares = {k: v.iloc[kept].reset_index(drop=True) for k, v in shares.items()}
    if jitter is True:
        jitter = 0.02
    if jitter > 0:
        spread = cc.max(axis=0) - cc.min(axis=0)
        cc = cc + (rng.uniform(-0.5, 0.5, size=cc.shape)) * jitter * spread

    args["nodes"] = B
    args["layout"] = pd.DataFrame(cc, columns=cc_names)
    if "color" not in args and g_main is not None:
        args["color"] = g_main
    args["frequencies"] = True
    args["show_axes"] = True
    if "axes_labels" not in args:
        args["axes_labels"] = [f"{axis_title(ax)} ({pctvar[ax['factor'] - 1]:.1f}%)" if ax["variable"] is None
                               else axis_title(ax) for ax in (x_axis, y_axis)]
    args["show_coordinates"] = True
    args["statistics"] = True
    args.pop("degree_filter", None)
    args.setdefault("label", "")
    args.setdefault("controls", (1, 2, 4))
    net = netCoin(**args)
    net.scat = {}
    if not patterns and "ntext" not in args:
        net.options["nodeTextDynamic"] = "scat"
    if x_axis["variable"] is None and y_axis["variable"] is None and x_axis["factor"] == 1 and y_axis["factor"] == 2:
        net.scat["plane_name"] = type.upper()
    if idx is not None and kept is None:
        net.scat["case_to_pattern"] = idx
        if case_weight is not None:
            net.scat["case_weight"] = case_weight
    if kept is not None:
        net.scat["sampled_nodes"] = kept + 1
        net.scat["sampled_from"] = n_all
    net.clusters, net.cluster_means = {}, None
    if g_cols:
        net.scat["cluster_columns"] = g_cols
        if shares is not None:
            net.scat["cluster_shares"] = shares
        net.clusters = current_clusters(net)
        net.cluster_means = current_cluster_means(net)
    # The planes and groups of looking4clusters are added last, so that asking for them leaves everything
    # above as it was, the sample drawn by max_n and the jitter included; the classes of the latent model
    # travel the same way
    if l4c_lay or rot_lay or l4c_meth or latent_m:
        from .rotations import rotation_axes, rotation_name
        found = (l4c_run(l4c_data, l4c_lay, l4c_meth, ks_used, seed=seed, force_execution=force_execution)
                 if (l4c_lay or l4c_meth) else {"axes": {}, "clusters": {}})
        text_layouts = [layouts] if isinstance(layouts, str) else (
            list(layouts) if isinstance(layouts, (list, tuple)) else None)
        rotated = rotation_axes(ff, vm, eig, case_weight, nfactors, rot_lay, type,
                                quiet=text_layouts is not None and "all" in [str(x).lower() for x in text_layouts],
                                asked=rotation_name(text_layouts) if text_layouts is not None else None)
        found["axes"] = {**rotated, **found["axes"]}    # the rotated planes come first, as they are of the very components
        if latent_m:
            from .latent import run_latent
            if latent_found is None or cluster_index != "bic":
                latent_found = run_latent(latent_data, latent_m[0], ks_used, seed=seed)
            width = len(str(max(ks_used)))
            for k in ks_used:
                if k in latent_found["clusters"]:
                    found["clusters"][f"{latent_m[0]}({k:0{width}d})"] = latent_found["clusters"][k]
        net = add_l4c(net, found, idx=idx, kept=None if kept is None else kept + 1, case_weight=case_weight,
                      sort=sort_clusters, suffix=suffix, align_to=align_to,
                      axes_labels="keep" if stated_labels else None)
    if k_index is not None:
        net.scat["cluster_index"] = k_index
    return net
