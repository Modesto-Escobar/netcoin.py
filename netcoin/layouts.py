"""Node coordinates computed from coincidences (ports of netCoin's layouts)."""
from __future__ import annotations

import math
import warnings
from typing import Sequence

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

    Parameters
    ----------
    incidences : DataFrame
        One row per case, one 0/1 column per variable.

    Returns
    -------
    DataFrame
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

    Parameters
    ----------
    coincidences : Coin or DataFrame
        The coincidences, or a correlation matrix.

    Returns
    -------
    DataFrame
        ``F1`` and ``F2``, one row per variable.
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


def _sort_loadings(loadings: np.ndarray) -> np.ndarray:
    """Order the factors by their sum of squares and make each one mostly positive (R's ``sortLoadings``)."""
    order = np.argsort(-(loadings ** 2).sum(axis=0), kind="stable")
    loadings = loadings[:, order]
    flip = loadings.sum(axis=0) < 0
    loadings[:, flip] = -loadings[:, flip]
    return loadings


def _factanal_ml(corr: np.ndarray, factors: int) -> np.ndarray:
    """Maximum likelihood factor loadings of a correlation matrix (R's ``factanal.fit.mle``)."""
    from scipy.optimize import minimize

    p = corr.shape[0]

    def eig(psi):
        scaled = corr / np.sqrt(np.outer(psi, psi))
        values, vectors = np.linalg.eigh(scaled)
        order = np.argsort(values)[::-1]
        return values[order], vectors[:, order]

    def loadings_of(psi):
        values, vectors = eig(psi)
        load = vectors[:, :factors] * np.sqrt(np.maximum(values[:factors] - 1, 0))
        return np.sqrt(psi)[:, None] * load

    def objective(psi):
        values = eig(psi)[0][factors:]
        return -(np.sum(np.log(values) - values) - factors + p)

    def gradient(psi):
        load = loadings_of(psi)
        g = load @ load.T + np.diag(psi) - corr
        return np.diag(g) / psi ** 2

    scale = 0.01                                    # R optimises with parscale = 0.01
    start = (1 - 0.5 * factors / p) / np.diag(np.linalg.inv(corr))
    result = minimize(lambda u: objective(u * scale), start / scale, jac=lambda u: gradient(u * scale) * scale,
                      method="L-BFGS-B", bounds=[(0.005 / scale, 1 / scale)] * p, options={"gtol": 0, "maxiter": 1000})
    if not result.success and result.status != 2:
        raise ValueError("unable to optimize the factor analysis")
    return loadings_of(result.x * scale)


def _varimax(x: np.ndarray, eps: float = 1e-5) -> np.ndarray:
    """Varimax rotation with Kaiser normalisation (R's ``stats::varimax``)."""
    p, k = x.shape
    scale = np.sqrt((x ** 2).sum(axis=1))
    x = x / scale[:, None]
    rotation, d = np.eye(k), 0.0
    for _ in range(1000):
        z = x @ rotation
        b = x.T @ (z ** 3 - z @ np.diag((z ** 2).sum(axis=0)) / p)
        u, s, vt = np.linalg.svd(b)
        rotation = u @ vt
        previous, d = d, s.sum()
        if d < previous * (1 + eps):
            break
    return (x @ rotation) * scale[:, None]


def _oblimin(a: np.ndarray, eps: float = 1e-5, maxit: int = 2000) -> np.ndarray:
    """Direct oblimin (quartimin) rotation by gradient projection (the ``GPArotation`` default)."""
    k = a.shape[1]

    def criterion(loadings):
        weights = (loadings ** 2).sum(axis=1, keepdims=True) - loadings ** 2
        return (loadings ** 2 * weights).sum() / 4, loadings * weights

    t = np.eye(k)
    inverse = np.linalg.inv(t)
    loadings = a @ inverse.T
    f, gq = criterion(loadings)
    grad = -(loadings.T @ gq @ inverse).T
    alpha = 1.0
    for _ in range(maxit + 1):
        gp = grad - t @ np.diag((t * grad).sum(axis=0))
        s = np.sqrt((gp ** 2).sum())
        if s < eps:
            break
        alpha *= 2
        for _ in range(11):
            x = t - alpha * gp
            candidate = x / np.sqrt((x ** 2).sum(axis=0))
            inverse_c = np.linalg.inv(candidate)
            loadings = a @ inverse_c.T
            f_new, gq = criterion(loadings)
            if f - f_new > 0.5 * s ** 2 * alpha:
                break
            alpha /= 2
        t, f = candidate, f_new
        inverse = inverse_c
        grad = -(loadings.T @ gq @ inverse).T
    return loadings


_FACT_NAMES = {"pc": "Principal Components", "nf": "Non-rotated Factorial",
               "vf": "Varimax factorial", "of": "Oblimin factorial"}


def layout_fact(data: pd.DataFrame, method="pc"):
    """Coordinates of the variables from their loadings on two components or factors.

    Port of ``layoutFact``. ``method`` is any of ``"pc"`` (principal components of the
    correlations), ``"nf"`` (maximum likelihood factors, not rotated), ``"vf"`` (varimax) and
    ``"of"`` (oblimin), or ``"null"`` for no layout. The cases with a missing value are left
    out. One method gives a DataFrame (``Comp.1``, ``Comp.2`` or ``Factor1``, ``Factor2``, one
    row per variable); several give a dict of them named as in R. The sign of each axis is
    arbitrary, as it is in R.

    Parameters
    ----------
    data : DataFrame
        One numeric column per variable.
    method : str or list of str
        ``"pc"``, ``"nf"``, ``"vf"``, ``"of"`` or ``"null"``.

    Returns
    -------
    DataFrame or dict of DataFrame
    """
    methods = [method] if isinstance(method, str) else list(method)
    unknown = sorted(set(methods) - {"null", "pc", "nf", "vf", "of"})
    if unknown:
        warnings.warn(", ".join(unknown) + " method is not implemented. ")
    if methods[0] == "null":
        return None
    clean = data.dropna()
    names = list(clean.columns)
    corr = np.corrcoef(clean.to_numpy(dtype=float), rowvar=False)
    out = {}
    if "pc" in methods:
        values, vectors = np.linalg.eigh(corr)
        order = np.argsort(values)[::-1][:2]
        out[_FACT_NAMES["pc"]] = pd.DataFrame(vectors[:, order], index=names, columns=["Comp.1", "Comp.2"])
    if {"nf", "vf", "of"} & set(methods):
        load = _sort_loadings(_factanal_ml(corr, 2))
        rotations = {"nf": None, "vf": _varimax, "of": _oblimin}
        for key in ("nf", "vf", "of"):
            if key in methods:
                rotated = load if rotations[key] is None else _sort_loadings(rotations[key](load))
                out[_FACT_NAMES[key]] = pd.DataFrame(rotated, index=names, columns=["Factor1", "Factor2"])
    if not out:
        return None
    return next(iter(out.values())) if len(out) == 1 else out


def layout_special(layout, special: dict):
    """Replace the layout names that only the calling function can compute.

    Port of ``layoutSpecial``. ``special`` maps a two-letter code (``"mc"``, ``"pc"``...) to a
    function without arguments that returns the coordinates. A single name comes back as the
    coordinates; among several (a list of names or a dict), the coordinates take the place of
    the name and keep it as the name of the plane. Anything else is left as it is, to be
    computed later (``"fr"``, ``"circle"``...).
    """
    if isinstance(layout, str):
        items, single = [(None, layout)], True
    elif isinstance(layout, (list, tuple)) and all(isinstance(x, str) for x in layout):
        items, single = [(None, x) for x in layout], len(layout) == 1
    elif isinstance(layout, dict):
        items, single = list(layout.items()), False
    else:
        return layout
    out = {}
    for label, value in items:
        if isinstance(value, str) and value[:2].lower() in special:
            label = label or value
            value = special[value[:2].lower()]()
        elif label is None:
            label = value if isinstance(value, str) else str(len(out) + 1)
        out[label] = value
    return next(iter(out.values())) if single else out


# ---------------------------------------------------------- circle and grid ----
def layout_circle(nodes: pd.DataFrame, items: Sequence | None = None, deg: float = 0,
                  name: str | None = None) -> pd.DataFrame:
    """Place nodes on a circle (rD3plot's ``layoutCircle``).

    Parameters
    ----------
    nodes : DataFrame
        The table of nodes (one row each).
    items : names or 0-based positions, optional
        The nodes to place (all by default); the rest get ``NaN``. They are spread
        evenly in the order given, the first at angle ``deg``.
    deg : float
        Angle in degrees of the first node.
    name : str, optional
        Column with the names of the nodes (otherwise ``"1"``, ``"2"``...).

    Returns
    -------
    DataFrame
        ``x`` and ``y`` (rounded to 3 decimals), one row per node.
    """
    names = [str(v) for v in nodes[name]] if name else [str(i + 1) for i in range(len(nodes))]
    out = pd.DataFrame({"x": np.nan, "y": np.nan}, index=names)
    if items is None:
        chosen = list(range(len(names)))
    elif all(isinstance(i, (int, np.integer)) for i in items):
        chosen = [i for i in items if 0 <= i < len(names)]
    elif all(isinstance(i, str) for i in items):
        chosen = [names.index(i) for i in items if i in names]
    else:
        warnings.warn("nodes: must be numeric or character")
        return out
    if len(chosen) < 1:
        warnings.warn("nodes: length of items present in 'N' must be greater than zero")
    elif len(chosen) == 1:
        out.iloc[chosen[0]] = 0.0
    else:
        k = len(chosen)
        angle = np.arange(1, k + 1) / (k / 2) * math.pi
        angle = angle - angle[0] + deg * math.pi / 180
        out.iloc[chosen, 0] = np.round(np.cos(angle), 3)
        out.iloc[chosen, 1] = np.round(np.sin(angle), 3)
    return out


def layout_grid(nodes: pd.DataFrame, string: str, name: str | None = None,
                byrow: bool = False) -> pd.DataFrame:
    """Place nodes on a grid described by a text (rD3plot's ``layoutGrid``).

    ``string`` lists the columns of the grid separated by ``"."`` and the nodes of each
    column from top to bottom separated by ``","``; ``";"`` leaves an empty cell and a
    column that ends with ``","`` leaves one at the bottom. A column with fewer nodes
    is centred. For example ``"a,b.c.d,e,f"``. With ``byrow`` the grid is read by rows.

    Returns a DataFrame ``x``, ``y`` in [0, 1] (``NaN`` for the nodes that are not in the
    text).

    Parameters
    ----------
    nodes : DataFrame
        The nodes (their ``name`` column or ``name`` is matched against the text).
    string : str
        The description of the grid.
    name : str, optional
        Column of ``nodes`` with the names.
    byrow : bool
        Read the grid by rows.
    """
    if not isinstance(string, str):
        warnings.warn("string: must be character")
        string = ""
    names = [str(v) for v in nodes[name]] if name else [str(i + 1) for i in range(len(nodes))]
    where = {n: i for i, n in enumerate(names)}
    x = np.full(len(names), np.nan)
    y = np.full(len(names), np.nan)

    levels = string.split(".") if string else []
    xlen = len(levels)
    ylen = 0
    parsed = []
    for i, level in enumerate(levels):
        a = level.replace("*", "").replace(";", ",,")
        add_to_end = a.endswith(",")
        parts = a.split(",") if a else []
        if parts and parts[-1] == "":          # strsplit() drops one trailing empty piece
            parts.pop()
        if add_to_end:
            parts.append("")
        ylen = max(ylen, len(parts))
        for j, item in enumerate(parts):
            if item != "" and item in where:
                x[where[item]] = i
                y[where[item]] = j
        parsed.append(parts)
    for parts in parsed:
        if len(parts) != ylen:
            for item in parts:
                if item != "" and item in where:
                    y[where[item]] += (ylen - len(parts)) / 2

    placed = ~np.isnan(x)
    with np.errstate(all="ignore"):
        gx = x / (xlen - 1)
        gy = 1 - y / (ylen - 1)
    if byrow:
        gx, gy = 1 - gy, 1 - gx
    for arr in (gx, gy):
        arr[placed & np.isnan(arr)] = 0.0          # 0/0 of a one-column or one-row grid
    return pd.DataFrame({"x": gx, "y": gy}, index=names)


# ------------------------------------------------------------ force atlas 2 ----
def force_atlas2(adjacency: np.ndarray, iterations: int = 100, linlog: bool = False,
                 pos: np.ndarray | None = None, nohubs: bool = False, k: float = 400,
                 gravity: float = 1, ks: float = 0.1, ksmax: float = 10, delta: float = 1,
                 center: Sequence[float] | None = None, tolerance: float = 0.1, dim: int = 2,
                 seed: int | None = None) -> np.ndarray:
    """ForceAtlas2 layout (port of rD3plot's ``layout.force.atlas.2``).

    ``adjacency`` is the (weighted) adjacency matrix. Without ``pos`` the start is
    random, so pass ``seed`` for a reproducible layout.
    """
    coincidences = np.asarray(adjacency, dtype=float)
    center = np.zeros(dim) if center is None else np.asarray(center, dtype=float)
    n = coincidences.shape[0]
    binary = (coincidences != 0).astype(float)
    deg = binary.sum(axis=1)

    rng = np.random.default_rng(seed)
    if pos is None:
        difference = 2000 / (n * dim)
        grid = np.arange(-1000, 1000 + difference / 2, difference)
        position = rng.choice(grid, size=(n, dim), replace=len(grid) < n * dim).astype(float)
    else:
        position = np.array(pos, dtype=float)
    at_center = np.all(position == center, axis=1)
    position[at_center] = center + 0.01

    forces1 = np.zeros((dim, n))
    for _ in range(iterations):
        forces2 = forces1
        diff = position[:, None, :] - position[None, :, :]
        distances = np.sqrt((diff ** 2).sum(axis=2))
        distances[distances < 0.01] = 0.01
        units = [diff[:, :, j] / distances for j in range(dim)]
        repulsion = k * ((deg[None, :] + 1) * (deg[:, None] + 1)) / distances
        attraction = np.log(1 + distances) if linlog else distances
        attraction = (coincidences ** delta) * attraction
        if nohubs:
            attraction = attraction / (deg[:, None] + 1)
        far = np.vstack([np.nansum(u * (repulsion - attraction), axis=1) for u in units])
        to_center = center[None, :] - position
        norm = np.sqrt((to_center ** 2).sum(axis=1))
        uv2 = (to_center / norm[:, None]).T
        gravity_force = uv2 * np.tile(gravity * (deg + 1), (dim, 1))
        forces1 = np.round(far + gravity_force, 2)

        swing = np.abs(np.sqrt(((forces1 - forces2) ** 2).sum(axis=0)))
        global_swing = np.sum((deg + 1) * swing)
        if np.all(swing == 0):
            break
        tra = np.abs(np.sqrt(((forces1 + forces2) ** 2).sum(axis=0))) / 2
        global_tra = np.sum((deg + 1) * tra)
        global_speed = tolerance * global_tra / global_swing
        speed = ks * global_speed / (1 + global_speed * np.sqrt(swing))
        speed_constrain = ksmax / np.abs(np.sqrt((forces1 ** 2).sum(axis=0)))
        speed = np.where(speed >= speed_constrain, speed_constrain, speed)
        position = position + (forces1 * np.tile(speed, (dim, 1))).T
    return position


# --------------------------------------------------- igraph: layouts by name ----
LAYOUT_CODES = ("da", "dr", "ci", "fo", "fr", "ge", "gr", "ka", "lg", "md", "ra", "re", "st", "su")
COMMUNITY_CODES = ("ed", "fa", "la", "le", "lo", "op", "sp", "wa")


def require_igraph():
    try:
        import igraph
    except ImportError as error:                      # pragma: no cover
        raise ImportError("this needs python-igraph: pip install netcoin[igraph]") from error
    return igraph


def layout_code(layout: str) -> str | None:
    """The two-letter code of a layout name (``"layout.fruchterman.reingold"``, ``"kamada"``...)."""
    code = str(layout).replace("layout.", "").replace("layout_", "").lower()[:2]
    if code in LAYOUT_CODES:
        return code
    warnings.warn(f"{layout} is not a valid layout")
    return None


def community_code(method: str) -> str | None:
    code = str(method).replace("cluster_", "").replace("community_", "").lower()[:2]
    if code in COMMUNITY_CODES:
        return code
    warnings.warn(f"{method} is not a valid community detection method")
    return None


def _flip_y(layout):
    """R's ``layout_as_tree`` flips y by default so that the root is on top."""
    coords = np.array(layout.coords, dtype=float)
    coords[:, 1] = coords[:, 1].max() - coords[:, 1] if len(coords) else coords[:, 1]
    return type("Flipped", (), {"coords": coords})()


def compute_layout(graph, code: str, criteria: str | None = None, seed: int | None = None) -> np.ndarray:
    """Coordinates of the vertices of an ``igraph.Graph`` with the layout of a two-letter code.

    ``da`` davidson-harel, ``dr`` drl, ``ci`` circle, ``fo`` force atlas 2 (weighted by the
    ``criteria`` edge attribute), ``fr`` fruchterman-reingold, ``gr`` grid, ``ka``
    kamada-kawai, ``lg`` lgl, ``md`` mds, ``ra`` random, ``re`` reingold-tilford, ``st``
    star, ``su`` sugiyama. (``ge`` gem is not in python-igraph.)
    """
    require_igraph()
    weights = "weight" if "weight" in graph.es.attributes() else None
    if seed is not None:
        import random
        random.seed(seed)
    if code == "fo":
        n = graph.vcount()
        adj = np.zeros((n, n))
        for e in graph.es:
            w = e[criteria] if criteria and criteria in graph.es.attributes() else 1.0
            adj[e.source, e.target] = adj[e.target, e.source] = 1.0 if w is None else w
        return force_atlas2(adj, seed=seed)
    if code == "ge":
        raise NotImplementedError("the 'gem' layout is not available in python-igraph")
    method = {
        "da": lambda: graph.layout_davidson_harel(),
        "dr": lambda: graph.layout_drl(weights=weights),
        "ci": lambda: graph.layout_circle(),
        "fr": lambda: graph.layout_fruchterman_reingold(weights=weights),
        "gr": lambda: graph.layout_grid(),
        "ka": lambda: graph.layout_kamada_kawai(weights=weights),
        "lg": lambda: graph.layout_lgl(),
        "md": lambda: graph.layout_mds(),
        "ra": lambda: graph.layout_random(),
        "re": lambda: _flip_y(graph.layout_reingold_tilford()),
        "st": lambda: graph.layout_star(),
        "su": lambda: graph.layout_sugiyama(),
    }[code]
    return np.array(method().coords, dtype=float)[:graph.vcount()]


def compute_membership(graph, code: str) -> list[int]:
    """Community of every vertex (1, 2, ...) of an ``igraph.Graph`` with a two-letter method."""
    require_igraph()
    weights = "weight" if "weight" in graph.es.attributes() else None
    if code == "ed":
        # R's result is the unweighted one: edge betweenness reads weights as distances
        clusters = graph.community_edge_betweenness().as_clustering()
    elif code == "fa":
        clusters = graph.community_fastgreedy(weights=weights).as_clustering()
    elif code == "la":
        clusters = graph.community_label_propagation(weights=weights)
    elif code == "le":
        clusters = graph.community_leading_eigenvector(weights=weights)
    elif code == "lo":
        clusters = graph.community_multilevel(weights=weights)
    elif code == "op":
        clusters = graph.community_optimal_modularity(weights=weights)
    elif code == "sp":
        clusters = graph.community_spinglass(weights=weights)
    else:
        clusters = graph.community_walktrap(weights=weights).as_clustering()
    return [m + 1 for m in clusters.membership]
