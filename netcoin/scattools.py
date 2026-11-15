"""What can be added to a scattergram: planes (``add_axes``), groups (``add_clusters``) and their alignment.

Ports of ``addAxes``, ``addClusters``, ``replaceClusters`` and ``layoutsAlign`` of netCoin. All of
them take a :class:`netCoin` built by :func:`sur_scat` and return a modified copy.
"""
from __future__ import annotations

import re
import warnings
from copy import deepcopy
from typing import Sequence

import numpy as np
import pandas as pd

from .core import netCoin
from .dichotomize import _collation_key
from .surscat import (cluster_shares, current_cluster_means, current_clusters, current_layouts)
from .utils import getByLanguage, groupList, restList


# ---------------------------------------------------------------- small helpers ----
def name_clash(asked: Sequence[str], given: Sequence[str], shown: int = 5) -> str:
    """``"a" (added again as "a.1"), ...`` for the warning of a name already taken."""
    pairs = [f'"{a}" (added again as "{g}")' for a, g in list(zip(asked, given))[:shown]]
    out = ", ".join(pairs)
    if len(asked) > shown:
        out += f" and {len(asked) - shown} more"
    return out


def sampled_note(net) -> str:
    """The sentence that explains a size mismatch when ``max_n`` drew a sample of the nodes."""
    sampled = net.scat.get("sampled_from")
    if sampled is None:
        return ""
    return (f" These {len(net.nodes)} nodes were drawn at random out of {sampled}, as max_n capped them, so anything "
            "computed on the whole data no longer lines up with them. Compute it on the node table of this object "
            "instead, which holds the sample itself, or raise max_n so that every case is drawn.")


def make_unique(existing: Sequence[str], name: str) -> str:
    """``make.unique(c(existing, name))`` for the last element."""
    if name not in existing:
        return name
    k = 1
    while f"{name}.{k}" in existing:
        k += 1
    return f"{name}.{k}"


def _copy(net):
    new = deepcopy(net)
    if not hasattr(new, "scat"):
        new.scat = {}
    return new


# ------------------------------------------------------------------- the planes ----
def as_plane(z, name: str = "plane", who: str = "plane") -> np.ndarray:
    z = np.asarray(z.to_numpy() if isinstance(z, pd.DataFrame) else z, dtype=float)
    if z.ndim != 2 or z.shape[1] == 0:
        raise ValueError(f"{who}: {name} holds no column of coordinates")
    return z


def norm_plane(z: np.ndarray, weights=None, who: str = "normPlane"):
    """Centre the plane and divide it by its root mean square radius (one divisor for the whole plane)."""
    if weights is None:
        centre = z.mean(axis=0)
        z = z - centre
        scale = np.sqrt((z ** 2).sum() / len(z))
    else:
        total = float(np.sum(weights))
        if total <= 0:
            raise ValueError(f"{who}: the weights add up to {total}")
        centre = (z * weights[:, None]).sum(axis=0) / total
        z = z - centre
        scale = np.sqrt(((z ** 2).sum(axis=1) * weights).sum() / total)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(f"{who}: the plane has no size (every node sits on the same point), so it cannot be brought "
                         "to a common scale")
    return z / scale, centre, scale


def procrustes_rotate(a: np.ndarray, b: np.ndarray, weights=None):
    """The orthogonal matrix ``r`` that brings ``a @ r`` closest to ``b`` (reflections allowed), and ``a @ r``."""
    m = a.T @ b if weights is None else (a * weights[:, None]).T @ b
    u, _, vt = np.linalg.svd(m)
    r = u @ vt
    return a @ r, r


def rotated_axis_titles(names, r: np.ndarray, suffix: str = ""):
    """Titles of the axes of a plane once turned by ``r``, followed by ``suffix``; percentages are recomputed."""
    if names is None:
        return None
    pattern = re.compile(r" *\(([-0-9.]+)%\) *$")
    base = [pattern.sub("", str(c)) for c in names]
    out = [b + suffix if b else "" for b in base]
    if len(names) == r.shape[1] and all(pattern.search(str(c)) for c in names):
        pct = np.array([float(pattern.search(str(c)).group(1)) for c in names])
        shares = (r ** 2 * pct[:, None]).sum(axis=0)
        out = [f"{t} ({s:.1f}%)" for t, s in zip(out, shares)]
    return out


def layout_weights(net, weights, n: int, who: str):
    """The weight of every node: what ``layoutWeights`` settles (cases per pattern by default)."""
    idx = net.scat.get("case_to_pattern")
    n_cases = None if idx is None else len(idx)
    case_w = net.scat.get("case_weight")
    if case_w is not None and n_cases is not None and len(case_w) != n_cases:
        case_w = None
    if weights is None:
        if idx is None:
            return None, "node"
        cw = np.ones(n_cases) if case_w is None else np.asarray(case_w, dtype=float)
        return np.bincount(np.asarray(idx) - 1, weights=cw, minlength=n), \
            "cases per pattern" if case_w is None else "cases per pattern (caseWeight)"
    if isinstance(weights, str):
        if weights not in net.nodes.columns:
            raise ValueError(f'{who}: there is no node column called "{weights}"')
        w, source = net.nodes[weights].to_numpy(dtype=float), f'column "{weights}"'
    else:
        w, source = np.asarray(weights, dtype=float), "weights"
    if w.size == 1:
        w = np.full(n, float(w))
    if len(w) != n:
        if n_cases is not None and len(w) == n_cases:
            w = np.bincount(np.asarray(idx) - 1, weights=w, minlength=n)
            source += ", added up by pattern"
        else:
            raise ValueError(f"{who}: weights holds {len(w)} values and there are {n} nodes"
                             + (f" (or {n_cases} cases)" if n_cases is not None else "") + sampled_note(net))
    if np.all(w == 1):
        return None, "node"
    return w, source


def plane_name(to, available: list[str], who: str) -> str:
    if isinstance(to, (int, np.integer)):
        if to < 1 or to > len(available):
            raise ValueError(f"{who}: to must be an index between 1 and {len(available)}")
        return available[to - 1]
    if to not in available:
        raise ValueError(f'{who}: there is no plane called "{to}". Available: ' + ", ".join(available))
    return to


def layouts_align(net, align: str = "procrustes", to=1, gpa: bool = False, weights=None, which=None,
                  suffix: str | None = None, tol: float = 1e-10, maxit: int = 100, who: str = "layoutsAlign"):
    """Lay the planes of a scattergram over one another (``layoutsAlign``).

    ``align`` is ``"procrustes"`` (centre, equal size and best rotation or reflection), ``"cloud"``
    (centre and size only), ``"canvas"`` (each axis to [0, 1]) or ``"none"``. ``to`` is the plane the
    rest are brought to (name or 1-based index); with ``gpa`` every plane is rotated towards the mean of
    them all. With ``suffix`` the originals stay and the aligned planes are added, right after their
    own, under their name and the suffix; otherwise they take the place of the originals.

    Parameters
    ----------
    net : netCoin
        A scattergram with more than one plane.
    align : {"procrustes", "cloud", "canvas", "none"}
        How the planes are laid over one another.
    to : int or str
        The plane the rest are brought to (1-based position or name).
    gpa : bool
        Rotate every plane towards the mean of them all (generalized Procrustes).
    weights : sequence or str, optional
        Weight of each node in the fit.
    which : sequence, optional
        The planes to align (all by default).
    suffix : str, optional
        Keep the originals and add the aligned planes under their name followed by ``suffix``.
    tol, maxit : float, int
        Convergence of the generalized Procrustes iteration.

    Returns
    -------
    netCoin
        A modified copy.
    """
    if align not in ("procrustes", "cloud", "canvas", "none"):
        raise ValueError(f"{who}: align must be 'procrustes', 'cloud', 'canvas' or 'none'")
    if suffix is not None and (not isinstance(suffix, str) or not suffix):
        raise ValueError(f"{who}: suffix must be a single non-empty string")
    if gpa and align != "procrustes":
        raise ValueError(f'{who}: gpa rotates the planes, so it needs align="procrustes"')
    net = _copy(net)
    planes = current_layouts(net)
    if len(planes) < 2:
        raise ValueError(f"{who}: at least two planes are needed to compare them, and there "
                         + ("is 1" if len(planes) == 1 else f"are {len(planes)}"))
    names = list(planes)
    if which is None:
        selected = names
    elif all(isinstance(w, (int, np.integer)) for w in np.atleast_1d(which)):
        selected = [names[i - 1] for i in np.atleast_1d(which)]
    else:
        missing = [w for w in np.atleast_1d(which) if w not in names]
        if missing:
            raise ValueError(f'{who}: there is no plane called "{missing[0]}". Available: ' + ", ".join(names))
        selected = list(np.atleast_1d(which))
    ref = plane_name(to, names, who)
    use = list(dict.fromkeys(selected + [ref]))
    L = {nm: as_plane(planes[nm], nm, who) for nm in use}
    rows = {z.shape[0] for z in L.values()}
    if len(rows) > 1:
        raise ValueError(f"{who}: the planes hold a different number of rows ({sorted(rows)}); every one of them "
                         "must hold one per node")
    cols = {z.shape[1] for z in L.values()}
    if len(cols) > 1:
        raise ValueError(f"{who}: the planes hold a different number of columns ({sorted(cols)}); they can only be "
                         "laid over one another with the same ones")
    n, p = rows.pop(), cols.pop()
    w_all, w_source = layout_weights(net, weights, n, who)
    ok = np.all([~np.isnan(z).any(axis=1) for z in L.values()], axis=0)
    if w_all is not None:
        ok &= ~np.isnan(w_all)
    if not ok.any():
        raise ValueError(f"{who}: no row is complete in every plane")
    w = None if w_all is None else w_all[ok]

    centre, scale, rot = {}, {}, {}
    for nm in use:
        z = L[nm][ok]
        if align == "canvas":
            d = z.max(axis=0) - z.min(axis=0)
            if not np.all(np.isfinite(d)) or np.any(d <= 0):
                raise ValueError(f'{who}: one axis of "{nm}" has no range, so it cannot be brought to the canvas')
            centre[nm], scale[nm] = z.min(axis=0), d
        elif align != "none":
            _, centre[nm], scale[nm] = norm_plane(z, w, who)
        else:
            centre[nm], scale[nm] = np.zeros(p), np.ones(p)
        rot[nm] = np.eye(p)

    def common(nm, z=None):
        z = L[nm] if z is None else z
        return (z - centre[nm]) / scale[nm]

    consensus = None
    if align == "procrustes":
        Z = {nm: common(nm)[ok] for nm in use}
        if not gpa:
            for nm in use:
                if nm != ref:
                    rot[nm] = procrustes_rotate(Z[nm], Z[ref], w)[1]
        else:
            M, previous, converged = Z[ref].copy(), np.inf, False
            for _ in range(maxit):
                for nm in selected:
                    rot[nm] = procrustes_rotate(Z[nm], M, w)[1]
                Y = [Z[nm] @ rot[nm] for nm in selected]
                M = sum(Y) / len(Y)
                ss = sum(((y - M) ** 2).sum(axis=1).dot(np.ones(len(y)) if w is None else w) for y in Y)
                if previous - ss <= tol * max(1, ss):
                    converged = True
                    break
                previous = ss
            if not converged:
                warnings.warn(f"{who}: the generalized Procrustes fit did not settle in {maxit} iterations; raise maxit")
            q = procrustes_rotate(M, Z[ref], w)[1]
            for nm in selected:
                rot[nm] = rot[nm] @ q
            consensus = np.full((n, p), np.nan)
            consensus[ok] = M @ q

    def to_ref(u):
        return u * scale[ref] + centre[ref]

    aligned = {}
    for nm in selected:
        a = to_ref(common(nm) @ rot[nm])
        frame = pd.DataFrame(a, columns=rotated_axis_titles(list(planes[nm].columns), rot[nm]) or list(planes[nm].columns))
        aligned[nm] = frame
    lay = dict(planes)
    if suffix is None:
        lay.update(aligned)
        net.nodes["fx"] = lay[next(iter(lay))].iloc[:, 0].to_numpy()
        net.nodes["fy"] = lay[next(iter(lay))].iloc[:, 1].to_numpy()
    else:
        add = aligned if gpa else {k: v for k, v in aligned.items() if k != ref}
        for nm, frame in add.items():
            frame.columns = rotated_axis_titles(list(frame.columns), np.eye(frame.shape[1]), suffix)
        asked = [nm + suffix for nm in add]
        given, taken = [], list(lay)
        for a in asked:
            g = make_unique(taken, a)
            given.append(g)
            taken.append(g)
        if any(g != a for g, a in zip(given, asked)):
            clash = [(a, g) for a, g in zip(asked, given) if a != g]
            warnings.warn(f"{who}: the object already held " + name_clash([c[0] for c in clash], [c[1] for c in clash])
                          + ". Choose another suffix to tell them apart.")
        name_of = dict(zip(add, given))
        merged = {}
        for nm in lay:
            merged[nm] = lay[nm]
            if nm in add:
                merged[name_of[nm]] = add[nm]
        lay = merged
    net.layouts = lay
    if len(lay) > 1 and net.options.get("axesLabels") is None:
        planes_axes_labels(net)
    if net.scat.get("cluster_columns"):
        net.clusters = current_clusters(net)
    net.scat["alignment"] = {"align": align, "to": ref, "gpa": gpa, "rows": int(ok.sum()), "dropped": int((~ok).sum()),
                             "center": {k: centre[k] for k in selected}, "scale": {k: scale[k] for k in selected},
                             "rotation": {k: rot[k] for k in selected}, "weights": w_all,
                             **({"consensus": consensus} if consensus is not None else {})}
    return net


# ---------------------------------------------------------------------- add_axes ----
def planes_axes_labels(net):
    """The titles of the axes of every plane, as the viewer takes them (``options['layoutsAxesLabels']``)."""
    if net.options.get("axesLabels") is not None or not net.layouts or len(net.layouts) < 2:
        net.options.pop("layoutsAxesLabels", None)
        return net
    net.options["layoutsAxesLabels"] = {k: ([str(c) for c in m.columns] if len(m.columns) == 2 else ["", ""])
                                        for k, m in net.layouts.items()}
    return net


def split_axis_pairs(axes) -> dict:
    """A table of an even number of columns holds one plane per consecutive pair, named ``"c1-c2"``."""
    if isinstance(axes, pd.DataFrame):
        frame = axes
    else:
        frame = pd.DataFrame(np.asarray(axes))
        if isinstance(axes, np.ndarray) or not hasattr(axes, "columns"):
            frame.columns = [f"axis{i + 1}" for i in range(frame.shape[1])]
    if not all(pd.api.types.is_numeric_dtype(frame[c]) for c in frame.columns):
        raise ValueError("axes must hold numeric coordinates")
    if frame.shape[1] < 2 or frame.shape[1] % 2:
        raise ValueError("axes must have an even number of columns, a pair of coordinates per plane, but it has "
                         f"{frame.shape[1]}")
    out = {}
    for j in range(0, frame.shape[1], 2):
        pair = frame.iloc[:, [j, j + 1]].reset_index(drop=True)
        pair.columns = [str(c) for c in frame.columns[[j, j + 1]]]
        out["-".join(pair.columns)] = pair
    return out


def select_l4c_reductions(obj, which=None) -> dict:
    entries = obj.reductions
    if not entries:
        raise ValueError("looking4clusters object holds no dimensionality reduction")
    if which is None:
        keys = list(entries)
    elif isinstance(which, (int, np.integer)) or (isinstance(which, (list, tuple)) and all(isinstance(w, int) for w in which)):
        idx = list(np.atleast_1d(which))
        if any(i < 1 or i > len(entries) for i in idx):
            raise ValueError(f"which must be indices between 1 and {len(entries)}")
        keys = [list(entries)[i - 1] for i in idx]
    else:
        wanted = [which] if isinstance(which, str) else list(which)
        unknown = [w for w in wanted if w not in entries]
        if unknown:
            raise ValueError(f'which="{unknown[0]}" is not among the reductions of the looking4clusters object: '
                             + ", ".join(entries))
        keys = wanted
    return {k: pd.DataFrame(entries[k]).reset_index(drop=True) for k in keys}


def extract_axes_list(axes, which=None) -> dict:
    """Any accepted input as a dict name -> two-column DataFrame (an empty name: the caller supplies it)."""
    from .clusters import Looking4Clusters

    if isinstance(axes, Looking4Clusters):
        return select_l4c_reductions(axes, which)
    if isinstance(axes, dict):
        if not axes:
            raise ValueError("axes holds no coordinates")
        out = {}
        for own, value in axes.items():
            pairs = split_axis_pairs(value)
            if own:
                pairs = {own: next(iter(pairs.values()))} if len(pairs) == 1 else {f"{own}.{k}": v for k, v in pairs.items()}
            out.update(pairs)
        return out
    if isinstance(axes, (list, tuple)) and axes and all(isinstance(a, (pd.DataFrame, np.ndarray)) for a in axes):
        out = {}
        for a in axes:
            out.update(split_axis_pairs(a))
        return out
    return split_axis_pairs(axes)


def align_axes(net, pair: pd.DataFrame, weight=None, label: str = "") -> pd.DataFrame:
    """Bring a pair of coordinates to the rows of the node table (cases are collapsed into patterns)."""
    what = f'axes "{label}"' if label else "axes"
    n = len(net.nodes)
    if len(pair) == n:
        return pair
    idx = net.scat.get("case_to_pattern")
    if idx is None or len(pair) != len(idx):
        raise ValueError(f"{what} has {len(pair)} rows but scatObj has {n} nodes. Coordinates must hold one row per "
                         "node, or one per case when the nodes are patterns collapsed by v_patterns." + sampled_note(net))
    if weight is None:
        weight = net.scat.get("case_weight")
    if weight is not None and len(weight) != len(pair):
        raise ValueError(f"weight must have one value per row of {what}")
    values = pair.to_numpy(dtype=float)
    idx = np.asarray(idx)
    out = np.empty((n, values.shape[1]))
    for g in range(1, n + 1):
        sel = idx == g
        out[g - 1] = values[sel].mean(axis=0) if weight is None else (values[sel] * np.asarray(weight)[sel][:, None]).sum(axis=0) / np.asarray(weight)[sel].sum()
    return pd.DataFrame(out, columns=pair.columns)


def add_axes(net, axes, name: str | None = None, which=None, weight=None, axes_labels=None):
    """Add one or several planes to a scattergram (``addAxes``).

    ``axes`` is a table with an even number of columns (a pair per plane), a list of them, a dict or a
    :class:`Looking4Clusters`, whose reductions become planes. ``name`` names the plane (a prefix if
    several); ``weight`` weighs the cases when they are collapsed into patterns. The coordinates already
    drawn become the first plane. ``axes_labels``: ``None`` (the default) lets every plane keep its own
    axis titles, a pair puts the same titles on every plane, and ``"keep"`` leaves what the object had.

    Parameters
    ----------
    net : netCoin
        A scattergram made by :func:`sur_scat`.
    axes : DataFrame, list, dict or Looking4Clusters
        The coordinates, one pair of columns per plane.
    name : str, optional
        Name of the plane (a prefix when there are several).
    which : optional
        For a :class:`Looking4Clusters`, which reductions to add.
    weight : sequence, optional
        Weight of each case, used when the cases are collapsed into patterns.
    axes_labels : None, "keep" or pair of str
        Titles of the axes; see above.

    Returns
    -------
    netCoin
        A modified copy; the planes are in ``net.layouts``.
    """
    if not isinstance(net, netCoin):
        raise TypeError("net must be a netCoin object returned by sur_scat")
    if axes is None:
        raise ValueError("axes must be provided")
    net = _copy(net)
    sets = extract_axes_list(axes, which)
    several = len(sets) > 1
    layouts = dict(current_layouts(net))
    renamed, renamed_to = [], []
    for own, plane in sets.items():
        own = own or ""
        key = (own or "axes") if name is None else (f"{name}.{own}" if several else name)
        asked = key
        key = make_unique(list(layouts), key)
        if key != asked:
            renamed.append(asked)
            renamed_to.append(key)
        layouts[key] = align_axes(net, plane, weight, key)
    if renamed:
        warnings.warn("the object already held " + name_clash(renamed, renamed_to) + ". State name to tell them apart, "
                      "or drop the previous planes to avoid repeating a set of axes.")
    net.layouts = layouts
    net.clusters = current_clusters(net)
    if not (isinstance(axes_labels, str) and axes_labels == "keep") and len(layouts) > 1:
        net.options.pop("axesLabels", None)
        if axes_labels is not None:
            net.options["axesLabels"] = [str(a) for a in axes_labels]
        planes_axes_labels(net)
    elif len(layouts) > 1 and net.options.get("axesLabels") is None:
        planes_axes_labels(net)
    return net


# ------------------------------------------------------------------- clusters ----
def _codes(values) -> np.ndarray:
    """``as.integer(as.factor(x))``: the groups numbered by their sorted values."""
    from .clusters import as_factor

    return np.asarray(as_factor(values).codes, dtype=int) + 1


def select_l4c_clusters(obj, which=None) -> np.ndarray:
    entries = obj.clusters
    if not entries:
        raise ValueError("looking4clusters object has no clusters")
    k = None
    if which is None:
        key = obj.options.get("myGroups")
        if key is None or key not in entries:
            key = next(iter(entries))
    elif isinstance(which, (int, np.integer)):
        if which < 1 or which > len(entries):
            raise ValueError(f"which must be a single index between 1 and {len(entries)}")
        key = list(entries)[which - 1]
    else:
        parts = str(which).split(":")
        key = parts[0]
        if len(parts) > 1:
            try:
                k = int(parts[1])
            except ValueError:
                raise ValueError(f'the number of groups in which="{which}" must be an integer')
        if key not in entries:
            raise ValueError(f'which="{key}" is not among the clusterizations of the looking4clusters object: '
                             + ", ".join(entries))
    entry = entries[key]
    if not isinstance(entry, pd.DataFrame):
        if k is not None:
            raise ValueError(f'clusterization "{key}" holds a single solution, so a number of groups cannot be chosen')
        return _codes(entry)
    if k is None:
        k = obj.optim.get(key)
        if k is None:
            return _codes(entry.iloc[:, 0])
    column = f"levels_{k}"
    if column not in entry.columns:
        raise ValueError(f'clusterization "{key}" has no solution with {k} groups. Available: '
                         + ", ".join(c.replace("levels_", "") for c in entry.columns))
    return _codes(entry[column])


def list_l4c_clusters(obj) -> dict:
    """Every clusterization of a looking4clusters object as ``{"method(k)": codes}``."""
    entries = obj.clusters
    if not entries:
        raise ValueError("looking4clusters object has no clusters")
    out = {}
    for key, entry in entries.items():
        if not isinstance(entry, pd.DataFrame):
            codes = _codes(entry)
            out[f"{key}({len(set(codes))})"] = codes
        else:
            for column in entry.columns:
                out[f"{key}({column.replace('levels_', '')})"] = _codes(entry[column])
    return out


def extract_cluster_list(x, which=None) -> dict:
    from .clusters import Looking4Clusters

    if isinstance(x, Looking4Clusters) and which is None:
        return list_l4c_clusters(x)
    if isinstance(x, Looking4Clusters):
        return {"": select_l4c_clusters(x, which)}
    if isinstance(x, pd.DataFrame):
        return {str(c): _codes(x[c].astype(str)) for c in x.columns}
    if isinstance(x, dict):
        return {str(k): _codes(np.asarray(v).astype(str)) for k, v in x.items()}
    if hasattr(x, "labels_"):                                          # a fitted scikit-learn clusterer
        return {"": np.asarray(x.labels_, dtype=int) + 1}
    return {"": _codes(np.asarray(x).astype(str))}


def add_clusters(net, clusters, name: str | None = None, sort: bool = True, weight=None, which=None,
                 max_groups: int | None = 25, rest: int | None = None):
    """Add one or several cluster columns to a scattergram (``addClusters``).

    ``clusters`` is a vector (one group per node, or per case when the nodes are patterns), a table of
    one clusterization per column, a fitted scikit-learn clusterer or a :class:`Looking4Clusters`
    (every clusterization it holds, or the one ``which`` chooses). The groups are renumbered from left
    to right on the first axis when ``sort``. ``max_groups`` skips clusterizations with more groups
    (a continuous variable passed by mistake); ``rest`` is the number of a residual group.

    Parameters
    ----------
    net : netCoin
        A scattergram made by :func:`sur_scat`.
    clusters : sequence, DataFrame, clusterer or Looking4Clusters
        The groups.
    name : str, optional
        Name of the column (a prefix when there are several).
    sort : bool
        Renumber the groups from left to right on the first axis.
    weight : sequence, optional
        Weight of each node or case, for the centroids and the order of the groups.
    which : optional
        For a :class:`Looking4Clusters`, which clusterizations to add.
    max_groups : int, optional
        Skip clusterizations with more groups than this.
    rest : int, optional
        Number of the group that gathers what is left over.

    Returns
    -------
    netCoin
        A modified copy; the columns are added to ``net.nodes``, and the centroids of every group on every plane
        to ``net.clusters``.
    """
    if not isinstance(net, netCoin):
        raise TypeError("net must be a netCoin object returned by sur_scat")
    if clusters is None:
        raise ValueError("clusters must be provided")
    net = _copy(net)
    sets = extract_cluster_list(clusters, which)
    several = len(sets) > 1
    renamed, renamed_to = [], []
    for own, values in sets.items():
        n_groups = len(set(values))
        if max_groups is not None and n_groups > max_groups:
            warnings.warn("clusterization " + (f'"{own}" ' if own else "") + f"has {n_groups} groups, more than "
                          f"max_groups={max_groups}, and was not added. Raise max_groups to add it anyway.")
            continue
        column = None
        if name is None:
            column = own or None
        else:
            column = f"{name}.{own}" if several else name
        if column is not None:
            asked = column
            column = make_unique(list(net.nodes.columns), column)
            if column != asked:
                renamed.append(asked)
                renamed_to.append(column)
        _add_one_cluster(net, values, column, sort, weight, rest)
    if renamed:
        warnings.warn("the node table already held " + name_clash(renamed, renamed_to) + ". State name to tell them "
                      "apart, or drop the previous columns to avoid repeating a clusterization.")
    net.clusters = current_clusters(net)
    net.cluster_means = current_cluster_means(net)
    return net


def _add_one_cluster(net, clusters, name, sort, weight, rest):
    language = net.options.get("language") or "en"
    clusters = np.asarray(clusters)
    case_clusters = case_idx = display = None
    n_nodes = len(net.nodes)
    if len(clusters) != n_nodes:
        idx = net.scat.get("case_to_pattern")
        if idx is not None and len(clusters) == len(idx):
            case_clusters = clusters.astype(str)
            case_idx = np.asarray(idx)
            collapsed, display = [], []
            for g in range(1, n_nodes + 1):
                values = case_clusters[case_idx == g]
                counts = pd.Series(values).value_counts()
                ordered = sorted(counts.index, key=lambda s: (-counts[s], _collation_key(s)))
                collapsed.append(ordered[0])
                display.append("|".join(ordered))
            clusters = np.array(collapsed)
        else:
            raise ValueError(f"clusters has length {len(clusters)} but scatObj has {n_nodes} nodes. Clusters must hold "
                             "one assignment per node, or one per case when the nodes are patterns collapsed by "
                             "v_patterns." + sampled_note(net))
    if name is None:
        name = make_unique(list(net.nodes.columns), "clusters")
    cl = np.array([str(c) for c in clusters])
    unique_cl = list(dict.fromkeys(cl))
    all_values = list(unique_cl)
    if display is not None:
        shown = list(dict.fromkeys(v for d in display for v in d.split("|")))
        all_values = list(dict.fromkeys(all_values + shown))
    rest_num = None if rest is None else int(rest)
    first = net.nodes["fx"].to_numpy(dtype=float) if "fx" in net.nodes.columns else None
    sorted_values = None
    if sort and first is not None and len(first) == len(cl):
        means = {}
        for u in unique_cl:
            sel = cl == u
            means[u] = np.average(first[sel], weights=None if weight is None else np.asarray(weight)[sel])
        mu = {v: means.get(v, np.nan) for v in all_values}
        if rest_num is not None and str(rest_num) in mu:
            mu[str(rest_num)] = np.inf
        order = sorted(range(len(all_values)), key=lambda i: (np.isnan(mu[all_values[i]]), mu[all_values[i]]))
        sorted_values = [all_values[i] for i in order]
        renum = {v: str(i + 1) for i, v in enumerate(sorted_values)}
        cl = np.array([renum[c] for c in cl])
        if display is not None:
            display = ["|".join(renum.get(v, "") for v in d.split("|")) for d in display]
        unique_cl = [str(i) for i in range(1, len(sorted_values) + 1)]
        if rest_num is not None:
            rest_num = len(sorted_values)
    values = display if display is not None else list(cl)
    unique_values = list(dict.fromkeys(values))

    group_word = re.sub(r"\([0-9]+\)(\.[0-9]+)?$", "", name)
    group_word = (group_word[:1].upper() + group_word[1:]) if group_word else getByLanguage(groupList, language)

    def split(val):
        return [int(x) for x in val.split("|")]

    width = max(len(str(g)) for v in unique_values for g in split(v))

    def padded(val):
        return "|".join(f"{g:0{width}d}" for g in split(val))

    alone = [v for v in unique_values if "|" not in v]
    singles = list(dict.fromkeys(alone + [str(g) for v in unique_values for g in split(v)]))
    by_number = lambda vals: sorted(vals, key=padded)
    ordered_values = by_number(singles) + by_number([v for v in unique_values if v not in singles])
    rest_word = getByLanguage(restList, language)

    def shown_label(val):
        return [rest_word if rest_num is not None and g == rest_num else f"{g:0{width}d}" for g in split(val)]

    labels = ["|".join(f"{group_word}: {s}" for s in shown_label(v)) for v in ordered_values]
    mapping = dict(zip(ordered_values, labels))
    net.nodes[name] = pd.Categorical([mapping[v] for v in values], categories=labels, ordered=True)

    if case_clusters is not None:
        numbers = ([str(sorted_values.index(c) + 1) for c in case_clusters] if sorted_values is not None
                   else list(case_clusters))
        nums = sorted(set(numbers), key=int)
        case_label = pd.Categorical([f"{group_word}: " + "|".join(shown_label(n)) for n in numbers],
                                    categories=[f"{group_word}: " + "|".join(shown_label(n)) for n in nums])
        cw = net.scat.get("case_weight")
        if cw is not None and len(cw) != len(case_idx):
            cw = None
        net.scat.setdefault("cluster_shares", {})[name] = cluster_shares(case_idx, n_nodes, case_label, cw)
        net.scat[f"collapse_{name}"] = {"n_patterns": n_nodes, "n_cases": len(case_clusters)}
    columns = list(net.scat.get("cluster_columns", []))
    if name not in columns:
        columns.append(name)
    net.scat["cluster_columns"] = columns


def replace_clusters(net, clusters, name: str | None = None, sort: bool = True, weight=None, which=None,
                     max_groups: int | None = 25):
    """Replace the cluster columns of a scattergram by new ones (``replaceClusters``).

    Parameters
    ----------
    net : netCoin
        A scattergram made by :func:`sur_scat`.
    clusters, name, sort, weight, which, max_groups
        As in :func:`add_clusters`; the cluster columns the object held are removed first.

    Returns
    -------
    netCoin
        A modified copy.
    """
    from .clusters import Looking4Clusters

    net = _copy(net)
    old = [c for c in net.scat.get("cluster_columns", []) if c in net.nodes.columns]
    if not old:
        old = [c for c in net.nodes.columns if re.match(r"^Groups?\([0-9]", c)]
    if old:
        net.nodes = net.nodes.drop(columns=old)
        net.scat["cluster_columns"] = [c for c in net.scat.get("cluster_columns", []) if c not in old]
        shares = {k: v for k, v in net.scat.get("cluster_shares", {}).items() if k not in old}
        net.scat["cluster_shares"] = shares
    several = (isinstance(clusters, Looking4Clusters) and which is None) or (isinstance(clusters, pd.DataFrame) and clusters.shape[1] > 1)
    if name is None and old and not several:
        name = old[0]
    return add_clusters(net, clusters, name=name, sort=sort, weight=weight, which=which, max_groups=max_groups)
