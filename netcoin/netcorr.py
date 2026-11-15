"""Networks of correlations: ``net_corr`` and ``d_net_corr`` (ports of ``netCorr`` and ``d_netCorr``)."""
from __future__ import annotations

import math
import warnings
from functools import lru_cache
from itertools import permutations
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .core import _normalize_legacy, name_by_language, netCoin
from .layouts import layout_fact, layout_pca, layout_special
from .utils import frequencyNames, getByLanguage

_KNOWN = ["pearson", "spearman", "kendall", "partial", "dcor"]


# ------------------------------------------------------------ the arguments ----
def _corr_methods(method) -> list[str]:
    """The methods asked for, by their full names (``"p"`` is Pearson; other names can be shortened)."""
    methods = [method] if isinstance(method, str) else list(method)
    if not methods or any(m is None for m in methods):
        raise ValueError("method must name at least one of: " + ", ".join(_KNOWN))
    out = []
    for m in methods:
        m = str(m).lower()
        if m == "p":
            full = "pearson"
        else:
            hits = [k for k in _KNOWN if k == m] or [k for k in _KNOWN if k.startswith(m)]
            if len(hits) != 1:
                raise ValueError(f"method '{m}' is not available. Choose among: " + ", ".join(_KNOWN))
            full = hits[0]
        if full not in out:
            out.append(full)
    return out


def _corr_criteria(criteria, methods: list[str]) -> str:
    """The column of the links a criterion refers to: ``value``, ``z``, ``p``, a method or ``z(method)``/``p(method)``."""
    cr = str(criteria if isinstance(criteria, str) else list(criteria)[0]).lower()
    cols = ["value", "z", "p"]
    if len(methods) > 1:
        cols += methods[1:] + [f"z({m})" for m in methods[1:]] + [f"p({m})" for m in methods[1:]]
    if cr == methods[0]:
        return "value"
    if cr[:2] in ("z(", "p(") and cr.endswith(")") and cr[2:-1] == methods[0]:
        return cr[0]
    if cr not in cols:
        raise ValueError(f"criteria '{cr}' is not available. Choose among: "
                         + ", ".join(dict.fromkeys([methods[0]] + cols)))
    return cr


# ------------------------------------------------------ coefficients of association ----
def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    with np.errstate(divide="ignore", invalid="ignore"):
        x, y = x - x.mean(), y - y.mean()
        return float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum()))


def _pair_matrix(X: np.ndarray, w, f) -> np.ndarray:
    """A function of two variables over the cases observed in both, for every pair."""
    m = X.shape[1]
    seen = ~np.isnan(X)
    out = np.eye(m)
    for j in range(m - 1):
        for i in range(j + 1, m):
            idx = seen[:, i] & seen[:, j]
            value = np.nan if idx.sum() < 2 else f(X[idx, i], X[idx, j], None if w is None else w[idx])
            out[i, j] = out[j, i] = value
    return out


def _wpearson(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted Pearson correlations, each pair resting on the cases observed in both variables."""
    seen = ~np.isnan(X)
    present = seen.astype(float)
    X = np.where(seen, X, 0.0)
    wc = w[:, None]
    mu = (X * wc).sum(axis=0) / (present * wc).sum(axis=0)
    X = (X - mu) * present
    W = (present * wc).T @ present
    A = (X * wc).T @ present
    B = (X ** 2 * wc).T @ present
    C = (X * wc).T @ X
    with np.errstate(divide="ignore", invalid="ignore"):
        V = B / W - (A / W) ** 2
        R = (C / W - A * A.T / W ** 2) / np.sqrt(V * V.T)
    np.fill_diagonal(R, 1.0)
    return R


def _wcor2(x, y, w) -> float:
    total = w.sum()
    x = x - (w * x).sum() / total
    y = y - (w * y).sum() / total
    return float((w * x * y).sum() / np.sqrt((w * x * x).sum() * (w * y * y).sum()))


def _wrank(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Ranks weighted by ``w``; the tied values share the middle of their block."""
    order = np.argsort(x, kind="stable")
    xs = x[order]
    group = np.concatenate([[0], np.cumsum(xs[1:] != xs[:-1])])
    mass = np.bincount(group, weights=w[order])
    out = np.empty(len(x))
    out[order] = (np.cumsum(mass) - mass / 2)[group]
    return out


def _ucenter(x: np.ndarray) -> np.ndarray:
    """Distances centred the unbiased way (Szekely and Rizzo 2014), as a vector."""
    n = len(x)
    a = np.abs(x[:, None] - x[None, :])
    rs = a.sum(axis=1)
    a = a - (rs[:, None] + rs[None, :]) / (n - 2) + a.sum() / ((n - 1) * (n - 2))
    np.fill_diagonal(a, 0.0)
    return a.ravel()


def _dcor_matrix(X: np.ndarray) -> np.ndarray:
    n, m = X.shape
    if n < 4:
        raise ValueError("distance correlation needs at least four cases")
    if float(n) ** 2 * m > 6e7:
        raise ValueError("too many cases for distance correlation (it keeps n^2 distances per variable)")
    U = np.column_stack([_ucenter(X[:, j]) for j in range(m)])
    G = U.T @ U / (n * (n - 3))
    d = np.diag(G)
    R = G / np.sqrt(np.outer(d, d))
    np.fill_diagonal(R, 1.0)
    return R


def _dcor2(x, y, w=None) -> float:
    if len(x) < 4:
        return np.nan
    ux, uy = _ucenter(x), _ucenter(y)
    return float((ux * uy).sum() / np.sqrt((ux * ux).sum() * (uy * uy).sum()))


def _rank_columns(X: np.ndarray) -> np.ndarray:
    return np.column_stack([stats.rankdata(X[:, j]) for j in range(X.shape[1])])


def _correlations(X: np.ndarray, method: str, pairwise: bool) -> np.ndarray:
    """Pearson, Spearman or Kendall (tau-b) between the columns, without weights."""
    m = X.shape[1]
    if method == "kendall":
        def tau(x, y, w=None):
            return float(stats.kendalltau(x, y).statistic)
        return _pair_matrix(X, None, tau)
    if method == "spearman":
        if pairwise:
            return _pair_matrix(X, None, lambda x, y, w: _pearson(stats.rankdata(x), stats.rankdata(y)))
        X = _rank_columns(X)
    with np.errstate(divide="ignore", invalid="ignore"):
        if pairwise:
            return _pair_matrix(X, None, lambda x, y, w: _pearson(x, y))
        return np.atleast_2d(np.corrcoef(X, rowvar=False)) if m > 1 else np.ones((1, 1))


# ----------------------------------------------------------- exact contrasts ----
@lru_cache(maxsize=None)
def _kendall_cdf(n: int) -> np.ndarray:
    """P(inversions <= k) of a permutation of ``n`` (Mahonian numbers), the null of Kendall's S."""
    f = np.array([1.0])
    for k in range(2, n + 1):
        g = np.zeros(k * (k - 1) // 2 + 1)
        for j in range(k):
            g[j:j + len(f)] += f
        f = g
    return np.cumsum(f) / f.sum()


@lru_cache(maxsize=None)
def _spearman_null(n: int):
    """Distribution of the sum of squared rank differences D of ``n`` cases without ties."""
    perms = np.array(list(permutations(range(n))), dtype=np.int16)
    d = ((perms - np.arange(n)) ** 2).sum(axis=1)
    values, counts = np.unique(d, return_counts=True)
    return values.astype(float), np.cumsum(counts) / counts.sum()


def _spearman_exact_p(D: np.ndarray, n: int, bipolar: bool) -> np.ndarray:
    values, cum = _spearman_null(n)
    cum = np.concatenate([[0.0], cum])
    le = cum[np.searchsorted(values, D, side="right")]          # P(D0 <= D)
    if not bipolar:
        return le
    lt = cum[np.searchsorted(values, D, side="left")]           # P(D0 < D)
    return np.minimum(1.0, 2 * np.minimum(le, 1 - lt))


def _spearman_d(X: np.ndarray, seen: np.ndarray, ij, pairwise: bool) -> np.ndarray:
    """D of the pairs for samples of at most nine cases; NaN where there are ties."""
    rows, cols = ij
    if not pairwise:
        tied = np.array([len(np.unique(X[:, j])) < X.shape[0] for j in range(X.shape[1])])
        ranks = _rank_columns(X)
        G = ranks.T @ ranks
        D = np.diag(G)[rows] + np.diag(G)[cols] - 2 * G[rows, cols]
        D = np.where(tied[rows] | tied[cols], np.nan, D)
        return np.round(D)
    out = np.empty(len(rows))
    for k, (i, j) in enumerate(zip(rows, cols)):
        idx = seen[:, i] & seen[:, j]
        x, y = X[idx, i], X[idx, j]
        if len(np.unique(x)) < len(x) or len(np.unique(y)) < len(y):
            out[k] = np.nan
        else:
            out[k] = ((stats.rankdata(x) - stats.rankdata(y)) ** 2).sum()
    return out


def _tie_sums(x: np.ndarray) -> np.ndarray:
    _, t = np.unique(x, return_counts=True)
    t = t.astype(float)
    return np.array([(t * (t - 1)).sum(), (t * (t - 1) * (2 * t + 5)).sum(), (t * (t - 1) * (t - 2)).sum()])


def _kendall_ties(X, seen, ij, pairwise):
    rows, cols = ij
    if not pairwise:
        C = np.column_stack([_tie_sums(X[:, j]) for j in range(X.shape[1])])
        return C[:, cols], C[:, rows]
    x, y = [], []
    for i, j in zip(rows, cols):
        idx = seen[:, i] & seen[:, j]
        x.append(_tie_sums(X[idx, j]))
        y.append(_tie_sums(X[idx, i]))
    return np.array(x).T.reshape(3, -1), np.array(y).T.reshape(3, -1)


# ------------------------------------------------------------ the statistics ----
def corr_stat(X: np.ndarray, w: np.ndarray | None, method: str, pairwise: bool, bipolar: bool) -> dict:
    """Coefficients and contrasts of one method (R's ``corrStat``).

    Returns ``R`` (matrix of coefficients) and, for the pairs of the lower triangle taken
    column by column, the statistic ``z``, its p-value ``p`` and the number of cases ``n`` it
    rests on (Kish's effective number with weights). The contrast is one-sided (positive
    association) unless the links are bipolar.
    """
    m = X.shape[1]
    seen = ~np.isnan(X)
    present = seen.astype(float)
    cols, rows = np.triu_indices(m, 1)           # lower triangle, column by column
    ij = (rows, cols)
    if w is None:
        Ne = present.T @ present
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            Ne = ((present * w[:, None]).T @ present) ** 2 / ((present * (w ** 2)[:, None]).T @ present)
    df = None

    if method == "pearson":
        R = _correlations(X, "pearson", pairwise) if w is None else _wpearson(X, w)
    elif method == "spearman":
        if w is None:
            R = _correlations(X, "spearman", pairwise)
        elif not pairwise:
            R = _wpearson(np.column_stack([_wrank(X[:, j], w) for j in range(m)]), w)
        else:
            R = _pair_matrix(X, w, lambda x, y, ww: _wcor2(_wrank(x, ww), _wrank(y, ww), ww))
    elif method == "kendall":
        R = _correlations(X, "kendall", pairwise)
    elif method == "partial":
        complete = seen.all(axis=1)
        if not complete.all():
            warnings.warn("partial correlations rest on the cases complete in all the variables")
        Xc = X[complete]
        wc = None if w is None else w[complete]
        Rc = np.corrcoef(Xc, rowvar=False) if wc is None else _wpearson(Xc, wc)
        nn = len(Xc) if wc is None else wc.sum() ** 2 / (wc ** 2).sum()
        if nn - m <= 0:
            raise ValueError("partial correlations need more complete cases than variables")
        try:
            theta = np.linalg.inv(Rc)
        except np.linalg.LinAlgError:
            raise ValueError("partial correlations cannot be computed: the correlation matrix is singular")
        R = -theta / np.sqrt(np.outer(np.diag(theta), np.diag(theta)))
        np.fill_diagonal(R, 1.0)
        Ne = np.full((m, m), float(nn))
        df = np.full((m, m), float(nn - m))
    elif method == "dcor":
        R = _dcor_matrix(X) if not pairwise else _pair_matrix(X, None, _dcor2)
    else:
        raise ValueError(f"method '{method}' is not available")

    R = np.clip(R, -1, 1)
    r, n = R[ij], Ne[ij]
    df = n - 2 if df is None else df[ij]

    def p_t(t, dof):
        return 2 * stats.t.cdf(-np.abs(t), dof) if bipolar else stats.t.sf(t, dof)

    def p_z(z):
        return 2 * stats.norm.cdf(-np.abs(z)) if bipolar else stats.norm.sf(z)

    with np.errstate(divide="ignore", invalid="ignore"):
        if method in ("pearson", "spearman", "partial"):
            z = r * np.sqrt(df / (1 - r ** 2))
            p = p_t(z, df)
            if method == "spearman" and w is None:
                small = np.where((n >= 3) & (n <= 9))[0]
                if len(small):
                    D = _spearman_d(X, seen, (rows[small], cols[small]), pairwise)
                    ok = ~np.isnan(D)
                    chosen = small[ok]
                    for nk in np.unique(n[chosen]):
                        take = n[chosen] == nk
                        p[chosen[take]] = _spearman_exact_p(D[ok][take], int(nk), bipolar)
        elif method == "kendall":
            tx, ty = _kendall_ties(X, seen, ij, pairwise)
            n0 = n * (n - 1) / 2
            sc = r * np.sqrt((n0 - tx[0] / 2) * (n0 - ty[0] / 2))
            v = (n * (n - 1) * (2 * n + 5) - tx[1] - ty[1]) / 18 + tx[0] * ty[0] / (2 * n * (n - 1)) \
                + tx[2] * ty[2] / (9 * n * (n - 1) * (n - 2))
            z = sc / np.sqrt(v)
            z[n < 3] = np.nan
            p = p_z(z)
            exact = np.where((n >= 3) & (n < 50) & (tx[0] == 0) & (ty[0] == 0) & ~np.isnan(sc))[0]
            for nk in np.unique(n[exact]):
                kk = exact[n[exact] == nk]
                cdf = _kendall_cdf(int(nk))
                s = np.round(np.abs(sc[kk]) if bipolar else sc[kk])
                index = np.clip((nk * (nk - 1) / 2 - s) / 2, 0, len(cdf) - 1).astype(int)
                pr = cdf[index]
                p[kk] = np.minimum(1.0, 2 * pr) if bipolar else pr
        else:                                      # dcor
            v = n * (n - 3) / 2
            z = np.sqrt(v - 1) * r / np.sqrt(1 - r ** 2)
            z[n < 4] = np.nan
            p = p_t(z, v - 1)
    return {"R": R, "z": z, "p": p, "n": n}


# ---------------------------------------------------------------------- net_corr ----
def net_corr(
    variables: pd.DataFrame,
    weight: str | Sequence[float] | None = None,
    pairwise: bool = False,
    sort: bool = False,
    decreasing: bool = True,
    frequency: bool | None = None,
    means: bool = True,
    method: str | Sequence[str] = "pearson",
    criteria: str = "p",
    bonferroni: bool = False,
    min_link: float = 0,
    max_link: float = np.inf,
    sort_links: str | None = None,
    decreasing_links: bool = True,
    igraph: bool = False,
    **netcoin_args,
):
    """A network of the correlations among quantitative variables (port of ``netCorr``).

    Every column is a node (with its mean, standard deviation, minimum and maximum) and every
    pair a link with the coefficient (``value``), the statistic of its contrast (``z``) and its
    p-value (``p``). The links are bipolar: they are coloured by their sign and the contrast is
    two-sided, unless you say ``link_bipolar=False``.

    Parameters
    ----------
    variables : DataFrame
        One numeric column per variable (at least two) and one row per case.
    weight : str or sequence, optional
        The weight of each case, or the name of the column that holds it (not for
        ``"kendall"`` and ``"dcor"``). The cases with a missing weight are left out.
    pairwise : bool
        Each pair uses the cases that have both values; otherwise only the complete cases are
        used. The number of cases of each link is then kept in the column ``n``.
    sort, decreasing : bool
        Order the nodes by number of valid cases and, within them, by mean.
    frequency : bool, optional
        Put the number of valid cases in the nodes (by default only with ``pairwise``).
    means : bool
        Put the mean, standard deviation, minimum and maximum in the nodes.
    method : str or list of str
        ``"pearson"``, ``"spearman"``, ``"kendall"``, ``"partial"`` or ``"dcor"`` (distance
        correlation, bias corrected); they can be abbreviated. With several, the first gives
        ``value``, ``z`` and ``p`` and the others their own columns (``spearman``,
        ``z(spearman)``, ``p(spearman)``).
    criteria : str
        What ``min_link`` and ``max_link`` refer to: ``"value"``, ``"z"``, ``"p"`` (default),
        another method, or ``"z(method)"``/``"p(method)"``. For a p-value ``max_link`` is 0.5
        (1 with bipolar links) unless you say otherwise; for a coefficient of a bipolar link its
        absolute value counts.
    bonferroni : bool
        Divide a p-value ``max_link`` by the number of pairs.
    min_link, max_link : float
        Range of ``criteria`` kept as links.
    sort_links, decreasing_links
        Column and direction to order the links by.
    igraph : bool
        Return an ``igraph.Graph``.
    **netcoin_args
        The rest goes to :class:`netCoin` (``nodes``, ``name``, ``language``, ``layout``,
        ``lwidth``...). ``layout`` can also be ``"pc"``: principal components of the
        correlations. By default ``lwidth``, ``lweight`` and ``lcolor`` are ``"value"``.

    Returns
    -------
    netCoin or igraph.Graph
    """
    args = _normalize_legacy(netcoin_args, strict=False)
    language = args.get("language", "en")
    language = language[0] if isinstance(language, (list, tuple)) else language
    nodes = args.get("nodes")
    if isinstance(nodes, dict):
        nodes = pd.DataFrame(nodes)
    name = args["name"] = name_by_language(args.get("name"), language, nodes)
    methods = _corr_methods(method)
    frequency = pairwise if frequency is None else frequency

    args.setdefault("link_bipolar", True)
    args.setdefault("lcolor", "value")
    args.setdefault("lwidth", "value")
    args.setdefault("lweight", "value")
    bipolar = args["link_bipolar"] is True

    variables = pd.DataFrame(variables).copy()
    w = None
    if weight is not None:
        wrong = "weight must be the name of a column of variables or a numeric vector with a weight per case"
        if isinstance(weight, str):
            if weight not in variables.columns:
                raise ValueError(wrong)
            w = variables[weight].to_numpy(dtype=float)
            variables = variables.drop(columns=weight)
        else:
            w = np.asarray(weight, dtype=float)
        if w.ndim != 1 or len(w) != len(variables):
            raise ValueError(wrong)
        if (w[~np.isnan(w)] < 0).any():
            raise ValueError("weights cannot be negative")
        if {"kendall", "dcor"} & set(methods):
            raise ValueError("weights are not available for the methods kendall and dcor")
    not_numeric = [c for c in variables.columns
                   if not pd.api.types.is_numeric_dtype(variables[c]) or pd.api.types.is_bool_dtype(variables[c])]
    if not_numeric:
        raise ValueError("variables must be quantitative: " + ", ".join(map(str, not_numeric)))

    keep = np.ones(len(variables), bool) if pairwise else variables.notna().all(axis=1).to_numpy()
    if w is not None:
        keep &= ~np.isnan(w)
        w = w[keep]
    variables = variables[keep]
    X = variables.to_numpy(dtype=float)
    m = X.shape[1]
    if m < 2:
        raise ValueError("at least two variables are needed")
    if X.shape[0] < 3:
        raise ValueError("at least three complete cases are needed")
    names = [str(c) for c in variables.columns]

    # nodes: weighted statistics when there are weights
    def node_stats(j):
        ok = ~np.isnan(X[:, j])
        x = X[ok, j]
        if w is None:
            mu, s = x.mean(), (x.std(ddof=1) if len(x) > 1 else np.nan)
        else:
            ww = w[ok]
            mu = (ww * x).sum() / ww.sum()
            s = np.sqrt((ww * (x - mu) ** 2).sum() / (ww.sum() - (ww ** 2).sum() / ww.sum()))
        return len(x), mu, s, (x.min() if len(x) else np.nan), (x.max() if len(x) else np.nan)

    S = np.array([node_stats(j) for j in range(m)], dtype=float)
    if sort:
        sign = -1 if decreasing else 1
        order = sorted(range(m), key=lambda i: (sign * S[i, 0], sign * S[i, 1]))
        X, S = X[:, order], S[order]
        names = [names[i] for i in order]
    statistics = pd.DataFrame({name: names})
    if frequency:
        statistics[getByLanguage(frequencyNames, language)] = S[:, 0]
    if means:
        statistics["mean"] = np.round(S[:, 1], 2)
        statistics["std"] = np.round(S[:, 2], 2)
        statistics["min"] = S[:, 3]
        statistics["max"] = S[:, 4]
    if nodes is not None:
        user = nodes.copy()
        user[name] = user[name].astype(str)
        merged = statistics.merge(user, on=name, how="left", sort=False)
        # R's merge(sort = FALSE) lists the variables found in nodes first, then the others
        found = merged[name].isin(user[name]).to_numpy()
        args["nodes"] = pd.concat([merged[found], merged[~found]], ignore_index=True)
    else:
        args["nodes"] = statistics

    # links: the statistics of every method over the pairs of the lower triangle
    pw = bool(pairwise and np.isnan(X).any())
    results = [corr_stat(X, w, mt, pw, bipolar) for mt in methods]
    cols, rows = np.triu_indices(m, 1)
    first = results[0]
    edges = {"Source": [names[j] for j in cols], "Target": [names[i] for i in rows],
             "value": first["R"][rows, cols], "z": first["z"], "p": first["p"]}
    if pw:
        edges["n"] = np.round(first["n"], 2)
    for method_name, result in list(zip(methods, results))[1:]:
        edges[method_name] = result["R"][rows, cols]
        edges[f"z({method_name})"] = result["z"]
        edges[f"p({method_name})"] = result["p"]
    links = pd.DataFrame(edges)
    R = pd.DataFrame(first["R"], index=names, columns=names)

    criteria = _corr_criteria(criteria, methods)
    is_p = criteria == "p" or criteria.startswith("p(")
    if is_p and max_link == np.inf:
        max_link = 0.5 * (1 + bipolar)
    if is_p and bonferroni:
        max_link = max_link / math.comb(m, 2)
    selected = links[criteria].to_numpy(dtype=float)
    if bipolar and not is_p:
        selected = np.abs(selected)
    with np.errstate(invalid="ignore"):
        links = links[~np.isnan(selected) & (selected >= min_link) & (selected <= max_link)]
    if sort_links is not None:
        key = links[sort_links].to_numpy(dtype=float) * (-1 if decreasing_links else 1)
        links = links.iloc[np.argsort(key, kind="stable")]
    links = links.reset_index(drop=True)
    args["links"] = links if len(links) else None

    if args.get("layout") is not None:
        node_names = list(args["nodes"][name])
        args["layout"] = layout_special(args["layout"], {"pc": lambda: layout_pca(R).reindex(node_names)})
    net = netCoin(**args)
    if igraph:
        from .igraph_io import to_igraph

        return to_igraph(net)
    return net


# -------------------------------------------------------------------- d_net_corr ----
def _component_labels(membership: Sequence[int]) -> list[str]:
    """``C-1``, ``C-2``...: the components by decreasing size, numbered as R's ``compon`` does."""
    ids = sorted(set(membership))
    sizes = {i: sum(1 for m in membership if m == i) for i in ids}
    by_size = sorted(ids, key=lambda i: -sizes[i])
    width = max(len(str(abs(int(i)))) for i in ids)
    label = {old: "C-" + f"{rank:>{width}d}" for rank, old in zip(ids, by_size)}
    return [label[m] for m in membership]


def d_net_corr(
    variables: pd.DataFrame,
    nodes: pd.DataFrame | None = None,
    weight: str | Sequence[float] | None = None,
    pairwise: bool = False,
    frequency: bool | None = None,
    means: bool = True,
    method: str | Sequence[str] = "pearson",
    criteria: str = "value",
    bonferroni: bool = False,
    min_link: float = 0,
    max_link: float = np.inf,
    sort_links: str | None = None,
    decreasing_links: bool = True,
    factorial: str | Sequence[str] = "null",
    components: bool = True,
    backcomponents: bool = False,
    sequence: Sequence[float] | None = None,
    text_filter: Sequence[float] = (1, 0.99),
    speed: float = 50,
    directory: str | None = None,
    **netcoin_args,
):
    """A network of correlations that evolves as the threshold of the correlation rises.

    Port of ``d_netCorr``. It builds :func:`net_corr` and then one frame for every value of
    ``sequence``: the network of the links with ``abs(value)`` above it, with the centrality of
    each node, the component it belongs to (``C-1`` is the largest) and the ``ratio`` of each
    link to the threshold. Only the frames where the number of components grows are kept
    (every frame if ``components=False``; also those where it falls with ``backcomponents``).

    Parameters
    ----------
    variables, nodes, weight, pairwise, frequency, means, method, criteria, bonferroni, \
min_link, max_link, sort_links, decreasing_links
        As in :func:`net_corr` (``criteria`` is ``"value"`` here).
    factorial : str or list of str
        Layout from the loadings: ``"pc"``, ``"nf"``, ``"vf"``, ``"of"`` (see
        :func:`layout_fact`) or ``"null"`` for none. Not compatible with ``layout``.
    components, backcomponents : bool
        Which frames are kept (see above).
    sequence : sequence of float
        The thresholds (0.20 to 1 by 0.01 by default).
    text_filter : (float, float)
        The links whose ``abs(value)`` reaches ``min(text_filter[0], threshold + text_filter[1])``
        show their value as text.
    speed : float
        Speed of the animation (0-100).
    directory : str, optional
        Write the page there.
    **netcoin_args
        The rest goes to :func:`net_corr`.

    Returns
    -------
    netCoin
        The frames merged into one evolving network (see :func:`multigraphCreate`).

    Notes
    -----
    The eigenvector centrality of a graph that is not connected is not unique, so its values
    for the frames with several components can differ from other programs'.
    """
    from .igraph_io import cal_centr, to_igraph
    from .multigraph import evolvingNetwork, frameGraph

    args = _normalize_legacy(netcoin_args, strict=False)
    if sequence is None:
        sequence = np.round(np.arange(20, 101) / 100, 2)
    factorial_list = [factorial] if isinstance(factorial, str) else list(factorial)
    user_layout = args.get("layout")
    if user_layout is not None and factorial_list != ["null"]:
        warnings.warn("Argument factorial is incompatible with layout")
    G = net_corr(variables, weight=weight, pairwise=pairwise, frequency=frequency, means=means,
                 method=method, criteria=criteria, bonferroni=bonferroni, min_link=min_link,
                 max_link=max_link, sort_links=sort_links, decreasing_links=decreasing_links,
                 **({"nodes": nodes} if nodes is not None else {}), **args)
    if user_layout is not None:
        C = user_layout
    else:
        raw = variables.drop(columns=weight) if isinstance(weight, str) else variables
        C = layout_fact(raw, factorial_list)
        if C is not None:
            order = list(G.nodes[G.options["nodeName"]])
            C = ({k: v.reindex(order) for k, v in C.items()} if isinstance(C, dict) else C.reindex(order))

    if G.links is None:
        raise ValueError("there are no links to follow")
    text_filter = list(text_filter)
    if len(text_filter) < 2:
        text_filter.append(0.99)
    bipolar = G.options.get("linkBipolar") is True
    frames: dict = {}
    previous = 0
    for threshold in sequence:
        links = G.links.copy()
        shown = np.abs(links["value"]) >= min(text_filter[0], threshold + text_filter[1])
        links["text"] = np.where(shown, links["value"].map(lambda v: f"{v:.2f}"), "")
        rule = f"{'abs(value)' if bipolar else 'value'} > {threshold:g}"
        changes = dict(links=links, link_filter=rule, lwidth="value", ltext="text",
                       main=f"Correlation: {threshold:g}")
        if C is not None:
            changes["layout"] = C
        if means:
            changes["size"] = "mean"
        H = G.add(**changes)
        if H.links is not None and "hidden" in H.links.columns:
            H.links = H.links[~H.links["hidden"]].reset_index(drop=True)
        elif H.links is not None:
            H.links["hidden"] = False
        name = H.options["nodeName"]
        if H.links is not None and len(H.links):
            ends = set(H.links["Source"]) | set(H.links["Target"])
            H.nodes = H.nodes[H.nodes[name].isin(ends)].reset_index(drop=True)
            membership = [k + 1 for k in to_igraph(H).connected_components().membership]
            count = sum(1 for k in set(membership) if membership.count(k) > 1)
        else:
            membership, count = [], 0
        if ((count > previous or not components) or (count < previous and components and backcomponents)) \
                and H.links is not None and len(H.links):
            with warnings.catch_warnings():       # the frames are often not connected
                warnings.filterwarnings("ignore", message="Some eigenvector centralities")
                central = cal_centr(H)
            H.nodes = pd.concat([H.nodes.reset_index(drop=True), central.reset_index(drop=True)], axis=1)
            H.nodes["component"] = _component_labels(membership)
            H.nodes = H.nodes.rename(columns={"degree": "Degree"})
            H.links["ratio"] = np.abs(H.links["value"]) / threshold
            frames[f"C{threshold:g}"] = H
        previous = count
    lineplots = ["component", "Degree", "closeness", "betweenness", "eigen", "ratio"]
    merged = evolvingNetwork(frames, speed=speed, lineplots=lineplots)
    if directory is not None:
        frameGraph(frames, directory, 0, speed, False, lineplots)
    return merged
