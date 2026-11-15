"""Log-linear networks: ``log_coin`` (port of ``logCoin``)."""
from __future__ import annotations

import math
import warnings
from itertools import combinations
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .core import _normalize_legacy, name_by_language, netCoin
from .dichotomize import _as_character, _collation_key
from .utils import getByLanguage, labelList, nameList


# ------------------------------------------------------------------ formulas ----
def _comb(endogenous: list[str], exogenous: list[str] | None, order: int) -> list[list[str]]:
    """The products of variables ``comb`` puts in the formula (each product is a list of names)."""
    exogenous = list(exogenous or [])
    endogenous = [v for v in endogenous if v not in exogenous]
    order_endo = min(order, len(endogenous))
    products = [list(c) for c in combinations(endogenous, order_endo)]
    terms: list[list[str]] = []
    if not exogenous or len(endogenous) >= order:
        terms += products
    for exo in exogenous:
        if order < len(endogenous) + len(exogenous):
            products = [list(c) for c in combinations(endogenous, min(order_endo, order - 1))]
        terms += [p + [exo] for p in products] if products and products != [[]] else [[exo]]
    return terms


def _expand(products: list[list[str]], variables: list[str]) -> list[tuple[str, ...]]:
    """The terms of ``a*b + a*c``: every subset of each product, by degree and first appearance."""
    seen: list[tuple[str, ...]] = []
    index = {v: i for i, v in enumerate(variables)}
    for product in products:
        product = sorted(dict.fromkeys(product), key=index.get)
        subsets = []
        for v in product:                          # R expands a*b*c as a, b, a:b, c, a:c, b:c, a:b:c
            subsets = subsets + [(v,)] + [s + (v,) for s in subsets]
        for s in subsets:
            if s not in seen:
                seen.append(s)
    return sorted(seen, key=len)


def deparse_formula(response: str, products: list[list[str]], cutoff: int = 60) -> str:
    """What ``Reduce(paste, deparse(formula))`` gives in R for ``response ~ a * b + c * d``.

    R breaks the line after an operator when the line is longer than 60 characters, and the
    lines that follow are indented four spaces; ``paste`` joins them with one more space.
    """
    lines: list[str] = []
    state = {"line": "", "indent": 0, "start": True}

    def write(text: str):
        if state["start"]:
            state["line"] += "    " * min(state["indent"], 4) + "  " * max(state["indent"] - 4, 0)
            state["start"] = False
        state["line"] += text

    def check(flag: dict):
        if len(state["line"]) > cutoff:
            if not flag["on"]:
                flag["on"] = True
                state["indent"] += 1
            lines.append(state["line"])
            state["line"], state["start"] = "", True

    def binary(left, op, right_writer, flag):
        left()
        write(f" {op} ")
        check(flag)
        right_writer()

    # the tree: response ~ ((t1 + t2) + t3) ..., each t a left-nested product
    def product(vars_: list[str]):
        def build(k):
            if k == 0:
                return lambda: write(vars_[0])
            inner = build(k - 1)

            def go():
                flag = {"on": False}
                binary(inner, "*", lambda: write(vars_[k]), flag)
                if flag["on"]:
                    state["indent"] -= 1
            return go
        return build(len(vars_) - 1)

    def total(k):
        if k == 0:
            return product(products[0])
        inner = total(k - 1)

        def go():
            flag = {"on": False}
            binary(inner, "+", product(products[k]), flag)
            if flag["on"]:
                state["indent"] -= 1
        return go

    flag = {"on": False}
    binary(lambda: write(response), "~", total(len(products) - 1), flag)
    lines.append(state["line"])
    return " ".join(lines)


# -------------------------------------------------------------------- tables ----
def _ipf(table: np.ndarray, margins: list[tuple[int, ...]], eps: float = 0.1, iterations: int = 40) -> np.ndarray:
    """Fitted table of a hierarchical log-linear model, by the iterative proportional fitting of ``loglin``."""
    fit = np.ones_like(table, dtype=float)
    axes = list(range(table.ndim))
    observed = [table.sum(axis=tuple(a for a in axes if a not in m), keepdims=True) for m in margins]
    for _ in range(iterations):
        worst = 0.0
        for m, obs in zip(margins, observed):
            current = fit.sum(axis=tuple(a for a in axes if a not in m), keepdims=True)
            worst = max(worst, float(np.abs(obs - current).max()))
            with np.errstate(divide="ignore", invalid="ignore"):
                factor = np.where(current > 0, obs / current, 0.0)
            fit = fit * factor
        if worst < eps:
            break
    return fit


def _maximal(terms: list[tuple[str, ...]]) -> list[tuple[str, ...]]:
    sets = [set(t) for t in terms]
    return [t for t, s in zip(terms, sets) if not any(s < other for other in sets)]


def _loglm(table: np.ndarray, names: list[str], terms: list[tuple[str, ...]]) -> tuple[float, int]:
    """L2 and degrees of freedom of the model whose terms are given (as ``MASS::loglm``)."""
    position = {n: i for i, n in enumerate(names)}
    margins = [tuple(position[v] for v in t) for t in _maximal(terms)]
    fit = _ipf(table, margins)
    keep = (table * fit) > 0
    lrt = 2 * float((table[keep] * np.log(table[keep] / fit[keep])).sum())
    used = set()
    for m in margins:
        for r in range(1, len(m) + 1):
            for sub in combinations(sorted(m), r):
                used.add(sub)
    df = table.size - 1 - sum(int(np.prod([table.shape[i] - 1 for i in sub])) for sub in used)
    return lrt, df


# ---------------------------------------------------------------- the glm ----
def _contrast(levels: int, last: bool) -> np.ndarray:
    """``contr.first`` (or ``contr.last``): the first (last) level is coded -1 in every column."""
    k = levels - 1
    if last:
        return np.vstack([np.eye(k), -np.ones((1, k))])
    return np.vstack([-np.ones((1, k)), np.eye(k)])


def _design(dt: pd.DataFrame, names: list[str], levels: dict[str, list[str]], terms, last: dict[str, bool]):
    """Model matrix and column names for ``terms`` with the contrasts chosen in ``last``."""
    columns, labels = [np.ones(len(dt))], ["(Intercept)"]
    coded = {}
    for v in names:
        k = len(levels[v])
        codes = dt[v].map({lv: i for i, lv in enumerate(levels[v])}).to_numpy()
        matrix = _contrast(k, last[v])[codes]
        kept = levels[v][:-1] if last[v] else levels[v][1:]
        coded[v] = (matrix, [v + lv for lv in kept])
    for term in terms:
        mats, labs = [np.ones((len(dt), 1))], [""]
        for v in term:                                    # the first factor varies fastest
            m, n = coded[v]
            mats = [a * m[:, [j]] for j in range(m.shape[1]) for a in mats]
            labs = [(x + ":" if x else "") + n[j] for j in range(len(n)) for x in labs]
        columns += [mats[i][:, 0] for i in range(len(mats))]
        labels += labs
    return np.column_stack(columns), labels


def _poisson(X: np.ndarray, y: np.ndarray, labels: list[str]) -> pd.DataFrame:
    """``summary(glm(y ~ X - 1, family = poisson))``: estimates, errors, z and p of the estimable columns.

    It follows ``glm.fit``: it starts at ``y + 0.1`` and stops when the deviance changes by less
    than 1e-8 (relative), and the covariance comes from the weights of the last iteration.
    """
    def deviance(mu):
        with np.errstate(divide="ignore", invalid="ignore"):
            return 2 * (np.where(y > 0, y * np.log(y / mu), 0.0) - (y - mu)).sum()

    mu = y + 0.1
    eta = np.log(mu)
    old = deviance(mu)
    for _ in range(25):
        w = mu                                     # weights of this iteration (from the previous fit)
        z = eta + (y - mu) / mu
        keep = _estimable(X * np.sqrt(w)[:, None])
        Xk = X[:, keep]
        A = Xk.T @ (Xk * w[:, None])
        beta = np.linalg.solve(A, Xk.T @ (w * z))
        eta = Xk @ beta
        mu = np.exp(eta)
        new = deviance(mu)
        if abs(new - old) / (abs(new) + 0.1) < 1e-8:
            break
        old = new
    cov = np.linalg.inv(A)
    se = np.sqrt(np.diag(cov))
    zval = beta / se
    names = [l for l, k in zip(labels, keep) if k]
    return pd.DataFrame({"Estimate": beta, "Std. Error": se, "z value": zval,
                         "Pr(>|z|)": 2 * stats.norm.sf(np.abs(zval))}, index=names)


def _estimable(X: np.ndarray, tol: float = 1e-7) -> np.ndarray:
    """The columns that are not a combination of the previous ones (as the QR of ``lm.fit`` finds)."""
    keep = np.zeros(X.shape[1], bool)
    basis: list[np.ndarray] = []
    for j in range(X.shape[1]):
        v = X[:, j].astype(float)
        norm = np.linalg.norm(v)
        r = v.copy()
        for b in basis:
            r -= (b @ r) * b
        if norm > 0 and np.linalg.norm(r) > tol * norm:
            basis.append(r / np.linalg.norm(r))
            keep[j] = True
    return keep


# -------------------------------------------------------------------- log_coin ----
# the cascade of contrasts of logCoin: (use the last category?, position of the variable) by number of pivots
_CASCADE = {
    1: [(1, 1)],
    2: [(1, 1), (1, 2), (0, 1)],
    3: [(1, 1), (1, 2), (0, 1), (1, 3), (1, 1), (0, 2), (0, 1)],
    4: [(1, 1), (1, 2), (0, 1), (1, 3), (1, 1), (0, 2), (0, 1), (1, 4), (1, 1), (1, 2), (0, 1), (0, 3), (1, 1),
        (0, 2), (0, 1)],
    5: [(1, 1), (1, 2), (0, 1), (1, 3), (1, 1), (0, 2), (0, 1), (1, 4), (1, 1), (1, 2), (0, 1), (0, 3), (1, 1),
        (0, 2), (0, 1), (1, 5), (1, 1), (1, 2), (0, 1), (1, 3), (1, 1), (0, 2), (0, 1), (0, 4), (1, 1), (1, 2),
        (0, 1), (0, 3), (1, 1), (0, 2), (0, 1)],
}


def _factor_levels(column: pd.Series) -> list[str]:
    """The levels ``as_factor`` gives: numbers sorted, text in order of appearance."""
    if pd.api.types.is_numeric_dtype(column) and not pd.api.types.is_bool_dtype(column):
        return _as_character(pd.Series(sorted(column.dropna().unique())))
    if isinstance(column.dtype, pd.CategoricalDtype):
        return [str(c) for c in column.cat.categories if (column == c).any()]
    return list(dict.fromkeys(_as_character(column)))


def log_coin(
    data: pd.DataFrame,
    variables: Sequence[str] | None = None,
    exogenous: str | Sequence[str] | None = None,
    no_first_cat: str | Sequence[str] | None = None,
    weight: str | Sequence[float] | None = None,
    order: int = 2,
    pairwise: bool = False,
    twotails: bool = False,
    pmax: float = 0.05,
    frequency: bool = False,
    percentage: bool = False,
    directed: bool = False,
    igraph: bool = False,
    **netcoin_args,
):
    """A network of the associations among categorical variables from a log-linear model.

    Port of ``logCoin``. A model with all the interactions up to ``order`` is fitted to the
    table of the variables (Poisson regression on the counts of the cells, with contrasts that
    make the effects deviations from the mean). Every category is a node, its estimate being the
    main effect; the interactions between two categories whose estimate is significant are links,
    and those of three or more categories are joined through a node ``u01``, ``u02``...

    Parameters
    ----------
    data : DataFrame
        One row per case. Spaces in the names of columns are replaced by ``_``.
    variables : list of str, optional
        The variables (all the columns by default).
    exogenous : str or list of str, optional
        Variables that are only used to explain the others: their own interactions are not
        in the model.
    no_first_cat : str or list of str, optional
        Variables whose first category is not shown as a node (the pivot of the contrasts).
        At most five variables can have their first category included.
    weight : str or sequence, optional
        Weight of each case or the column that holds it. The counts are rounded.
    order : int
        Order of the interactions of the model.
    twotails : bool
        Two-sided contrast: significant negative links are kept too and the links are bipolar,
        coloured by their estimate. Otherwise only the positive ones are kept, with the
        one-sided p-value ``Pr(>z)``.
    pmax : float
        Significance level of the links.
    frequency, percentage : bool
        Put the number and/or the proportion of cases of each category in the nodes (and
        size them by it).
    pairwise, directed
        Accepted for compatibility with R; they have no effect.
    igraph : bool
        Return an ``igraph.Graph``.
    **netcoin_args
        The rest goes to :class:`netCoin` (``language``, ``color``, ``lwidth``, ``main``...).

    Returns
    -------
    netCoin or igraph.Graph
    """
    args = _normalize_legacy(netcoin_args, strict=False)
    fix = lambda s: str(s).replace(" ", "_")
    data = data.copy()
    data.columns = [fix(c) for c in data.columns]
    variables = [fix(v) for v in (list(data.columns) if variables is None else _as_list(variables))]
    exogenous = [fix(v) for v in _as_list(exogenous)]
    no_first_cat = [fix(v) for v in _as_list(no_first_cat)]
    variables = list(dict.fromkeys(
        [v for v in dict.fromkeys(variables + no_first_cat) if v not in exogenous] + exogenous))

    # the weights
    if weight is None:
        w = np.ones(len(data))
    elif isinstance(weight, str):
        weight = fix(weight)
        variables = [v for v in variables if v != weight]
        w = data[weight].to_numpy(dtype=float)
    else:
        w = np.asarray(weight, dtype=float)
        if len(w) != len(data):
            raise ValueError("Weights have not the correct dimensions")
    missing = [v for v in variables if v not in data.columns]
    if missing:
        raise ValueError("not in data: " + ", ".join(missing))
    frame = data[variables].copy()
    frame["weight"] = w
    endogenous = [v for v in variables if v not in exogenous]
    products = _comb(endogenous, exogenous, order)
    if not products:
        raise ValueError("there are no interactions to fit: give at least two variables")

    language = args.setdefault("language", "en")
    language = language[0] if isinstance(language, (list, tuple)) else language
    args["language"] = language
    nodes_arg = args.get("nodes")
    name = args["name"] = name_by_language(args.get("name"), language, nodes_arg)

    pivots = [v for v in variables if v not in no_first_cat]
    if len(pivots) > 5:
        raise ValueError("This function doesn't support more than 5 variables with first categories included")
    args["scenarios"] = int((frame.notna().sum(axis=1) > 0).sum())

    # the table of the cells that occur
    appearance = list(dict.fromkeys(v for p in products for v in p))
    complete = frame[appearance + ["weight"]].dropna()
    keys = {v: complete[v] for v in appearance}
    dt = complete.groupby(appearance, sort=False, dropna=True)["weight"].sum().reset_index()
    sort_columns = {}
    for v in reversed(appearance):                     # R's aggregate sorts with the first variable fastest
        if pd.api.types.is_numeric_dtype(dt[v]) and not pd.api.types.is_bool_dtype(dt[v]):
            sort_columns[v] = dt[v].to_numpy(dtype=float)
        else:
            order_of = {t: i for i, t in enumerate(sorted(set(_as_character(dt[v])), key=_collation_key))}
            sort_columns[v] = np.array([order_of[t] for t in _as_character(dt[v])])
    index = np.lexsort(tuple(sort_columns[v] for v in appearance))
    dt = dt.iloc[index].reset_index(drop=True)
    dt["weight"] = np.round(dt["weight"].to_numpy(dtype=float))
    n_cells = len(dt)
    levels = {v: _factor_levels(dt[v]) for v in appearance}
    for v in appearance:
        dt[v] = _as_character(dt[v])

    # the full table (xtabs): levels in alphabetical order
    table_levels = {v: sorted(set(dt[v]), key=_collation_key) if not _is_numeric(complete[v])
                    else sorted(set(dt[v]), key=float) for v in appearance}
    table = np.zeros([len(table_levels[v]) for v in appearance])
    positions = np.column_stack([dt[v].map({lv: i for i, lv in enumerate(table_levels[v])}).to_numpy()
                                 for v in appearance]).astype(int)
    for row, count in zip(positions, dt["weight"]):
        table[tuple(row)] += count
    terms = _expand(products, appearance)
    lrt, df = _loglm(table, appearance, terms)

    cells = 1
    for v in appearance:
        cells *= len(levels[v])
    formula = deparse_formula("weight", products)
    p_value = 1 - stats.chi2.cdf(lrt, df)
    args["note"] = (f"<p>L2= {lrt:3.2f}; d.f.= {df}; p(>X^2)= {p_value:3.3g}.</p><p>Covariance structures: "
                    f"{n_cells:1.0f} ({cells:1.0f}).</p><p>{formula}</p>")

    # coefficients under the cascade of contrasts
    last = {v: False for v in appearance}
    y = dt["weight"].to_numpy(dtype=float)
    X, labels = _design(dt, appearance, levels, terms, last)
    coefs = _poisson(X, y, labels)
    for use_last, position in _CASCADE.get(len(pivots), []):
        last[pivots[position - 1]] = bool(use_last)
        X, labels = _design(dt, appearance, levels, terms, last)
        more = _poisson(X, y, labels)
        coefs = pd.concat([coefs, more[~more.index.isin(coefs.index)]])
    coefs = coefs.loc[sorted(coefs.index, key=_collation_key)]

    # the nodes
    info = []
    for v in variables:
        if v not in levels:
            continue
        for k, lv in enumerate(levels[v]):
            info.append((v + lv, f"{v}:{lv}", v, lv, appearance.index(v) + 1, k + 1))
    informa = pd.DataFrame(info, columns=["name", "eti", "var", "cat", "varn", "labn"]).set_index("name", drop=False)
    main = coefs[~coefs.index.str.contains(":", regex=False) & (coefs.index != "(Intercept)")]
    nodes = pd.DataFrame({"Estimate": np.round(main["Estimate"].to_numpy(), 3), "name": list(main.index)})
    nodes = nodes[nodes["name"].isin(informa.index)]
    nodes = nodes.merge(informa[["name", "eti", "var", "cat"]].reset_index(drop=True), on="name")
    nodes = nodes.sort_values("name", key=lambda s: s.map(_collation_key), kind="stable")
    nodes = nodes[["Estimate", "name", "eti", "var", "cat"]].reset_index(drop=True)
    if frequency or percentage:
        marg = _margins(table, appearance, table_levels, frequency, percentage)
        nodes = nodes.merge(marg, on="name", how="left")
        nodes = nodes[["name"] + [c for c in nodes.columns if c != "name"]]
        nodes = nodes.sort_values("name", key=lambda s: s.map(_collation_key), kind="stable").reset_index(drop=True)
    order_of_vars = {v: i for i, v in enumerate(variables)}
    ordernodes = list(nodes.sort_values("var", key=lambda s: s.map(lambda v: order_of_vars.get(v, 1e9)),
                                        kind="stable")["name"])

    # the links
    threshold = abs(stats.norm.ppf(pmax / 2 if twotails else pmax))
    has_colon = coefs.index.str.contains(":", regex=False)
    if not twotails:
        selected = coefs[(coefs["z value"] > threshold) & has_colon].copy()
        selected["Pr(>|z|)"] = selected["Pr(>|z|)"] / 2
        selected = selected.rename(columns={"Pr(>|z|)": "Pr(>z)"})
    else:
        selected = coefs[(coefs["z value"].abs() > threshold) & has_colon].copy()
        args["link_bipolar"] = True
        args["lcolor"] = "Estimate"
    selected["order"] = selected.index.str.count(":") + 1
    links = selected[selected["order"] == 2].copy()
    links["Source"] = [i.split(":")[0] for i in links.index]
    links["Target"] = [i.split(":")[-1] for i in links.index]
    if len(links):
        pairs = [_orient(s, t, ordernodes) for s, t in zip(links["Source"], links["Target"])]
        links["Source"], links["Target"] = [p[0] for p in pairs], [p[1] for p in pairs]
    links = links.reset_index(drop=True)

    ulinks = selected[selected["order"] > 2]
    if len(ulinks):
        width = math.ceil((len(ulinks) + 1) / 10)
        codes = [f"u{i + 1:0{width}d}" for i in range(len(ulinks))]
        members = [i.split(":") for i in ulinks.index]
        joint = pd.DataFrame({"Source": [c for c, m in zip(codes, members) for _ in m],
                              "Target": [t for m in members for t in m]})
        unodes = pd.DataFrame({"name": codes, "Estimate": math.floor(nodes["Estimate"].min()), "eti": "",
                               "var": [f".I{o}" for o in ulinks["order"]], "cat": ""})
        if frequency or percentage:
            joint_counts = [_joint_frequency(m, informa, table, appearance, table_levels) for m in members]
            if frequency:
                unodes["freq"] = joint_counts
            if percentage:
                unodes["prop"] = np.array(joint_counts) / table.sum()
        nodes = pd.concat([nodes, unodes], ignore_index=True, sort=False)
        extra = ulinks.assign(Source=codes).reset_index(drop=True)
        merged = joint.merge(extra, on="Source", how="left")
        links = pd.concat([links, merged], ignore_index=True, sort=False)

    nodes = nodes.rename(columns={"name": getByLanguage(nameList, language), "eti": getByLanguage(labelList, language)})
    args["nodes"] = nodes
    if len(links):
        links = links[["Source", "Target"] + [c for c in links.columns if c not in ("Source", "Target")]] \
            if False else links
        args["links"] = links
        args.setdefault("degree_filter", 1)
    if "size" not in args:
        args["size"] = "prop" if percentage else "freq" if frequency else "Estimate"
    args.setdefault("color", "var")
    args.setdefault("lwidth", "Estimate")
    args.setdefault("label", getByLanguage(labelList, language))
    net = netCoin(**args)
    if igraph:
        from .igraph_io import to_igraph

        return to_igraph(net)
    return net


def _as_list(value) -> list:
    if value is None:
        return []
    return [value] if isinstance(value, str) or not hasattr(value, "__iter__") else list(value)


def _is_numeric(series: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series)


def _orient(source: str, target: str, order: list[str]) -> tuple[str, str]:
    """The node that comes later in ``order`` is the source (``olinks`` in R)."""
    if order.index(source) < order.index(target):
        return target, source
    return source, target


def _margins(table, names, table_levels, frequency, proportion) -> pd.DataFrame:
    rows = []
    total = table.sum()
    for axis, v in enumerate(names):
        counts = table.sum(axis=tuple(a for a in range(table.ndim) if a != axis))
        for lv, count in zip(table_levels[v], counts):
            rows.append((v + lv, count, count / total))
    out = pd.DataFrame(rows, columns=["name", "freq", "prop"])
    if not frequency:
        out = out.drop(columns="freq")
    if not proportion:
        out = out.drop(columns="prop")
    return out


def _joint_frequency(members, informa, table, names, table_levels) -> float:
    """Cases that have all the categories of an interaction."""
    index = [slice(None)] * table.ndim
    for m in members:
        v, lv = informa.loc[m, "var"], informa.loc[m, "cat"]
        index[names.index(v)] = table_levels[v].index(lv)
    return float(table[tuple(index)].sum())
