"""A network from a survey table: ``sur_coin`` and ``cob_coin`` (ports of ``surCoin`` and ``cobCoin``)."""
from __future__ import annotations

import warnings
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from .coincidence import as_nodes
from .coincidence import coin as _coin
from .core import _normalize_legacy, name_by_language, netCoin
from .dichotomize import _as_character, _collation_key, dichotomize
from .layouts import layout_circle, layout_mca, layout_pca, layout_special
from .similarity import edge_column_name, edge_list
from .utils import frequencyList, getByLanguage, labelList

_METRIC_METHODS = ["Pearson", "Haberman", "Z"]
_LINK_ATTRIBUTES = ("lwidth", "lweight", "lcolor", "ltext")


def _as_list(value) -> list:
    if value is None:
        return []
    if isinstance(value, str) or not hasattr(value, "__iter__"):
        return [value]
    return list(value)


def _prefix(label: str) -> str:
    """What comes before the first colon (``sub(":.*", "", x)`` in R)."""
    return label.partition(":")[0]


def _union(*lists) -> list:
    return list(dict.fromkeys(x for items in lists for x in items))


# ------------------------------------------------------------ correlations ----
def _corr(a: np.ndarray, b: np.ndarray, weight: np.ndarray | None = None) -> np.ndarray:
    """Weighted correlations between the columns of ``a`` and those of ``b`` (R's ``corr``).

    The cases with a missing value in any of ``a``, ``b`` or ``weight`` are left out.
    """
    n = a.shape[0]
    weight = np.full(n, 1.0 / n) if weight is None else np.asarray(weight, dtype=float)
    ok = ~(np.isnan(a).any(axis=1) | np.isnan(b).any(axis=1) | np.isnan(weight))
    a, b, w = a[ok], b[ok], weight[ok] / weight[ok].sum()
    a = a - (a * w[:, None]).sum(axis=0)
    b = b - (b * w[:, None]).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.sqrt(np.outer((a ** 2 * w[:, None]).sum(axis=0), (b ** 2 * w[:, None]).sum(axis=0)))
        return (a * w[:, None]).T @ b / scale


def _corr_pairwise(a: pd.DataFrame, b: pd.DataFrame, weight: np.ndarray | None) -> np.ndarray:
    """The same, each pair with the cases that have both values (R's ``corrp``)."""
    out = np.full((a.shape[1], b.shape[1]), np.nan)
    for i, ci in enumerate(a.columns):
        for j, cj in enumerate(b.columns):
            if ci == cj:
                out[i, j] = 1.0
            else:
                out[i, j] = _corr(b[[cj]].to_numpy(dtype=float), a[[ci]].to_numpy(dtype=float), weight)[0, 0]
    return out


def _metric_edges(matrix: pd.DataFrame) -> dict[str, np.ndarray]:
    """The cells of a (metrics x variables) matrix as links, as ``mats2edges`` does.

    One link per cell, going column by column: ``Source`` is the variable of the column and
    ``Target`` the metric of the row. Of the pairs of metrics only the lower triangle is kept.
    """
    rows, cols = list(matrix.index), list(matrix.columns)
    where = {r: i for i, r in enumerate(rows)}
    source, target, value = [], [], []
    for j, col in enumerate(cols):
        for i, row in enumerate(rows):
            if col in where and where[col] >= i:       # a pair of metrics, once and no diagonal
                continue
            source.append(col)
            target.append(row)
            value.append(matrix.iloc[i, j])
    return {"Source": np.array(source, dtype=object), "Target": np.array(target, dtype=object),
            "value": np.array(value, dtype=float)}


# ---------------------------------------------------------------- dichotomies ----
def _dicho(data: pd.DataFrame, variables: list[str], value) -> pd.DataFrame:
    """1 where the variable takes ``value`` and 0 elsewhere (also where it is missing)."""
    values = _as_list(value)
    if len(values) != len(variables):
        values = [values[0]] * len(variables)
    out = {}
    for variable, target in zip(variables, values):
        wanted = _as_character(pd.Series([target]))[0]
        text = np.array(_as_character(data[variable]), dtype=object)
        out[variable] = ((text == wanted) & data[variable].notna().to_numpy()).astype(int)
    return pd.DataFrame(out, index=data.index)


# --------------------------------------------------------------------- sur_coin ----
def sur_coin(
    data: pd.DataFrame,
    variables: Sequence[str] | None = None,
    commonlabel: str | Sequence[str] | None = None,
    dichotomies: str | Sequence[str] | None = None,
    value_dicho=1,
    metric: str | Sequence[str] | None = None,
    exogenous: str | Sequence[str] | None = None,
    weight: str | Sequence[float] | None = None,
    subsample: bool = False,
    pairwise: bool = False,
    minimum: float = 1,
    maximum: float | None = None,
    sort: bool = False,
    decreasing: bool = True,
    frequency: bool = False,
    percentages: bool = True,
    procedures: str | Sequence[str] = "Haberman",
    criteria: str = "Z",
    bonferroni: bool = False,
    support: float = -np.inf,
    min_link: float = -np.inf,
    max_link: float = np.inf,
    directed: bool = False,
    diagonal: bool = False,
    sort_links: str | None = None,
    decreasing_links: bool = True,
    igraph: bool = False,
    coin: bool = False,
    directory: str | None = None,
    **netcoin_args,
):
    """A network of the categories of a survey (port of ``surCoin`` of the R package).

    Every category of the ``variables`` becomes a node (``"variable:category"``), and the
    links tell which categories go together more (or less) than expected: the pairs are
    counted with :func:`coin` and described with :func:`edge_list`.

    Parameters
    ----------
    data : DataFrame
        One row per case; the variables can be text, categories or numbers. A missing value
        is a case that has none of the categories of that variable.
    variables : list of str, optional
        The qualitative variables (all the columns by default).
    commonlabel : str or list, optional
        Variables whose name stays in the label of their categories (``"variable:category"``);
        in the others only the category is shown. Creates the column ``label``.
    dichotomies : str or list, optional
        Variables to take as yes/no: 1 where they equal ``value_dicho``. They make one node
        each, named after the variable.
    value_dicho
        The value that counts as yes in the ``dichotomies`` (one for all, or one per variable).
    metric : str or list, optional
        Quantitative variables. Each makes a node (its mean, rescaled to 0-100, is its
        ``"%"``) and is correlated with the rest, with the Pearson coefficient, the Haberman
        residual and the p-value of the correlation (``procedures`` and ``criteria`` are
        limited to those three).
    exogenous : str or list, optional
        Variables (or ``metric``/``dichotomies``) that cannot be explained one by another:
        the links between two of them are hidden (``chaine == "No"``).
    weight : str or sequence, optional
        The weight of each case, or the name of the column that holds it. With it, the cases
        with a missing value are left out (unless ``pairwise``).
    subsample, pairwise, minimum, maximum, sort, decreasing
        As in :func:`coin`; ``minimum`` also limits the categories. ``maximum`` is the number
        of cases by default.
    frequency, percentages : bool
        Put the frequency and/or the percentage of each category in the nodes. With
        ``percentages`` the nodes are sized by ``"%"``.
    procedures, criteria, bonferroni, support
        As in :func:`edge_list`; ``procedures`` also takes the link attributes that
        ``lwidth``, ``lweight``, ``lcolor`` and ``ltext`` refer to.
    min_link, max_link : float
        Range of ``criteria`` kept as links (``minL`` and ``maxL`` in R). For a p-value
        ``max_link`` is 0.5 unless you say otherwise.
    directed, diagonal : bool
        As in :func:`edge_list`.
    sort_links, decreasing_links
        Order of the links (``sortL`` and ``decreasingL`` in R).
    igraph : bool
        Return an ``igraph.Graph`` instead of a :class:`netCoin`.
    coin : bool
        Return the :class:`Coin` object of the categories and stop there.
    directory : str, optional
        Write the web page there.
    **netcoin_args
        The rest goes to :class:`netCoin` (``nodes``, ``name``, ``language``, ``level``,
        ``color``, ``layout``, ``lwidth``, ``main``...). ``nodes`` restricts and orders the
        categories and adds its columns; ``layout`` can also be ``"mc"`` (correspondence
        analysis) or ``"pc"`` (principal components of the correlations), or a matrix with one
        row per row of ``nodes``.

    Returns
    -------
    netCoin, Coin or igraph.Graph
    """
    args = _normalize_legacy(netcoin_args, strict=False)
    language = args.setdefault("language", "en")
    language = language[0] if isinstance(language, (list, tuple)) else language
    args["language"] = language
    level = args.pop("level", 0.95)
    if directory is None:
        directory = args.pop("directory", None)
    else:
        args.pop("directory", None)

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a DataFrame")
    data = data.reset_index(drop=True)
    variables = list(data.columns) if variables is None else _as_list(variables)
    var_order = list(variables)
    metric_vars = _as_list(metric)
    exogenous = _as_list(exogenous)

    if criteria in ("Z", "hyp") and max_link == np.inf:
        max_link = 0.5
    procedures = _union(_as_list(procedures), [args[k] for k in _LINK_ATTRIBUTES if args.get(k)])

    if metric_vars:
        variables = [v for v in variables if v not in metric_vars]
        procedures = [p for p in _union(procedures, ["Z"]) if p in _METRIC_METHODS]
        criteria = [c for c in _union(_as_list(criteria), ["Z"]) if c in _METRIC_METHODS][0]

    nodes = args.get("nodes")
    if isinstance(nodes, dict):
        nodes = pd.DataFrame(nodes)
    name = args["name"] = name_by_language(args.get("name"), language, nodes)
    if nodes is not None:
        if name not in nodes.columns:
            raise ValueError(f"name: '{name}' column missing in nodes")
        nodes = nodes.reset_index(drop=True)

    dichotomies = _as_list(dichotomies)
    if isinstance(weight, str):
        variables = [v for v in variables if v != weight]
    missing = [v for v in _union(metric_vars, dichotomies, variables) + ([weight] if isinstance(weight, str) else [])
               if v not in data.columns]
    if missing:
        raise ValueError("not in data: " + ", ".join(map(str, missing)))

    # a dichotomy that never changes tells nothing
    dichotomies = [d for d in dichotomies if data[d].nunique(dropna=False) != 1]
    all_vars = _union(metric_vars, dichotomies, variables)

    # weights and incomplete cases
    weights = None
    if weight is not None:
        if isinstance(weight, str):
            weights = data[weight].to_numpy(dtype=float)
        else:
            weights = np.asarray(weight, dtype=float)
            if len(weights) != len(data):
                raise ValueError("Weights have not the correct dimensions")
        if not pairwise:
            keep = data[all_vars].notna().all(axis=1).to_numpy() & ~np.isnan(weights)
            data, weights = data[keep].reset_index(drop=True), weights[keep]
    data = data[all_vars]
    args["scenarios"] = int((data.notna().sum(axis=1) > 0).sum())

    if "size" not in args and percentages:
        args["size"] = "%"

    # the table of incidences: dichotomies first, then the categories of the rest
    incidences = dichos = None
    if dichotomies:
        dichos = _dicho(data, dichotomies, value_dicho)
        variables = [v for v in variables if v not in dichotomies]
    if variables:
        incidences = dichotomize(data, variables, sep="", min_freq=minimum, length=0, values=None,
                                 sparse=False, add=False, sort=sort, nas=None)
        if dichos is not None:
            incidences = pd.concat([dichos, incidences], axis=1)
    elif dichos is not None:
        incidences = dichos

    if nodes is not None and incidences is not None:
        wanted = [str(v) for v in nodes[name]]
        absent = [v for v in wanted if v not in incidences.columns and v not in metric_vars]
        incidences = incidences[[v for v in wanted if v in incidences.columns]]
        metric_vars = [m for m in metric_vars if m in wanted]
        if absent:
            warnings.warn(", ".join(absent) + " is/are not present in the data frame.")

    if incidences is None or incidences.shape[1] < 2:
        raise ValueError("There are no more than 1 qualitative category. Try net_corr.")
    order = sorted(range(incidences.shape[1]),
                   key=lambda i: var_order.index(_prefix(incidences.columns[i]))
                   if _prefix(incidences.columns[i]) in var_order else len(var_order))
    incidences = incidences.iloc[:, order]

    cooc = _coin(incidences, minimum, maximum, sort, decreasing, weight=weights,
                 subsample=subsample, pairwise=pairwise)
    if coin:
        return cooc

    # ---- the nodes
    nodes_out = as_nodes(cooc, frequency, percentages, language)
    nodes_out = nodes_out.rename(columns={nodes_out.columns[0]: name})
    nodes_out["variable"] = nodes_out[name].map(_prefix)
    if nodes is not None:
        user = nodes[[c for c in nodes.columns if c not in frequencyList]].copy()
        user[name] = user[name].astype(str)
        nodes_out = nodes_out.merge(user, on=name, how="left", sort=False)
        nodes_out = nodes_out.iloc[sorted(range(len(nodes_out)),
                                          key=lambda i: _collation_key(nodes_out[name].iloc[i]))
                                   ].reset_index(drop=True)
    elif commonlabel is not None:
        label = args["label"] = getByLanguage(labelList, language)
        keep_prefix = set(_as_list(commonlabel))
        full = nodes_out[name].astype(str)
        nodes_out[label] = [t if (t[:t.find(":")] if ":" in t else "") in keep_prefix
                            else (t[t.find(":") + 1:] if ":" in t else t) for t in full]

    # ---- the links
    links = edge_list(cooc, procedures, criteria, level, bonferroni, min_link, max_link, support,
                      directed, diagonal, sort_links, decreasing_links, pairwise)
    links = None if links.empty else links
    for key in _LINK_ATTRIBUTES:
        if args.get(key):
            args[key] = edge_column_name(args[key])

    # ---- the layout
    layout = args.get("layout")
    layout_matrix = layout is not None and not isinstance(layout, (str, list, tuple, dict))
    plain_matrix = layout_matrix and not metric_vars
    if plain_matrix:
        args["layout"] = _matrix_layout(layout, nodes, name, nodes_out)
    elif layout is not None:
        if isinstance(layout, (str, list, tuple, dict)):
            names = list(nodes_out[name])
            args["layout"] = layout_special(layout, {
                "mc": lambda: layout_mca(incidences).reindex(names),
                "pc": lambda: layout_pca(cooc).reindex(names),
            })
        elif metric_vars:
            args["layout"] = None

    # ---- quantitative variables
    if metric_vars:
        nodes_out, links = _add_metric(nodes_out, links, data, incidences, metric_vars, nodes, name,
                                       percentages, procedures, criteria, min_link, max_link,
                                       weights, pairwise)
        if layout_matrix:
            args["layout"] = _matrix_layout(layout, nodes, name, nodes_out)

    if exogenous and links is not None:
        exogenous2 = set(exogenous) & set(_union(metric_vars, dichotomies))
        known = set(exogenous)

        def outside(column):
            return column.map(lambda s: (_prefix(s) in known if ":" in s else "" in known)
                              or s in exogenous2)
        links["chaine"] = np.where(outside(links["Source"]) & outside(links["Target"]), "No", "Yes")
        args["link_filter"] = _only_chained(args.get("link_filter"))

    args["nodes"], args["links"] = nodes_out, links
    net = netCoin(directory=directory, **args)
    if igraph:
        from .igraph_io import to_igraph

        return to_igraph(net)
    return net


def _only_chained(rule):
    """``rule`` and also ``chaine == "Yes"``."""
    if rule is None:
        return "chaine == 'Yes'"
    if isinstance(rule, str):
        return f"({rule}) & chaine == 'Yes'"
    if callable(rule):
        return lambda frame: np.asarray(rule(frame), dtype=bool) & (frame["chaine"] == "Yes").to_numpy()
    values = np.asarray(rule, dtype=bool)
    return lambda frame: values & (frame["chaine"] == "Yes").to_numpy()


def _matrix_layout(layout, nodes, name, nodes_out):
    """A matrix with one row per row of ``nodes``, in the order of the nodes of the network."""
    if nodes is None:
        warnings.warn("layout must be applied to the nodes variable")
        return layout
    frame = pd.DataFrame(np.asarray(layout, dtype=float) if not isinstance(layout, pd.DataFrame)
                         else layout.to_numpy(dtype=float))
    if len(frame) != len(nodes):
        warnings.warn("layout must have a coordinate per node")
        return layout
    frame.index = [str(v) for v in nodes[name]]
    return frame.reindex(list(nodes_out[name])).to_numpy()


def _add_metric(nodes_out, links, data, incidences, metric_vars, nodes, name, percentages,
                procedures, criteria, min_link, max_link, weights, pairwise):
    """The metric variables as nodes, and their correlations as links."""
    if percentages:
        nodes_out["mean"] = nodes_out["%"] / 100
    nodes_out["min"] = 0
    nodes_out["max"] = 1
    complete = data[metric_vars].dropna()
    means, mins, maxs = complete.mean(), complete.min(), complete.max()
    extra = pd.DataFrame({name: list(means.index), "mean": means.to_numpy(), "min": mins.to_numpy(),
                          "max": maxs.to_numpy(), "%": ((means - mins) / (maxs - mins) * 100).to_numpy(),
                          "variable": list(means.index)})
    if nodes is not None:
        extra = extra[extra[name].isin([str(v) for v in nodes[name]])].reset_index(drop=True)
        for col in extra[name]:
            row = nodes[nodes[name].astype(str) == col].iloc[0]
            for column in nodes.columns:
                if column not in extra.columns:
                    extra[column] = nodes[column].iloc[:0].reindex(extra.index)   # empty, of the dtype of the column
                extra.loc[extra[name] == col, column] = row[column]
    nodes_out = pd.concat([nodes_out, extra], ignore_index=True, sort=False)

    # correlations of each metric with the categories and with the other metrics
    metrics = data[metric_vars].astype(float)
    others = pd.concat([incidences, metrics], axis=1)
    if pairwise:
        values = _corr_pairwise(metrics, others, weights)
    else:
        values = _corr(metrics.to_numpy(), others.to_numpy(dtype=float), weights)
    r = pd.DataFrame(values, index=metric_vars, columns=list(others.columns))
    universe = _union([str(v) for v in nodes[name]] if nodes is not None else [],
                      list(incidences.columns), metric_vars)
    r = r.loc[[x for x in universe if x in r.index], [x for x in universe if x in r.columns]]

    n = len(incidences)
    edges = _metric_edges(r)
    r_values = edges["value"]
    methods = _union(procedures, [criteria])
    table = {"Source": edges["Source"], "Target": edges["Target"], "Pearson": r_values}
    with np.errstate(divide="ignore", invalid="ignore"):
        if "Haberman" in methods:
            table["Haberman"] = r_values * np.sqrt(n)
        if "Z" in methods:
            # the usual t of a correlation; R divides by sqrt((1 - r) / (n - 2)) (see README)
            t = r_values * np.sqrt((n - 2) / (1 - np.minimum(1, r_values) ** 2))
            table["Z"] = student_t.sf(t, n - 2)
    new = pd.DataFrame(table)[["Source", "Target"] + methods]
    new = new[(new[criteria] > min_link) & (new[criteria] < max_link)].rename(columns={"Z": "p(Z)"})
    links = new.reset_index(drop=True) if links is None else pd.concat([links, new], ignore_index=True, sort=False)
    return nodes_out, links


# --------------------------------------------------------------------- cob_coin ----
def cob_coin(
    data: pd.DataFrame,
    variables: Sequence[str] | None = None,
    degree: float = 0,
    significance: float = 0.05,
    **kwargs,
):
    """A CoWeb network (Upton 2000): the categories on a circle, linked if they go together.

    Port of ``cobCoin``. It builds the network with :func:`sur_coin`, keeps the links between
    categories of different variables whose p-value is below ``significance / 2`` (they go
    together more than expected) or above ``1 - significance / 2`` (less), and puts the nodes
    on a circle.

    Parameters
    ----------
    data : DataFrame
        One row per case.
    variables : list of str, optional
        The variables to use (all the columns by default).
    degree : float
        Angle in degrees of the first node of the circle.
    significance : float
        Two-sided significance level of the links kept.
    **kwargs
        As in :func:`sur_coin`. By default nodes are coloured and grouped by ``variable``, and
        links are as wide and coloured as their ``Haberman`` residual (``link_bipolar``).
    """
    args = _normalize_legacy(kwargs, strict=False)
    args.setdefault("color", "variable")
    args.setdefault("group", "variable")
    args.setdefault("lwidth", "Haberman")
    args.setdefault("lcolor", "Haberman")
    args.setdefault("link_bipolar", True)
    args["group_text"] = True
    args["max_link"] = 1
    args["commonlabel"] = ""
    directory = args.pop("directory", None)
    to_graph = args.pop("igraph", False)
    args.pop("coin", None)
    variables = list(data.columns) if variables is None else _as_list(variables)

    net = sur_coin(data[variables], **args)
    links = net.links
    if links is not None:
        if "p(Z)" not in links.columns:
            raise ValueError("cob_coin needs the p-value of the links: use criteria='Z' and "
                             "keep 'Z' among the procedures")
        different = links["Source"].map(_prefix) != links["Target"].map(_prefix)
        links = links[different]
        p = links["p(Z)"]
        links = links[(p < significance / 2) | (p > 1 - significance / 2)].reset_index(drop=True)
        links = links if len(links) else None
    changes = {"layout": layout_circle(net.nodes, deg=degree, name=net.options["nodeName"]),
               "links": links}
    if directory is not None:
        changes["directory"] = directory
    net = net.add(**changes)
    if to_graph:
        from .igraph_io import to_igraph

        return to_igraph(net)
    return net
