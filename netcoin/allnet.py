"""A complete network from a table of incidences (port of netCoin's ``allNet``)."""
from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from .coincidence import as_nodes, coin
from .core import netCoin
from .layouts import layout_mca, layout_pca, layout_special
from .similarity import edge_column_name, edge_list
from .utils import frequencyList, getByLanguage, nameList


def all_net(
    incidences: pd.DataFrame,
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
    level: float = 0.95,
    support: float = -np.inf,
    min_link: float = -np.inf,
    max_link: float = np.inf,
    directed: bool = False,
    diagonal: bool = False,
    sort_links: str | None = None,
    decreasing_links: bool = True,
    igraph: bool = False,
    **netcoin_args,
):
    """Build a :class:`netCoin` network from a table of incidences.

    Port of ``allNet`` from the R package netCoin: it counts the coincidences
    (:func:`coin`), turns the variables into nodes (:func:`as_nodes`), the
    pairs into links (:func:`edge_list`) and draws them with :class:`netCoin`.

    Parameters
    ----------
    incidences : DataFrame
        One row per case, one 0/1 column per variable (see :func:`dichotomize`).
    weight : str or sequence, optional
        Weight of each case, or the name of the column that holds it.
    subsample, pairwise, minimum, maximum, sort, decreasing
        As in :func:`coin`.
    frequency, percentages : bool
        Put the frequency and/or the percentage of each variable in the nodes.
        With ``percentages`` the nodes are sized by ``"%"``.
    procedures, criteria, bonferroni, level, support
        As in :func:`edge_list`; ``procedures`` also gets the link attributes
        that ``lwidth``, ``lweight``, ``lcolor`` and ``ltext`` refer to.
    min_link, max_link : float
        Range of ``criteria`` kept as links (``minL`` and ``maxL`` in R).
    directed, diagonal : bool
        As in :func:`edge_list`.
    sort_links, decreasing_links
        Order of the links (``sortL`` and ``decreasingL`` in R).
    igraph : bool
        Return an ``igraph.Graph`` instead (needs the optional ``igraph`` package).
    **netcoin_args
        Anything else is passed to :class:`netCoin` (``name``, ``nodes``,
        ``language``, ``color``, ``layout``...). ``layout`` can also be
        ``"mc"`` (correspondence analysis) or ``"pc"`` (principal components of
        the correlations); ``directory`` writes the web page.
    """
    if not isinstance(incidences, pd.DataFrame):
        raise TypeError("incidences must be a DataFrame")

    args = dict(netcoin_args)
    language = args.setdefault("language", "en")
    name = args.get("name") or getByLanguage(nameList, language)
    args["name"] = name
    user_nodes = args.get("nodes")
    if isinstance(user_nodes, dict):
        user_nodes = pd.DataFrame(user_nodes)
    if user_nodes is not None:
        if name not in user_nodes.columns:
            raise ValueError(f"name: '{name}' column missing in nodes")
        if user_nodes[name].duplicated().any():
            raise ValueError(f"name: '{name}' column values must be unique")
    if "size" not in args and percentages:
        args["size"] = "%"

    data = incidences
    if isinstance(weight, str):
        weight, data = data[weight].to_numpy(dtype=float), data.drop(columns=weight)
    elif weight is not None:
        weight = np.asarray(weight, dtype=float)
    if not pairwise:  # incomplete cases are dropped, with their weights
        complete = ~data.isna().any(axis=1).to_numpy()
        data = data[complete]
        weight = weight[complete] if weight is not None else None
    data = data.loc[:, data.sum(axis=0, skipna=True) > 0]

    cooc = coin(data, minimum, maximum, sort, decreasing, weight=weight,
                subsample=subsample, pairwise=pairwise)
    if args.get("size") == "frequency":
        frequency = True

    nodes = as_nodes(cooc, frequency, percentages, language)
    nodes = nodes.rename(columns={nodes.columns[0]: name})
    if user_nodes is None:
        args["nodes"] = nodes
    else:
        user = user_nodes.copy()
        user[name] = user[name].astype(str)
        extra = user[[c for c in user.columns if c not in frequencyList]]
        args["nodes"] = nodes.merge(extra, on=name, how="right", sort=False).reset_index(drop=True)

    link_args = [args[k] for k in ("lwidth", "lweight", "lcolor", "ltext") if args.get(k)]
    procedures = [procedures] if isinstance(procedures, str) else list(procedures)
    procedures = list(dict.fromkeys(procedures + link_args))
    args["links"] = edge_list(cooc, procedures, criteria, level, bonferroni, min_link, max_link,
                              support, directed, diagonal, sort_links, decreasing_links)
    for key in ("lwidth", "lweight", "lcolor", "ltext"):
        if args.get(key):
            args[key] = edge_column_name(args[key])

    if args.get("layout") is not None:
        names = list(args["nodes"][name])
        args["layout"] = layout_special(args["layout"], {
            "mc": lambda: layout_mca(data[cooc.labels]).reindex(names),
            "pc": lambda: layout_pca(cooc).reindex(names),
        })
    args["scenarios"] = cooc.n
    net = netCoin(**args)
    if igraph:
        from .igraph_io import to_igraph

        return to_igraph(net)
    return net
