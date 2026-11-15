"""Conversion between netCoin objects and igraph graphs (needs ``python-igraph``)."""
from __future__ import annotations

import colorsys
from typing import Sequence

import numpy as np
import pandas as pd

from .layouts import require_igraph
from .webgraph import category_colors


def to_color_scale(items) -> list[str]:
    """Colours for a vector: shades of red for numbers, the category palette otherwise."""
    values = pd.Series(items)
    if pd.api.types.is_numeric_dtype(values):
        low, high = values.min(), values.max()
        scaled = (values - low) / (high - low) if high != low else values * 0.0
        out = []
        for v in scaled:
            r, g, b = colorsys.hsv_to_rgb(1.0, 1.0, 0.0 if pd.isna(v) else float(v))
            out.append("#%02X%02X%02X" % (round(r * 255), round(g * 255), round(b * 255)))
        return out
    return category_colors(items)


def to_igraph(net):
    """The network of a :class:`netCoin` object as an ``igraph.Graph``.

    Port of ``toIgraph``. The vertices carry the node table (``name``, ``label``,
    ``size``, ``color``, ``shape``... taken from the columns the object uses for them;
    ``fx``/``fy`` become ``x``/``y``) and the edges the link table (``width``,
    ``weight``, ``color``, ``label`` likewise). The graph is directed if
    ``show_arrows`` was set.
    """
    ig = require_igraph()
    options = net.options
    nodes = net.nodes.copy()
    links = net.links.copy() if net.links is not None else pd.DataFrame(
        {"Source": pd.Series(dtype=str), "Target": pd.Series(dtype=str)})
    directed = bool(options.get("showArrows", False))

    for target, option in (("name", "nodeName"), ("label", "nodeLabel"), ("label.cex", "nodeLabelSize"),
                           ("size", "nodeSize"), ("color", "nodeColor"), ("shape", "nodeShape")):
        column = options.get(option)
        if column is not None and column in nodes.columns:
            nodes[target] = nodes[column]
    nodes = nodes.rename(columns={"fx": "x", "fy": "y"})
    nodes = nodes[["name"] + [c for c in nodes.columns if c != "name"]]
    nodes["name"] = nodes["name"].astype(str)

    source, target_col = options.get("linkSource", "Source"), options.get("linkTarget", "Target")
    ordered = [source, target_col] + [c for c in links.columns if c not in (source, target_col)]
    links = links[ordered].copy()
    links.columns = ["from", "to"] + ordered[2:]
    for key, option in (("width", "linkWidth"), ("weight", "linkWeight"), ("color", "linkColor"),
                        ("label", "linkText")):
        column = options.get(option)
        if column is not None and column in links.columns:
            links[key] = links[column]

    if "color" in nodes.columns:
        nodes["color"] = to_color_scale(nodes["color"])
    if "color" in links.columns and len(links):
        links["color"] = to_color_scale(links["color"])

    for frame in (nodes, links):
        for col in frame.columns:
            if isinstance(frame[col].dtype, pd.CategoricalDtype):
                frame[col] = frame[col].astype(object)
    return ig.Graph.DataFrame(links, directed=directed, vertices=nodes, use_vids=False)


def from_igraph(graph, **kwargs):
    """A :class:`netCoin` object from an ``igraph.Graph``.

    Port of ``fromIgraph``. The vertex attributes become node columns and the edge
    attributes link columns; those called ``label``, ``label.cex``, ``size``, ``color``
    and ``shape`` (vertices) or ``width``, ``weight``, ``color`` and ``label`` (edges)
    are used for the corresponding argument unless you give it. Other keyword
    arguments go to :class:`netCoin`.
    """
    from .core import netCoin

    require_igraph()
    arguments = dict(kwargs)
    if "main" not in arguments and "name" in graph.attributes():
        arguments["main"] = graph["name"]
    if "show_arrows" not in arguments and "showArrows" not in arguments and graph.is_directed():
        arguments["show_arrows"] = True

    names = graph.vs["name"] if "name" in graph.vs.attributes() else [str(i + 1) for i in range(graph.vcount())]
    name_column = arguments.get("name") or "name"
    arguments["name"] = name_column
    nodes = pd.DataFrame({name_column: [str(n) for n in names]})
    node_args = {"label": "label", "label.cex": "label_size", "size": "size", "color": "color",
                 "shape": "shape"}
    for attribute in graph.vs.attributes():
        nodes[attribute] = graph.vs[attribute]
        if attribute in node_args and node_args[attribute] not in arguments:
            arguments[node_args[attribute]] = attribute
    nodes[name_column] = [str(n) for n in names]

    edges = graph.get_edgelist()
    labels = [str(n) for n in names]
    links = pd.DataFrame({"Source": [labels[a] for a, _ in edges], "Target": [labels[b] for _, b in edges]})
    link_args = {"width": "lwidth", "weight": "lweight", "color": "lcolor", "label": "ltext"}
    for attribute in graph.es.attributes():
        links[attribute] = graph.es[attribute]
        if attribute in link_args and link_args[attribute] not in arguments:
            arguments[link_args[attribute]] = attribute

    return netCoin(nodes, links if len(links) else None, **arguments)


def cal_centr(graph, measures: Sequence[str] | str = ("degree", "wdegree", "closeness", "betweenness", "eigen"),
              order: str = "") -> pd.DataFrame:
    """Centrality of every node (port of ``calCentr``).

    Parameters
    ----------
    graph : netCoin or igraph.Graph
    measures : str or list of str
        ``"degree"``, ``"wdegree"`` (strength; weights in the edge attribute ``N``, as in
        R), ``"closeness"``, ``"betweenness"``, ``"eigen"`` or ``"all"`` (only the first
        letter counts). In directed graphs ``indegree`` and ``outdegree`` come with
        ``degree``.
    order : str
        Sort the result by this column, decreasing.

    Returns
    -------
    DataFrame
        One row per node: ``nodes`` and a column per measure. (R also returns the
        centralization of each measure; that is not computed here.)
    """
    require_igraph()
    if hasattr(graph, "nodes") and hasattr(graph, "options"):
        graph = to_igraph(graph)
    if isinstance(measures, str):
        measures = [measures]
    measures = list(measures)
    if any(str(m).lower() == "all" for m in measures):
        measures = ["degree", "wdegree", "closeness", "betweenness", "eigen"]
    letters = [str(m)[:1].lower() for m in measures]
    names = graph.vs["name"] if "name" in graph.vs.attributes() else list(range(graph.vcount()))
    out = pd.DataFrame({"nodes": names})
    # R: weights = E(graph)$N, which igraph replaces by the "weight" attribute when it is NULL
    attributes = graph.es.attributes()
    weights = "N" if "N" in attributes else ("weight" if "weight" in attributes else None)

    if "d" in letters:
        out["degree"] = graph.degree(loops=True)
        if "w" in letters:
            out["wdegree"] = graph.strength(weights=weights, loops=True)
        if graph.is_directed():
            out["indegree"] = graph.degree(mode="in", loops=True)
            if "w" in letters:
                out["windegree"] = graph.strength(mode="in", weights=weights, loops=True)
            out["outdegree"] = graph.degree(mode="out", loops=True)
            if "w" in letters:
                out["woutdegree"] = graph.strength(mode="out", weights=weights, loops=True)
    if "c" in letters:
        # R's centr_clo() looks at the paths that start at each node ("out") by default
        out["closeness"] = graph.closeness(mode="out" if graph.is_directed() else "all", normalized=True)
    if "b" in letters:
        out["betweenness"] = graph.betweenness(directed=True)
    if "e" in letters:
        undirected = graph.as_undirected(mode="collapse") if graph.is_directed() else graph
        out["eigen"] = undirected.eigenvector_centrality(scale=True)
    if order in measures or order in out.columns:
        out = out.sort_values(order, ascending=False, kind="stable").reset_index(drop=True)
    return out
