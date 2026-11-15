"""Save a netCoin object for other programs: Pajek and GraphML (ports of ``savePajek``, ``saveGhml``)."""
from __future__ import annotations

import os
from typing import Sequence

import numpy as np
import pandas as pd

from .dichotomize import _collation_key


def _number(x) -> str:
    """A value as R's ``paste`` writes it (up to 15 significant digits)."""
    if isinstance(x, (float, np.floating)):
        if np.isnan(x):
            return "NA"
        return f"{x:.15g}"
    return str(x)


def _levels(values: pd.Series) -> dict:
    """Code (1, 2...) of each value, as R's ``as.numeric(as.factor(x))``."""
    present = values.dropna().unique().tolist()
    if pd.api.types.is_numeric_dtype(values):
        ordered = sorted(present)
    else:
        ordered = sorted(map(str, present), key=_collation_key)
    return {v: i + 1 for i, v in enumerate(ordered)}


def save_pajek(net, file: str = "file.net", arcs: Sequence[str] | None = None,
               edges: Sequence[str] | None = None, partitions: Sequence[str] | None = None,
               vectors: Sequence[str] | None = None, encoding: str = "utf-8") -> str:
    """Write a netCoin object as a Pajek file.

    Port of ``savePajek``.

    Parameters
    ----------
    net : netCoin
    file : str
        Name of the file; ``.net`` is added if it has no extension, and a project
        file ``.paj`` is written instead when there are partitions or vectors.
    arcs, edges : list of str
        Link columns written as weights of directed (``*Arcs``) and undirected
        (``*Edges``) links, one block for each.
    partitions, vectors : list of str
        Node columns written as partitions (categories coded 1, 2...) and as vectors
        (missing values as 0).

    Returns
    -------
    str
        The path written.
    """
    nodes = net.nodes
    links = net.links if net.links is not None else pd.DataFrame()
    name = net.options["nodeName"]

    def check(items, frame, label):
        missing = set(items or []) - set(frame.columns)
        if missing:
            raise ValueError(f"At least one {label} is not amongst {', '.join(map(str, frame.columns))}.")

    check(partitions, nodes, "partition")
    check(vectors, nodes, "vector")
    check(arcs, links, "arc")
    check(edges, links, "edge")

    if "." not in os.path.basename(file):
        file += ".net"
    if partitions or vectors:
        file = file[:-4] + ".paj" if file.endswith(".net") else file

    number = {n: i + 1 for i, n in enumerate(nodes[name].astype(str))}
    lines = ["*Network " + str(net.options.get("main", "")),
             f"*Vertices {len(nodes)}"]
    lines += [f'{i} "{n}" ' for i, n in enumerate(nodes[name].astype(str), 1)]

    if len(links):
        ends = [(number[str(a)], number[str(b)]) for a, b in zip(links[net.options["linkSource"]],
                                                                 links[net.options["linkTarget"]])]
    else:
        ends = []

    count = 1
    for kind, columns in (("Arcs", arcs), ("Edges", edges)):
        for weights in columns or []:
            lines.append(f'*{kind} : {count} "{weights}"')
            lines += [f"{a} {b} {_number(w)}" for (a, b), w in zip(ends, links[weights])]
            count += 1
    for partition in partitions or []:
        codes = _levels(nodes[partition])
        lines += [f"*Partition {partition}", f"*Vertices {len(nodes)}"]
        lines += [str(codes.get(v, "NA")) if pd.notna(v) else "NA" for v in nodes[partition]]
    for vector in vectors or []:
        lines += [f"*Vector {vector}", f"*Vertices {len(nodes)}"]
        lines += ["0" if pd.isna(v) else _number(v) for v in nodes[vector]]

    with open(file, "w", encoding=encoding, newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")
    return file


def save_graphml(net, file: str = "netCoin.graphml") -> str:
    """Write a netCoin object as GraphML (port of ``saveGhml``). Needs python-igraph.

    Categorical columns are written as text, since GraphML only has numbers, text and
    logical values.

    Parameters
    ----------
    net : netCoin
        The network.
    file : str
        Path of the file.

    Returns
    -------
    str
        The path written.
    """
    from .core import netCoin
    from .igraph_io import to_igraph

    if not isinstance(net, netCoin):
        raise TypeError("This program only works with netCoin objects")
    if "." not in os.path.basename(file):
        file += ".graphml"
    graph = to_igraph(net)
    for attributes, sequence in ((graph.vs.attributes(), graph.vs), (graph.es.attributes(), graph.es)):
        for attribute in attributes:
            values = sequence[attribute]
            if any(v is None or (isinstance(v, float) and np.isnan(v)) for v in values) and \
                    not all(isinstance(v, (int, float, type(None))) for v in values):
                sequence[attribute] = ["" if v is None else str(v) for v in values]
            elif not all(isinstance(v, (int, float, bool, str, type(None))) for v in values):
                sequence[attribute] = [str(v) for v in values]
    graph.write_graphml(file)
    return file
