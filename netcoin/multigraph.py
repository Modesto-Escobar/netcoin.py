"""Several graphs in one page (port of rD3plot's ``rd3_multigraph`` and ``evolNetwork_rd3``).

Every graph is written to its own folder (``data/0``, ``data/1``...) and the main
page only holds a selector (``multigraph.js``); the ``parallel`` mode shows
several of them side by side, and the ``frame`` mode merges networks into one
that evolves from frame to frame.
"""
from __future__ import annotations

import os
import warnings
import webbrowser
from copy import deepcopy
from shutil import copytree
from typing import Sequence

import pandas as pd

from .core import netCoin
from .network import netCreate
from .utils import createHTML

_MODES = {"d": "default", "p": "parallel", "f": "frame"}


def _is_graph_folder(graph) -> bool:
  return isinstance(graph, (str, os.PathLike)) and os.path.isfile(os.path.join(graph, "index.html"))


def _clean_graphs(graphs, names) -> dict:
  """The graphs as an ordered dict name -> graph, skipping unsupported objects."""
  if isinstance(graphs, dict):
    names = list(graphs.keys())
    graphs = list(graphs.values())
  if not graphs or not isinstance(graphs, (list, tuple)):
    raise ValueError("graphs: must be a non-empty list")

  if names is not None and len(names) != len(graphs):
    warnings.warn("names: length must match 'graphs' length")
    names = None

  out = {}
  for i, graph in enumerate(graphs):
    if not (isinstance(graph, netCoin) or _is_graph_folder(graph) or hasattr(graph, "_create")):
      warnings.warn("Not supported object found (%s)" % type(graph).__name__)
      continue
    name = str(names[i]) if names is not None else "graph" + str(len(out) + 1)
    out[name] = graph
  if not out:
    raise ValueError("graphs: no supported graph found")
  return out


def _write_graph(graph, directory: str, title: str, idx: int, names: list[str]):
  """One graph of a multigraph, in its own folder."""
  if _is_graph_folder(graph):
    copytree(str(graph), directory)
    return
  graph = deepcopy(graph)                      # the caller's object is left untouched
  graph.options["multigraph"] = {"idx": idx, "names": names}
  graph._create(directory, title)


def multiGraph(graphs: dict, directory: str):
  """The selector page and one folder per graph."""
  names = list(graphs)
  options = {"names": names}
  external = [_is_graph_folder(g) for g in graphs.values()]
  if any(external):
    options["external"] = external
  createHTML(directory, ["reset.css"], ["multigraph.js"], {"options": options})
  os.makedirs(os.path.join(directory, "data"), exist_ok=True)
  for i, (name, graph) in enumerate(graphs.items()):
    _write_graph(graph, os.path.join(directory, "data", str(i)), name, i, names)


def polyGraph(graphs: dict, mfrow: Sequence[int], directory: str):
  """Several graphs shown at the same time in a grid of ``mfrow`` rows x columns."""
  createHTML(directory, [], ["polygraph.js"], {"options": {"mfrow": list(mfrow), "n": len(graphs)}})
  if mfrow[0] * mfrow[1] > 2:                  # too small for the controls
    graphs = {k: _without_controls(g) for k, g in graphs.items()}
  multiGraph(graphs, os.path.join(directory, "multiGraph"))


def _without_controls(graph):
  if isinstance(graph, netCoin):
    graph = deepcopy(graph)
    graph.options.pop("controls", None)
  return graph


# ---------------------------------------------------------------- frames ----
def _chr(value) -> str:
  """``as.character`` of one value; missing values are blank."""
  if value is None or (isinstance(value, float) and value != value):
    return ""
  if isinstance(value, float) and value.is_integer():
    return str(int(value))
  return str(value)


def evolvingNetwork(graphs: dict, frame: int = 0, speed: float = 50, loop: bool = False,
                    lineplots: Sequence[str] | None = None) -> netCoin:
  """Merge networks that share nodes into one that evolves, one frame per graph.

  Links get a ``_frame_`` column. A node attribute that changes between frames
  becomes the text ``"v1|v2|..."`` that ``network.js`` splits into one value per
  frame (a blank where the node is absent).
  """
  if len(graphs) < 2:
    raise ValueError("Cannot make an evolving network with only one graph")
  if not all(isinstance(g, netCoin) for g in graphs.values()):
    raise ValueError("All graphs must be 'netCoin' objects")
  if not (isinstance(frame, (int, float)) and frame >= 0):
    warnings.warn("frame: must be an integer greater than 0")
    frame = 0
  if not (isinstance(speed, (int, float)) and 0 <= speed <= 100):
    warnings.warn("speed: must be numeric between 0 and 100")
    speed = 50

  frames = list(graphs)
  multi = list(graphs.values())
  first = multi[0]

  name = first.options["nodeName"]
  if any(g.options["nodeName"] != name for g in multi):
    raise ValueError("name: all graphs must have the same name")
  node_cols = list(first.nodes.columns)
  for g in multi:
    if g.links is None:
      raise ValueError("links: all graphs must have links")
    if set(g.nodes.columns) != set(node_cols):
      raise ValueError("nodes: all graphs must have the same node columns")
  link_cols = list(first.links.columns)
  for g in multi:
    if set(g.links.columns) != set(link_cols):
      raise ValueError("links: all graphs must have the same link columns")
  source, target = first.options["linkSource"], first.options["linkTarget"]
  for g in multi:
    if (g.options["linkSource"], g.options["linkTarget"]) != (source, target):
      raise ValueError("source and target: all graphs must have the same names")

  # nodes: the union, in order of appearance
  node_names = list(dict.fromkeys(n for g in multi for n in g.nodes[name]))
  position = [{n: i for i, n in enumerate(g.nodes[name])} for g in multi]
  nodes = {name: node_names}
  for col in node_cols:
    if col == name:
      continue
    merged = []
    for n in node_names:
      values = [g.nodes[col].iloc[pos[n]] if n in pos else None for g, pos in zip(multi, position)]
      values = [None if (isinstance(v, float) and v != v) else v for v in values]
      distinct = list(dict.fromkeys(values))
      merged.append(distinct[0] if len(distinct) == 1 else "|".join(_chr(v) for v in values))
    nodes[col] = merged
  nodes = pd.DataFrame(nodes)

  # links, with the frame they belong to
  parts = []
  for i, g in enumerate(multi):
    part = g.links[link_cols].copy()
    part["_frame_"] = i
    parts.append(part)
  links = pd.concat(parts, ignore_index=True)

  options = deepcopy(first.options)
  for item in ("main", "note", "repulsion", "distance", "zoom"):
    values = [g.options.get(item) for g in multi]
    if len(set(map(str, values))) != 1:
      options[item] = values
  options["frames"] = frames
  options["frame"] = frame
  options["speed"] = speed
  if loop is True:
    options["loop"] = True
  if lineplots:
    allowed = set(node_cols) | set(link_cols) | {"degree"}
    keys = [k for k in dict.fromkeys(map(str, lineplots)) if k in allowed]
    if keys:
      options["lineplotsKeys"] = keys

  trees = []
  for i, g in enumerate(multi):
    if g.tree is not None:
      part = g.tree.copy()
      part["_frame_"] = i
      trees.append(part)
  tree = pd.concat(trees, ignore_index=True) if trees else None
  return netCoin._from_parts(nodes, links, options, tree)


def frameGraph(graphs: dict, directory: str, frame=0, speed=50, loop=False, lineplots=None):
  netCreate(evolvingNetwork(graphs, frame, speed, loop, lineplots), directory)


# --------------------------------------------------------------- public ----
def multigraphCreate(graphs, names=None, mode="default", directory="MultiGraph", show=True, *,
                     mfrow=(1, 2), frame=0, speed=50, loop=False, lineplots=None) -> str:
  """Write several graphs into one folder and (optionally) open it.

  Port of ``multigraphCreate`` from the R package netCoin.

  Parameters
  ----------
  graphs : list or dict
      The graphs: :class:`netCoin` objects or folders that already hold a graph
      (with an ``index.html``). A dict gives the names.
  names : list of str, optional
      Name of each graph (``graph1``, ``graph2``... by default).
  mode : {"default", "parallel", "frame"}
      ``default``: a selector to move from one graph to another;
      ``parallel``: all of them at the same time, in a grid of ``mfrow``;
      ``frame``: the networks are merged into one that evolves, one frame per graph.
  directory : str
      Folder to write. An existing folder is only replaced if it is a graph
      made by this package.
  show : bool
      Open it in the browser.
  mfrow : (rows, columns)
      Grid of the parallel mode.
  frame, speed, loop, lineplots
      Frame mode: the first frame shown, the speed (0-100), whether it repeats,
      and the variables to draw as line plots (a node or link column, or
      ``"degree"``).

  Returns
  -------
  str
      The folder written.
  """
  if isinstance(mode, (list, tuple)):
    mode = mode[0]
  letter = str(mode)[:1].lower()
  if letter not in _MODES:
    raise ValueError("mode: must be 'default', 'parallel' or 'frame'")

  graphs = _clean_graphs(graphs, names)
  directory = str(directory)

  if letter == "p":
    if not (len(mfrow) == 2 and all(isinstance(x, int) and x > 0 for x in mfrow)):
      raise ValueError("mfrow: must be two positive integers")
    polyGraph(graphs, mfrow, directory)
  elif letter == "f":
    frameGraph(graphs, directory, frame, speed, loop, lineplots)
  else:
    multiGraph(graphs, directory)

  if show:
    webbrowser.open(os.path.join(directory, "index.html"))
  return directory
