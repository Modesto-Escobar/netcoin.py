"""The page of a network: its JSON, its images and its files (port of rD3plot's ``networkJSON``/``netCreate``)."""
from __future__ import annotations

import os
import re
import shutil
import warnings
from copy import deepcopy

import pandas as pd

from .utils import createHTML, getLanguageScript, www
from .webgraph import columns, r_json


def getRawName(filepath: str) -> str:
  """The file name as hexadecimal text plus its extension (so any name is a safe one)."""
  filename = os.path.basename(filepath).split(".")
  ext = filename[-1]
  filename = ".".join(filename[:-1])
  return filename.encode("utf-8").hex() + "." + ext


def networkJSON(net) -> dict:
  """The data of the page of a network: what ``network.js`` reads from the ``data`` element."""
  nodes = net.nodes
  links = net.links
  options = deepcopy(net.options)
  name = [str(v) for v in nodes[options["nodeName"]]]
  index = {n: i for i, n in enumerate(name)}

  data = {"nodes": columns(nodes), "nodenames": [str(c) for c in nodes.columns]}
  if links is not None and len(links):
    out = links.copy()
    out[options["linkSource"]] = [index[str(v)] for v in out[options["linkSource"]]]
    out[options["linkTarget"]] = [index[str(v)] for v in out[options["linkTarget"]]]
    data["links"] = columns(out)
    data["linknames"] = [str(c) for c in out.columns]

  tree = getattr(net, "tree", None)
  if tree is not None and len(tree):
    targets = [str(v) for v in tree["Target"]]
    checkdup = targets
    if "_frame_" in tree.columns:
      checkdup = [t + str(f) for t, f in zip(targets, tree["_frame_"])]
    if len(set(checkdup)) == len(checkdup):
      indexed = tree.copy()
      indexed["Source"] = [index[str(v)] for v in tree["Source"]]
      indexed["Target"] = [index[str(v)] for v in tree["Target"]]
      data["tree"] = columns(indexed)
    else:
      warnings.warn("tree: there must be only one parent per node")

  layouts = getattr(net, "layouts", None)
  if layouts:
    data["layouts"] = {k: v.to_numpy().tolist() for k, v in layouts.items()}

  if options.get("statistics"):
    data["nodetypes"] = ["number" if pd.api.types.is_numeric_dtype(nodes[c]) else "string"
                         for c in nodes.columns]
    categories = {str(c): [str(x) for x in nodes[c].cat.categories]
                  for c in nodes.columns if isinstance(nodes[c].dtype, pd.CategoricalDtype)}
    if categories:
      options["categories"] = categories

  data["options"] = options
  return data


def _copy_image(filepath, image_dir) -> str:
  rawname = getRawName(filepath)
  if os.path.isfile(filepath):
    destination = os.path.join(image_dir, rawname)
    if not os.path.exists(destination):
      shutil.copyfile(filepath, destination)
  else:
    warnings.warn("missing image file (%s)" % filepath)
  return "images/" + rawname


def netWrapper(net, directory):
  """Copy the images, the background and the icons into ``directory`` and give the page data.

  The network itself is not modified.
  """
  net = deepcopy(net)
  options = net.options
  nodes = net.nodes
  image_dir = os.path.join(directory, "images")

  if options.get("imageItems"):
    os.makedirs(image_dir, exist_ok=True)
    items = options["imageItems"]
    if "imageNames" not in options:
      names = list(items)
      options["imageNames"] = names
      options["imageItems"] = items = [i + "_url" for i in names]
      regex = re.compile(r"\.[a-zA-Z0-9]+$")
      for item, original in zip(items, names):
        nodes[item] = nodes[original]
        nodes[original] = ["|".join(regex.sub("", os.path.basename(y)) for y in str(x).split("|"))
                           for x in nodes[original]]
    for item in options["imageItems"]:
      nodes[item] = ["|".join(_copy_image(f, image_dir) for f in str(x).split("|")) if pd.notna(x) else x
                     for x in nodes[item]]

  background = options.get("background")
  if background is not None and os.path.isfile(background):
    os.makedirs(image_dir, exist_ok=True)
    rawname = getRawName(background)
    shutil.copyfile(background, os.path.join(image_dir, rawname))
    options["background"] = 'url("images/%s")' % rawname

  if options.get("statistics"):
    source = os.path.join(www, "stats-icons")
    destination = os.path.join(directory, "stats-icons")
    os.makedirs(destination, exist_ok=True)
    for icon in os.listdir(source):
      shutil.copy(os.path.join(source, icon), destination)

  return networkJSON(net)


def netCreate(graph, directory, title=None):
  """Write the page of a network into ``directory``."""
  language = getLanguageScript(graph) or "en.js"

  styles = ["reset.css", "styles.css"]
  scripts = ["d3.min.js", "jspdf.min.js", "jszip.min.js", "iro.min.js"]
  if graph.options.get("statistics"):
    styles.append("statistics.css")
    scripts += ["stdlib.bundle.js", "statistics_" + language]
  scripts += ["images.js", "colorScales.js", "functions.js", language, "network.js"]

  createHTML(directory, styles, scripts, lambda: r_json(netWrapper(graph, directory)), title)
