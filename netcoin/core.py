"""The ``netCoin`` network object (port of netCoin's ``netCoin`` and rD3plot's ``network_rd3``)."""
from __future__ import annotations

import inspect
import math
import os
import warnings
import webbrowser
from copy import deepcopy
from tempfile import mkdtemp
from types import SimpleNamespace
from typing import Sequence

import numpy as np
import pandas as pd

from .dichotomize import _as_character
from .network import netCreate
from .utils import frequencyList, getByLanguage, nameList
from .webgraph import (capitalize, check_cex, check_default_color, check_item_value, check_language,
                       col2hex, column_option, get_shapes, is_color, is_shape)

# camelCase names of earlier versions (and of the R package) -> snake_case
_LEGACY = {
    "groupText": "group_text", "labelSize": "label_size", "imageNames": "image_names",
    "nodeBipolar": "node_bipolar", "nodeScaleLimits": "node_scale_limits",
    "nodeFilter": "node_filter", "degreeFilter": "degree_filter", "linkBipolar": "link_bipolar",
    "linkScaleLimits": "link_scale_limits", "linkFilter": "link_filter",
    "showCoordinates": "show_coordinates", "showArrows": "show_arrows", "showLegend": "show_legend",
    "showAxes": "show_axes", "axesLabels": "axes_labels", "nodesLabel": "nodes_label",
    "linksLabel": "links_label", "helpOn": "help_on", "roundedItems": "rounded_items",
    "defaultColor": "default_color", "nodeAttrNames": "node_attr_names",
    "linkAttrNames": "link_attr_names", "dir": "directory",
}


def _normalize_legacy(kwargs: dict, strict: bool = True) -> dict:
    """Turn camelCase keywords into the snake_case ones, with a deprecation warning.

    With ``strict`` (the constructor) any other keyword is an error; without it
    (``add``) the other keywords are kept as they are.
    """
    out: dict = {}
    for key, value in kwargs.items():
        if key in _LEGACY:
            new = _LEGACY[key]
            warnings.warn(f"'{key}' is deprecated, use '{new}'", DeprecationWarning, stacklevel=4)
        elif key in ("orderA", "orderD"):
            warnings.warn(f"'{key}' is deprecated, use 'sort' (and 'decreasing')",
                          DeprecationWarning, stacklevel=4)
            out["sort"] = value
            out["decreasing"] = key == "orderD"
            continue
        elif strict:
            raise TypeError(f"netCoin() got an unexpected keyword argument '{key}'")
        else:
            out[key] = value
            continue
        if new in out:
            raise TypeError(f"netCoin() got multiple values for '{new}'")
        out[new] = value
    return out


def _is_default(value, default) -> bool:
    try:
        return bool(value == default) if not isinstance(value, (list, np.ndarray, pd.DataFrame))             else False
    except (TypeError, ValueError):
        return False


def _as_frame(table, what: str):
    if table is None:
        return None
    frame = table.copy() if isinstance(table, pd.DataFrame) else pd.DataFrame(table)
    frame.columns = [str(c) for c in frame.columns]
    frame = frame.reset_index(drop=True)
    for col in frame.columns:
        if isinstance(frame[col].dtype, pd.CategoricalDtype):
            continue
    if frame.empty:
        warnings.warn(f"The {what} data frame is empty")
        return None
    return frame


def _make_unique(values: list[str]) -> list[str]:
    """R's ``make.unique``: the repeated names get ``.1``, ``.2``..."""
    seen: dict[str, int] = {}
    taken = set(values)
    out = []
    for v in values:
        if v not in seen:
            seen[v] = 0
            out.append(v)
        else:
            while True:
                seen[v] += 1
                candidate = f"{v}.{seen[v]}"
                if candidate not in taken:
                    break
            taken.add(candidate)
            out.append(candidate)
    return out


def name_by_language(name, language, nodes) -> str:
    """The column of names: the given one, or ``name``/``nombre``/``nom`` by language."""
    if name is None:
        name = getByLanguage(nameList, language)
    if nodes is not None:
        if name not in nodes.columns:
            warnings.warn(f"name: '{name}' column missing in nodes data frame")
        elif nodes[name].duplicated().any():
            warnings.warn(f"name: '{name}' column values must be unique")
    return name


def set_attr_by_value_key(label: str, item, items: pd.DataFrame | None):
    """A colour (or shape) given as ``{column: {value: colour}}``: one colour per row.

    Port of ``setAttrByValueKey``: each key must be a text column of ``items`` and its
    value a dict that maps the values of that column to colours; the result is a
    DataFrame with a column per key. Anything else is returned as it is.
    """
    if not isinstance(item, dict):
        return item
    checked = {}
    for key, mapping in item.items():
        if items is None or key not in items.columns or not (
                items[key].dtype == object or isinstance(items[key].dtype, pd.CategoricalDtype)):
            warnings.warn(f"{label}: the names in the list must match character columns of the items, "
                          f"but '{key}' doesn't")
        elif not isinstance(mapping, dict):
            warnings.warn(f"{label}: each item in the list must be a named character vector "
                          f"describing value-{label}, but '{key}' doesn't")
        else:
            checked[key] = [mapping.get(v) for v in items[key].astype(str)]
    return pd.DataFrame(checked) if checked else None


class netCoin:
    """An interactive network: nodes, links and the options of how to draw them.

    Port of ``netCoin`` (R). The tables are DataFrames: ``nodes`` and ``links`` (``None``
    when there are none), ``tree`` (``Source``, ``Target``) and ``options`` is the
    dict that the web page reads. Draw it with :meth:`plot` or ``directory=``.

    Parameters
    ----------
    nodes : DataFrame or dict, optional
        One row per node.
    links : DataFrame or dict, optional
        One row per link, with ``Source`` and ``Target`` (names of the nodes; if they are
        numbers, 1-based positions in ``nodes``). If there are no nodes they are taken
        from the links.
    tree : DataFrame or list of pairs, optional
        Parent -> child relations among the nodes.
    community : str, optional
        Detect communities with igraph (``"walktrap"``, ``"louvain"``, ``"fast_greedy"``,
        ``"edge_betweenness"``, ``"label_prop"``, ``"leading_eigen"``, ``"optimal"`` or
        ``"spinglass"``; two letters are enough) and put them in the node column
        ``community``.
    layout : str, matrix, DataFrame, list or dict, optional
        Coordinates: a name (``"fruchterman_reingold"``, ``"kamada_kawai"``, ``"circle"``,
        ``"grid"``, ``"star"``, ``"reingold_tilford"``, ``"sugiyama"``, ``"mds"``,
        ``"drl"``, ``"lgl"``, ``"davidson_harel"``, ``"random"``, ``"force_atlas2"``;
        two letters are enough), a matrix with two columns, or several of them in a dict
        (planes the page lets you switch between).
    name, label, group, group_text, label_size, size, color, shape, border, legend,
    sort, decreasing, ntext, info
        Node columns used for the name, label, groups, label size, size (or a number for
        all), colour, shape (``"Circle"``, ``"Square"``, ``"Diamond"``, ``"Triangle"``,
        ``"Cross"``, ``"Star"``, ``"Wye"``), border, legend, order, text of the tooltip
        and the info panel. ``color`` and ``shape`` can also be one value, one per
        node, or a dict ``{column: {value: colour}}``.
    image, image_names : list of str, optional
        Node columns with image files (and the names to show).
    centrality : list of str, optional
        Add the centralities ``"wdegree"``, ``"closeness"``, ``"betweenness"``,
        ``"eigen"`` as node columns (needs igraph).
    node_bipolar, link_bipolar : bool
        Treat negative and positive values as two poles of the colour scale.
    node_scale_limits, link_scale_limits : list of 2 or 3 numbers
        Limits of the numeric colour scale of ``color`` / ``lcolor``.
    node_filter, link_filter : str, callable, or sequence of bool
        Which nodes / links are shown. A text is a pandas expression over the columns
        (``"size > 10 and group == 'a'"``); a callable gets the table and returns the
        booleans.
    degree_filter : number or (min, max)
        Show only the nodes whose degree is in the range.
    lwidth, lweight, lcolor, ltext, intensity
        Link columns for the width, the weight (for layouts), the colour (also one
        value or one per link), the text and the intensity of the colour.
    repulsion, distance, zoom : number
        Behaviour of the force layout (0-100, 0-100, 0.1-10).
    fixed, limits, show_coordinates, show_arrows, show_axes, axes_labels
        Fix the positions, give the area (x1, y1, x2, y2), show the coordinates in the
        tables, draw arrows and axes (with their labels).
    main, note, help, help_on, scenarios, cex, background, default_color, language
        Titles, help text (``help_on`` shows it at the start), the number of scenarios,
        font size, background colour or image file, colour of the nodes, and
        ``"en"``/``"es"``/``"ca"``.
    show_legend, frequencies, statistics, nodes_label, links_label
        Show the legend, the frequency column, the statistics panel and the labels of the
        two tables.
    mode : list of str
        ``"network"``, ``"heatmap"`` or both.
    controls : list of int
        Which controls are shown: 1 sidebar, 2 selection buttons, 3 tables, 4 sliders and
        buttons, 5 export buttons, 6, 7, 8 legend (0 hides all, negatives remove).
    rounded_items : bool
        Rounded icons in the tables.
    directory : str, optional
        Write the page there.

    The camelCase names of earlier versions (``nodeFilter``, ``showArrows``,
    ``defaultColor``...) still work but are deprecated; ``orderA``/``orderD`` became
    ``sort`` and ``decreasing``.
    """

    def __init__(
        self,
        nodes=None,
        links=None,
        tree=None,
        community: str | None = None,
        layout=None,
        name: str | None = None,
        label: str | None = None,
        group: str | None = None,
        group_text: bool = False,
        label_size: str | None = None,
        size=None,
        color=None,
        shape=None,
        border: str | None = None,
        legend: str | None = None,
        sort: str | None = None,
        decreasing: bool = False,
        ntext: str | None = None,
        info: str | None = None,
        image: Sequence[str] | str | None = None,
        image_names: Sequence[str] | str | None = None,
        centrality=None,
        node_bipolar: bool = False,
        node_scale_limits: Sequence[float] | None = None,
        node_filter=None,
        degree_filter=None,
        lwidth: str | None = None,
        lweight: str | None = None,
        lcolor=None,
        ltext: str | None = None,
        intensity: str | None = None,
        link_bipolar: bool = False,
        link_scale_limits: Sequence[float] | None = None,
        link_filter=None,
        repulsion: float = 25,
        distance: float = 10,
        zoom: float = 1,
        fixed: bool | None = None,
        limits: Sequence[float] | None = None,
        main: str | None = None,
        note: str | None = None,
        show_coordinates: bool = False,
        show_arrows: bool = False,
        show_legend: bool = True,
        frequencies: bool = False,
        statistics: bool = False,
        show_axes: bool = False,
        axes_labels: Sequence[str] | None = None,
        nodes_label: str | None = None,
        links_label: str | None = None,
        scenarios: float | None = None,
        help: str | None = None,
        help_on: bool = False,
        mode: Sequence[str] | str = ("network", "heatmap"),
        rounded_items: bool = False,
        controls: Sequence[int] | int | None = tuple(range(1, 9)),
        cex: float = 1,
        background: str | None = None,
        default_color: str = "#1f77b4",
        language: str = "en",
        directory: str | None = None,
        *,
        source: str = "Source",
        target: str = "Target",
        node_attr_names: Sequence[str] | None = None,
        link_attr_names: Sequence[str] | None = None,
        **legacy,
    ):
        args = {k: v for k, v in locals().items() if k not in ("self", "legacy", "__class__")}
        defaults = inspect.signature(netCoin.__init__).parameters
        for key, value in _normalize_legacy(legacy).items():
            if key in args and not _is_default(args[key], defaults[key].default):
                raise TypeError(f"netCoin() got multiple values for '{key}'")
            args[key] = value
        self._args = {k: deepcopy(v) for k, v in args.items() if k != "directory"}
        self._build(SimpleNamespace(**args))
        if args["directory"]:
            self._create(args["directory"])

    # -------------------------------------------------------------- build ----
    def _build(self, a: SimpleNamespace):
        if a.links is None and a.nodes is None:
            raise ValueError("You must explicit a nodes or links data frame.")
        if isinstance(a.nodes, netCoin):
            raise TypeError("Using netCoin to change netCoin object attributes is deprecated. "
                            "Use the add() method instead.")

        nodes = _as_frame(a.nodes, "nodes")
        links = _as_frame(a.links, "links")
        if nodes is None and links is None:
            raise ValueError("You must explicit a nodes or links data frame.")
        name = name_by_language(a.name, a.language, nodes)

        color = set_attr_by_value_key("color", a.color, nodes)
        shape = set_attr_by_value_key("shape", a.shape, nodes)
        lcolor = set_attr_by_value_key("lcolor", a.lcolor, links)

        options: dict = {}

        # links and nodes
        source, target = a.source, a.target
        if links is not None:
            if source == "source":
                links = links.rename(columns={"source": "Source"})
                source = "Source"
            if target == "target":
                links = links.rename(columns={"target": "Target"})
                target = "Target"
            for column in (source, target):
                if column not in links.columns:
                    raise ValueError(f'"{column}" must be in links')
            options["linkSource"], options["linkTarget"] = source, target
            if nodes is not None and name in nodes.columns:
                names_for_index = _as_character(nodes[name])
            else:
                names_for_index = None
            for column in (source, target):
                if pd.api.types.is_numeric_dtype(links[column]) and names_for_index is not None:
                    links[column] = [names_for_index[int(v) - 1] if 0 < int(v) <= len(names_for_index)
                                     else str(v) for v in links[column]]
                else:
                    links[column] = _as_character(links[column])
            if nodes is None:
                union = list(dict.fromkeys(list(links[source]) + list(links[target])))
                nodes = pd.DataFrame({name: union})

        if name not in nodes.columns:
            raise ValueError(f'"{name}" must be in nodes')
        nodes[name] = _as_character(nodes[name])
        options["nodeName"] = name
        if nodes[name].duplicated().any():
            nodes[name] = _make_unique(list(nodes[name]))
            warnings.warn("name: must be unique in nodes")

        if links is not None:
            known = set(nodes[name])
            keep = links[source].isin(known) & links[target].isin(known) & (links[source] != links[target])
            links = links[keep].reset_index(drop=True)
            if links.empty:
                links = None
                warnings.warn(f"links: no row (source and target) matches the '{name}' column of the nodes")
        if links is None:
            options.pop("linkSource", None)
            options.pop("linkTarget", None)

        # graph options
        options["cex"] = check_cex(a.cex)
        distance, repulsion, zoom = a.distance, a.repulsion, a.zoom
        if not (isinstance(distance, (int, float)) and 0 <= distance <= 100):
            distance = 10
            warnings.warn("distance: must be numeric between 0 and 100")
        options["distance"] = distance
        if not (isinstance(repulsion, (int, float)) and 0 <= repulsion <= 100):
            repulsion = 25
            warnings.warn("repulsion: must be numeric between 0 and 100")
        options["repulsion"] = repulsion
        if not (isinstance(zoom, (int, float)) and 0.1 <= zoom <= 10):
            zoom = 1
            warnings.warn("zoom: must be numeric between 0.1 and 10")
        options["zoom"] = zoom

        if a.scenarios is not None:
            if isinstance(a.scenarios, (int, float, np.integer, np.floating)) and not isinstance(a.scenarios, bool):
                options["scenarios"] = int(a.scenarios) if float(a.scenarios).is_integer() else a.scenarios
            else:
                warnings.warn("scenarios: must be numeric")
        if a.limits is not None:
            if len(a.limits) != 4:
                warnings.warn("limits: must be a numeric list of length 4")
            else:
                options["limits"] = [float(x) for x in a.limits]

        for key, value in (("main", a.main), ("note", a.note), ("help", a.help)):
            if value is not None:
                options[key] = value
        if a.background is not None:
            if isinstance(a.background, str) and is_color(a.background):
                options["background"] = col2hex(a.background)
            elif isinstance(a.background, str) and os.path.isfile(a.background):
                options["background"] = a.background
            else:
                warnings.warn("background: must be a valid color or image file")

        options["language"] = check_language(a.language)
        default_color = check_default_color(a.default_color)
        if default_color is not None:
            options["defaultColor"] = default_color

        if a.help is None and options["language"] == "es":
            path = os.path.join(os.path.dirname(__file__), "www", "help_es.html")
            with open(path, encoding="utf-8") as handle:
                options["help"] = "".join(line.rstrip("\r\n") for line in handle)

        for key, flag in (("nodeBipolar", a.node_bipolar), ("linkBipolar", a.link_bipolar),
                          ("helpOn", a.help_on), ("frequencies", a.frequencies),
                          ("statistics", a.statistics)):
            if flag is True:
                options[key] = True
        if a.controls is not None:
            options["controls"] = [int(c) for c in ([a.controls] if np.isscalar(a.controls) else a.controls)]
        if a.mode is not None:
            letters = [str(m)[:1].lower() for m in ([a.mode] if isinstance(a.mode, str) else a.mode)]
            options["mode"] = letters if not isinstance(a.mode, str) else letters[0]
        for key, value in (("axesLabels", a.axes_labels), ("nodesLabel", a.nodes_label),
                           ("linksLabel", a.links_label)):
            if value is not None:
                options[key] = [str(v) for v in value] if not isinstance(value, str) else str(value)

        fixed = a.show_coordinates if a.fixed is None else a.fixed
        for key, flag in (("fixed", fixed), ("showCoordinates", a.show_coordinates),
                          ("showArrows", a.show_arrows), ("showAxes", a.show_axes),
                          ("roundedItems", a.rounded_items)):
            if flag is True:
                options[key] = True

        # the legend is control 8
        controls = options.get("controls")
        if a.show_legend:
            if controls is None:
                options["controls"] = [8]
            elif 0 not in controls and -8 not in controls and 8 not in controls:
                options["controls"] = controls + [8]
        elif controls is not None and 8 in controls:
            rest = [c for c in controls if c != 8]
            if rest:
                options["controls"] = rest
            else:
                del options["controls"]

        # node options
        if a.label is None:
            options["nodeLabel"] = name
        else:
            column_option(options, "nodeLabel", a.label)
        column_option(options, "nodeLabelSize", a.label_size)
        column_option(options, "nodeGroup", a.group)
        if a.group_text is True:
            options["groupText"] = True
        if isinstance(a.size, (int, float, np.integer, np.floating)) and not isinstance(a.size, bool):
            options["defaultNodeSize"] = a.size
        else:
            column_option(options, "nodeSize", a.size)
        column_option(options, "nodeLegend", a.legend)
        column_option(options, "nodeText", a.ntext)
        column_option(options, "nodeInfo", a.info)
        column_option(options, "nodeOrder", a.sort)
        if a.decreasing is True:
            options["decreasing"] = True
        column_option(options, "nodeBorder", a.border)

        if a.image is not None:
            images = [a.image] if isinstance(a.image, str) else list(a.image)
            if set(images) - set(nodes.columns):
                warnings.warn("image: names must match in nodes colnames.")
            else:
                options["imageItems"] = images
                if a.image_names is not None:
                    names = [a.image_names] if isinstance(a.image_names, str) else list(a.image_names)
                    if len(names) == len(images):
                        options["imageNames"] = names
                    else:
                        warnings.warn("imageNames: length must match with 'image' vector")

        # link options
        column_option(options, "linkWidth", a.lwidth)
        column_option(options, "linkWeight", a.lweight)
        column_option(options, "linkText", a.ltext)
        column_option(options, "linkIntensity", a.intensity)

        nodes, links = self._filters(nodes, links, name, source, target, a)

        # colours and shapes
        nodes, options = check_item_value(nodes, options, "nodeColor", color, "color")
        nodes, options = check_item_value(nodes, options, "nodeShape", shape, "shape",
                                          is_item=is_shape, auto_items=get_shapes, sanitize=capitalize)
        if links is not None:
            links, options = check_item_value(links, options, "linkColor", lcolor, "lcolor")

        self._scale_limits(nodes, links, options, a)

        self.nodes = nodes
        self.links = links
        self.options = options
        self.tree = self._check_tree(a.tree, nodes, name)
        self.layouts = None

        # column order
        if a.node_attr_names:
            first = [c for c in a.node_attr_names if c in nodes.columns]
            if name not in first:
                first = [name] + first
            self.nodes = nodes[first + [c for c in nodes.columns if c not in first]]
        if a.link_attr_names and links is not None:
            first = [c for c in a.link_attr_names if c in links.columns]
            self.links = links[first + [c for c in links.columns if c not in first]]

        if a.layout is not None:
            self._add_layout(a.layout, a.lweight)
        if a.community is not None:
            self._add_community(a.community)
        if a.centrality is not None:
            from .igraph_io import cal_centr
            columns = cal_centr(self, a.centrality)
            for col in columns.columns:
                if col not in ("nodes", "degree"):
                    self.nodes[col] = columns[col].to_numpy()

    @staticmethod
    def _filters(nodes, links, name, source, target, a):
        """``hidden`` columns from node_filter, link_filter and degree_filter."""
        if a.node_filter is None and a.link_filter is None and a.degree_filter is None:
            return nodes, links

        def evaluate(rule, frame):
            if callable(rule):
                rule = rule(frame)
            elif isinstance(rule, str):
                rule = frame.eval(rule)
            values = np.asarray(rule, dtype=bool)
            if values.shape != (len(frame),):
                raise ValueError("a filter must give one boolean per row")
            return values

        def hide_links():
            # In R this helper changes a local copy of the links and has no effect: node
            # filters never mark links as hidden in the object (the page hides the links of
            # hidden nodes itself). The same is done here so that both objects are equal.
            return None

        nodes = nodes.copy()
        nodes["hidden"] = False
        if links is not None:
            links = links.copy()
            links["hidden"] = False
        if a.node_filter is not None:
            nodes["hidden"] = ~evaluate(a.node_filter, nodes)
            hide_links()
        if links is not None and a.link_filter is not None:
            links["hidden"] = links["hidden"] | ~evaluate(a.link_filter, links)
        if a.degree_filter is not None:
            limits = [float(x) for x in np.atleast_1d(a.degree_filter)]
            if len(limits) == 1:
                limits.append(math.inf)
            if links is None:
                if min(limits) > 0:
                    nodes["hidden"] = True
            else:
                visible = links[~links["hidden"]]
                counts = pd.concat([visible[source], visible[target]]).value_counts()
                degree = nodes[name].map(counts).fillna(0)
                nodes["hidden"] = nodes["hidden"] | ~((degree >= limits[0]) & (degree <= limits[1]))
                hide_links()
        if not nodes["hidden"].any():
            nodes = nodes.drop(columns="hidden")
        if links is not None and not links["hidden"].any():
            links = links.drop(columns="hidden")
        return nodes, links

    @staticmethod
    def _scale_limits(nodes, links, options, a):
        for limits, key, frame, label in ((a.node_scale_limits, "nodeColor", nodes, "color"),
                                          (a.link_scale_limits, "linkColor", links, "lcolor")):
            if limits is None:
                continue
            values = sorted(set(float(x) for x in limits))
            column = options.get(key)
            if len(values) not in (2, 3):
                warnings.warn(f"{label.replace('lcolor', 'linkScaleLimits').replace('color', 'nodeScaleLimits')}"
                              ": must be a numeric list of length 2 or 3")
            elif column is None or frame is None or not pd.api.types.is_numeric_dtype(frame[column]):
                warnings.warn(f"parameter '{label}' must specify a numeric column")
            else:
                prefix = "nodeScaleLimits_" if key == "nodeColor" else "linkScaleLimits_"
                options[prefix + column] = values

    @staticmethod
    def _check_tree(tree, nodes, name):
        if tree is None:
            return None
        frame = tree.copy() if isinstance(tree, pd.DataFrame) else pd.DataFrame(list(tree))
        if frame.shape[1] < 2:
            raise ValueError("tree: needs two columns (Source and Target)")
        frame.columns = ["Source", "Target"] + [str(c) for c in frame.columns[2:]]
        frame["Source"] = _as_character(frame["Source"])
        frame["Target"] = _as_character(frame["Target"])
        known = set(nodes[name])
        frame = frame[frame["Source"].isin(known) & frame["Target"].isin(known)
                      & (frame["Source"] != frame["Target"])].reset_index(drop=True)
        if frame.empty:
            warnings.warn("tree: no row (Source and Target) matches the name column of the nodes")
            return None
        return frame[["Source", "Target"]]

    # ------------------------------------------------------ layout, groups ----
    def _layout_frame(self, value) -> pd.DataFrame | None:
        """A layout as a DataFrame of two columns, one row per node (None if invalid)."""
        n = len(self.nodes)
        frame = pd.DataFrame(value)
        if frame.shape != (n, 2) or not all(pd.api.types.is_numeric_dtype(t) for t in frame.dtypes):
            return None
        if isinstance(value, pd.DataFrame):
            index = [str(i) for i in value.index]
            names = list(self.nodes[self.options["nodeName"]])
            if set(index) == set(names) and index != names:
                frame = frame.set_axis(index, axis=0).loc[names]
        frame = frame.reset_index(drop=True)
        if all(isinstance(c, (int, np.integer)) for c in frame.columns):
            frame.columns = ["", ""]
        return frame

    def _add_layout(self, layout, lweight=None):
        from .igraph_io import to_igraph
        from .layouts import compute_layout, layout_code

        def compute(code):
            return pd.DataFrame(compute_layout(to_igraph(self), code, criteria=lweight),
                                columns=["", ""])

        bad = "layout: each layout must be a numeric matrix and have a pair of coordinates per node"
        if isinstance(layout, str):
            code = layout_code(layout)
            if code is not None:
                self._set_single_layout(compute(code))
            return
        if isinstance(layout, (list, tuple)) and layout and all(isinstance(x, str) for x in layout):
            layout = {x: x for x in layout}
        if isinstance(layout, dict):
            planes = {}
            for label, value in layout.items():
                if isinstance(value, str):
                    code = layout_code(value)
                    if code is None:
                        continue
                    planes[str(label)] = compute(code)
                else:
                    frame = self._layout_frame(value)
                    if frame is None:
                        warnings.warn(bad)
                        return
                    planes[str(label)] = frame
            if len(planes) == 1:
                self._set_single_layout(next(iter(planes.values())))
            elif planes:
                self.layouts = planes
                if "axesLabels" not in self.options:
                    self.options["layoutsAxesLabels"] = {k: [str(c) for c in v.columns]
                                                         for k, v in planes.items()}
            return
        frame = self._layout_frame(layout)
        if frame is None:
            warnings.warn(bad)
        else:
            self._set_single_layout(frame)

    def _set_single_layout(self, frame: pd.DataFrame):
        self.nodes["fx"] = frame.iloc[:, 0].to_numpy()
        self.nodes["fy"] = frame.iloc[:, 1].to_numpy()
        labels = [str(c) for c in frame.columns]
        if "axesLabels" not in self.options and any(labels):
            self.options["axesLabels"] = labels

    def _add_community(self, method):
        from .igraph_io import to_igraph
        from .layouts import community_code, compute_membership

        code = community_code(method)
        if code is None:
            return
        membership = compute_membership(to_igraph(self), code)
        width = len(str(max(membership)))
        self.nodes["community"] = [f"G.{m:0{width}d}" for m in membership]
        used = [self.options.get("nodeShape"), self.options.get("nodeColor")]
        if "community" not in used:
            self.options["nodeGroup"] = "community"

    # ---------------------------------------------------------- accessors ----
    def ncount(self) -> int:
        """Number of nodes."""
        return len(self.nodes)

    def lcount(self) -> int:
        """Number of links."""
        return 0 if self.links is None else len(self.links)

    def nodes_frame(self) -> pd.DataFrame:
        """A copy of the node table."""
        return self.nodes.copy()

    def links_frame(self) -> pd.DataFrame:
        """A copy of the link table (empty if there are no links)."""
        if self.links is None:
            return pd.DataFrame({"Source": [], "Target": []})
        return self.links.copy()

    def add(self, **kwargs):
        """A new netCoin object with some arguments changed or added.

        Port of ``addNetCoin``: every argument of the constructor can be given
        (``nodes``, ``links``, ``color``, ``main``...) and replaces the one the object was
        built with; the rest is kept. The original is not modified.
        """
        changes = _normalize_legacy(kwargs, strict=False)
        unknown = set(changes) - set(self._args) - {"directory"}
        if unknown:
            raise TypeError("unknown netCoin arguments: " + ", ".join(sorted(unknown)))
        arguments = deepcopy(self._args)
        arguments.update(changes)
        return netCoin(**arguments)

    @classmethod
    def _from_parts(cls, nodes: pd.DataFrame, links: pd.DataFrame | None, options: dict,
                    tree: pd.DataFrame | None = None, layouts: dict | None = None):
        """An object made of finished parts (no checks): what multigraphCreate builds."""
        net = cls.__new__(cls)
        net._args = {}
        net.nodes, net.links, net.options, net.tree, net.layouts = nodes, links, options, tree, layouts
        return net

    @classmethod
    def from_matrix(cls, matrix, link_attr_names: Sequence[str] | None = None, **kwargs):
        """A network from a matrix of links: one row per link, columns = source, target, ....

        The first two columns are the (1-based) positions in ``nodes`` of the source and
        target of each link (or names); the others are link attributes, named by
        ``link_attr_names`` (``Source``, ``Target``, ``attr2``... by default).
        """
        rows = [list(r) for r in matrix]
        if len({len(r) for r in rows}) != 1:
            raise ValueError("all rows in 'matrix' must have the same length")
        width = len(rows[0])
        if width < 2:
            raise ValueError("'matrix' must have at least two items")
        names = list(link_attr_names) if link_attr_names is not None and len(link_attr_names) >= 2 else None
        if names is None:
            names = ["Source", "Target"] + [f"attr{i}" for i in range(2, width)]
        names = names[:width]
        kwargs["source"], kwargs["target"] = names[:2]
        kwargs["links"] = pd.DataFrame(rows, columns=names)
        kwargs["link_attr_names"] = names
        return cls(**kwargs)

    @classmethod
    def fromMatrix(cls, *args, **kwargs):
        """Deprecated name of :meth:`from_matrix`."""
        warnings.warn("'fromMatrix' is deprecated, use 'from_matrix'", DeprecationWarning, stacklevel=2)
        if "linkAttrNames" in kwargs:
            kwargs["link_attr_names"] = kwargs.pop("linkAttrNames")
        return cls.from_matrix(*args, **kwargs)

    # ------------------------------------------------------------- output ----
    def _create(self, directory, title=None):
        """Write the web page of this graph (also used by multigraphCreate)."""
        netCreate(self, directory, title)

    def plot(self, directory: str | None = None):
        """Write the page (in a temporary folder by default) and open it in the browser."""
        directory = directory or mkdtemp()
        self._create(directory)
        webbrowser.open(os.path.join(directory, "index.html"))

    def __repr__(self) -> str:
        lines = ["<netCoin object>"]
        shown = [c for c in self.nodes.columns if c != self.options.get("nodeText")]
        lines.append(f"Nodes: {len(self.nodes)}  ({', '.join(map(str, shown))})")
        if self.links is not None:
            lines.append(f"Links: {len(self.links)}  ({', '.join(map(str, self.links.columns))})")
        else:
            lines.append("Links: 0")
        if {"fx", "fy"} <= set(self.nodes.columns):
            planes = list(self.layouts) if self.layouts else ["layout"]
            lines.append("Planes: " + ", ".join(planes))
        return "\n".join(lines)

    def __str__(self) -> str:
        output = []
        if "main" in self.options:
            output.append("Title: " + str(self.options["main"]))
            output.append("")
        for frame, count in ((self.nodes, self.ncount()), (self.links, self.lcount())):
            if count:
                output.append(frame.head(6).to_string(index=False))
                if count > 6:
                    output.append("...")
                output.append("")
        if "note" in self.options:
            output.append(str(self.options["note"]))
            output.append("")
        return "\n".join(output)

    def summary(self):
        """Print the number of nodes and links and how the frequency and link width are distributed."""
        lines = [f"{self.ncount()} nodes and {self.lcount()} links."]
        stats = ["Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max"]

        def describe(values):
            q = values.describe()
            return [f"{q[k]:.2f}" for k in ("min", "25%", "50%", "mean", "75%", "max")]

        row = "{:>10}" * len(stats)
        freq = [c for c in frequencyList if c in self.nodes.columns]
        if len(freq) == 1:
            lines += [f"{freq[0]} distribution of nodes:", row.format(*stats),
                      row.format(*describe(pd.to_numeric(self.nodes[freq[0]], errors="coerce")))]
        width = self.options.get("linkWidth")
        if width is None and self.links is not None and self.links.shape[1] > 2:
            width = self.links.columns[2]
        if width is not None and self.links is not None:
            lines += [f"{width}'s distribution:", row.format(*stats),
                      row.format(*describe(pd.to_numeric(self.links[width], errors="coerce")))]
        print("\n".join(lines))
