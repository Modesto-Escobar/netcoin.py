"""Pie charts of coincidences (ports of rD3plot's ``pie_rd3`` and netCoin's ``pieCoin``)."""
from __future__ import annotations

import warnings
from typing import Sequence

import numpy as np
import pandas as pd

from .coincidence import Coin
from .webgraph import (WebGraph, check_cex, check_item_value, check_language, col2hex, columns,
                       is_color)

DEFAULT_COLORS = ("#000000", "#8dc7e6", "#ffffff", "#005587")
DEFAULT_LABELS = ("XY", "X", "", "Y")


def pie_from_coin(c: Coin) -> tuple[np.ndarray, np.ndarray]:
    """For every pair, its 2x2 table (a, b, c, d) and the expected pair (x, n - x).

    Port of ``piefromCoin``: both arrays are indexed ``[row, column, part]``.
    """
    a = c.matrix.to_numpy(dtype=float)
    n = float(c.n)
    diag = np.diag(a)
    k = a.shape[0]
    v = np.empty((k, k, 4))
    v[:, :, 0] = a
    v[:, :, 1] = diag[:, None] - a
    v[:, :, 2] = diag[None, :] - a
    v[:, :, 3] = n + a - diag[:, None] - diag[None, :]
    expected = np.outer(diag, diag) / n
    w = np.stack([expected, n - expected], axis=2)
    return v, w


class PieCoin(WebGraph):
    """A matrix of pie charts, one per pair of variables.

    Port of rD3plot's ``pie_rd3``. Build it with :func:`pie_coin`.

    Attributes
    ----------
    data : ndarray
        ``[row, column, part]`` counts of each pie.
    dataw : ndarray or None
        The expected counts, the same way.
    rownames, colnames : list of str
    nodes, links : DataFrame or None
        Optional information of the variables and of the pairs.
    options : dict
    """

    _engine = "pie.js"
    _needs_jspdf = False

    def __init__(
        self,
        v: np.ndarray,
        w: np.ndarray | None = None,
        labels: Sequence[str] | None = None,
        colors: Sequence[str] | None = None,
        nodes: pd.DataFrame | None = None,
        links: pd.DataFrame | None = None,
        name: str | None = None,
        source: str | None = None,
        target: str | None = None,
        lcolor=None,
        abline_x: int | None = None,
        abline_y: int | None = None,
        hide_upper: bool = False,
        main: str | None = None,
        note: str | None = None,
        show_legend: bool = True,
        help: str | None = None,
        help_on: bool = False,
        cex: float = 1,
        language: str = "en",
        rownames: Sequence[str] | None = None,
        colnames: Sequence[str] | None = None,
        directory: str | None = None,
    ):
        v = np.asarray(v, dtype=float)
        if v.ndim != 3:
            raise ValueError("v: incorrect data format")
        dims = v.shape
        self.data = v
        self.dim = list(dims)
        self.rownames = list(rownames) if rownames is not None else None
        self.colnames = list(colnames) if colnames is not None else None

        self.dataw = None
        if w is not None:
            w = np.asarray(w, dtype=float)
            if w.ndim != 3 or w.shape[:2] != dims[:2] or w.shape[2] != 2:
                warnings.warn("w: incorrect dimensions")
            else:
                self.dataw = w

        options: dict = {}
        self.nodes = None
        if nodes is not None:
            name = name or str(nodes.columns[0])
            nodes = nodes.copy()
            nodes[name] = nodes[name].astype(str)
            universe = list(dict.fromkeys((self.colnames or []) + (self.rownames or [])))
            matches = [n for n in nodes[name] if n in universe]
            if matches:
                options["nodeName"] = name
                nodes = nodes.drop_duplicates(name).set_index(name, drop=False)
                self.nodes = nodes.loc[list(dict.fromkeys(matches))].reset_index(drop=True)
            else:
                warnings.warn("nodes: no match with 'v'")

        self.links = None
        if links is not None:
            if dims[0] == 1 and dims[1] == 1:
                warnings.warn("links: not enough dimensions to apply link information")
            else:
                links = links.copy()
                source = source or str(links.columns[0])
                links[source] = links[source].astype(str)
                links = links[links[source].isin(self.rownames or [])]
                options["linkSource"] = source
                target = target or str(links.columns[1])
                links[target] = links[target].astype(str)
                links = links[links[target].isin(self.colnames or [])]
                options["linkTarget"] = target
                if links.empty:
                    warnings.warn("links: no row (source and target) matches the rownames and "
                                  "colnames of data")
                else:
                    self.links = links.reset_index(drop=True)

        if labels is not None:
            if len(labels) != dims[2]:
                warnings.warn("labels: incorrect length")
            else:
                options["labels"] = list(labels)
        if colors is not None and is_color(list(colors)):
            if len(colors) != dims[2]:
                warnings.warn("colors: incorrect length")
            else:
                options["colors"] = [col2hex(c) for c in colors]

        for key, value, axis in (("ablineX", abline_x, 0), ("ablineY", abline_y, 1)):
            if value is not None:
                if value < 1 or value >= dims[axis]:
                    warnings.warn(f"{key}: incorrect values")
                else:
                    options[key] = int(value)

        for key, value in (("main", main), ("note", note), ("help", help)):
            if value is not None:
                options[key] = value
        options["cex"] = check_cex(cex)
        options["language"] = check_language(language)
        for key, flag in (("hideUpper", hide_upper), ("showLegend", show_legend), ("helpOn", help_on)):
            if flag is True:
                options[key] = True
        self.options = options

        if self.links is not None:
            self.links, self.options = check_item_value(self.links, self.options, "linkColor",
                                                        lcolor, "lcolor")
        if directory:
            self._create(directory)

    def _data(self) -> dict:
        data = {"data": self.data, "dim": self.dim}
        if self.colnames is not None:
            data["colnames"] = self.colnames
        if self.rownames is not None:
            data["rownames"] = self.rownames
        if self.dataw is not None:
            data["dataw"] = self.dataw
        if self.nodes is not None:
            data["nodes"] = columns(self.nodes)
            data["nodenames"] = [str(c) for c in self.nodes.columns]
        if self.links is not None:
            data["links"] = columns(self.links)
            data["linknames"] = [str(c) for c in self.links.columns]
        data["options"] = self.options
        return data


def pie_coin(
    x: Coin,
    colors: Sequence[str] = DEFAULT_COLORS,
    nodes: pd.DataFrame | None = None,
    links: pd.DataFrame | None = None,
    name: str | None = None,
    lcolor=None,
    expected: bool = True,
    abline: int | None = None,
    main: str | None = None,
    note: str | None = None,
    show_legend: bool = True,
    help: str | None = None,
    help_on: bool = False,
    cex: float = 1,
    language: str = "en",
    directory: str | None = None,
) -> PieCoin:
    """A matrix of pie charts: every pair of variables, split by who has what.

    Port of ``pieCoin`` from the R package netCoin. Each pie shows the cases that
    have both characteristics, only the row one, only the column one and
    neither; with ``expected`` the expected coincidences are marked.

    Parameters
    ----------
    x : Coin
        See :func:`coin`.
    colors : four colours
        For both, row only, neither and column only (in the order of the pie).
    nodes, links : DataFrame, optional
        Information about the variables and the pairs (links from ``Source`` rows to
        ``Target`` columns).
    lcolor
        Colour of the links: a column of ``links``, a colour, or one value per link.
    expected : bool
        Mark the expected coincidences.
    abline : int, optional
        Draw a line after that many rows and columns.
    main, note, show_legend, help, help_on, cex, language, directory
        Page options.

    Returns
    -------
    PieCoin
    """
    if not isinstance(x, Coin):
        raise TypeError("x must be a Coin object (see coin())")
    v, w = pie_from_coin(x)
    v = v[:, :, [0, 1, 3, 2]]                      # both, X only, neither... as R orders them
    if expected:
        for part in range(2):
            np.fill_diagonal(w[:, :, part], np.nan)
    else:
        w = None
    return PieCoin(v, w, list(DEFAULT_LABELS), list(colors), nodes, links, name, "Source", "Target",
                   lcolor, abline, abline, True, main, note, show_legend, help, help_on, cex, language,
                   rownames=x.labels, colnames=x.labels, directory=directory)
