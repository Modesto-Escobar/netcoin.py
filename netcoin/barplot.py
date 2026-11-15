"""Bar charts of coincidences (ports of rD3plot's ``barplot_rd3`` and netCoin's ``barCoin``)."""
from __future__ import annotations

import warnings
from typing import Sequence

import numpy as np
import pandas as pd

from .coincidence import as_nodes, coin
from .dichotomize import _collation_key, dichotomize
from .similarity import _check_level, edge_list
from .utils import frequencyList, getByLanguage, labelList, nameList
from .webgraph import (WebGraph, check_cex, check_default_color, check_language, column_option)


class BarCoin(WebGraph):
    """A bar chart: for each event, the bars of its coincidences with the others.

    Port of rD3plot's ``barplot_rd3``. Build it with :func:`bar_coin`.

    Attributes
    ----------
    nodes : DataFrame
        One row per event, with its name and number of incidences.
    links : DataFrame
        ``Source``, ``Target`` and the statistics of every pair (``coincidences``,
        ``expected``, confidence limits, p-value...).
    options : dict
        What the page reads (column names, colours, language...).
    """

    _engine = "barplot.js"

    def __init__(
        self,
        events: pd.DataFrame,
        links: pd.DataFrame,
        name: str | None = None,
        select: str | None = None,
        source: str | None = None,
        target: str | None = None,
        label: str | None = None,
        text: str | None = None,
        color: str | None = None,
        incidences: str | None = None,
        coincidences: str | None = None,
        expected: str | None = None,
        confidence: str | Sequence[str] | None = None,
        level: float = 0.95,
        significance: str | None = None,
        sort: str | None = None,
        decreasing: bool = False,
        scalebar: bool = False,
        default_color: str = "#1f77b4",
        note: str | None = None,
        cex: float = 1,
        language: str = "en",
        directory: str | None = None,
    ):
        events = events.copy()
        links = links.copy()
        name = name or str(events.columns[0])
        events[name] = events[name].astype(str)
        incidences = incidences or str(events.columns[1])
        source = source or str(links.columns[0])
        target = target or str(links.columns[1])
        coincidences = coincidences or str(links.columns[2])

        missing = (set(links[source].astype(str)) | set(links[target].astype(str))) - set(events[name])
        if missing:
            raise ValueError(f"{len(missing)} node link(s) not defined in '{name}' column "
                             "from events data frame.")
        rest = [c for c in links.columns if c not in (source, target)]
        links = pd.concat([links[[source]].astype(str).set_axis(["Source"], axis=1),
                           links[[target]].astype(str).set_axis(["Target"], axis=1),
                           links[rest]], axis=1)

        options = {
            "name": name, "incidences": incidences, "coincidences": coincidences,
            "defaultColor": check_default_color(default_color[0] if isinstance(default_color, (list, tuple))
                                                else default_color),
            "cex": check_cex(cex), "language": check_language(language),
        }
        if options["defaultColor"] is None:
            del options["defaultColor"]
        for key, value in (("note", note), ("label", label), ("text", text), ("color", color),
                           ("select", select)):
            if value is not None:
                options[key] = value
        options["rev"] = int(not decreasing)
        if sort is not None:
            if sort in events.columns:
                options["order"] = sort
                options["rev"] = int(pd.api.types.is_numeric_dtype(events[sort])) ^ int(bool(decreasing))
            else:
                warnings.warn("sort: must be a 'events' column")
        if scalebar:
            options["scalebar"] = True
        if expected is not None:
            options["expected"] = expected
            if confidence is not None:
                options["level"] = level
                options["confidence"] = confidence if isinstance(confidence, str) else list(confidence)
        if significance is not None:
            options["significance"] = significance

        self.nodes = events
        self.links = links
        self.options = options
        if directory:
            self._create(directory)

    def _data(self) -> dict:
        from .webgraph import records
        return {"nodes": records(self.nodes), "links": records(self.links), "options": self.options}

    def summary(self):
        """Print the number of nodes and links."""
        print(f"{len(self.nodes)} nodes and {len(self.links)} links.")


def bar_coin(
    data: pd.DataFrame,
    variables: Sequence[str] | None = None,
    commonlabel: str | Sequence[str] | None = None,
    dichotomies: str | Sequence[str] = "_all",
    value_dicho=1,
    weight: str | Sequence[float] | None = None,
    subsample: bool = False,
    sort: str | None = None,
    decreasing: bool = True,
    nodes: pd.DataFrame | None = None,
    name: str | None = None,
    select: str | None = None,
    scalebar: bool = False,
    note: str | None = None,
    label: str | None = None,
    text: str | None = None,
    color: str | None = None,
    default_color: str = "#1f77b4",
    expected: bool = False,
    confidence: bool = False,
    level: float = 0.95,
    significance: bool = False,
    minimum: float = 1,
    maximum: float | None = None,
    percentages: bool = False,
    criteria: str = "Z",
    bonferroni: bool = False,
    support: float = 1,
    min_link: float = -np.inf,
    max_link: float = 1,
    language: str = "en",
    cex: float = 1.0,
    directory: str | None = None,
) -> BarCoin:
    """Bar chart of the coincidences of every variable (category) with the others.

    Port of ``barCoin`` from the R package netCoin.

    Parameters
    ----------
    data : DataFrame
        The survey or table of cases.
    variables : list of str, optional
        Columns to use (all by default).
    commonlabel : str or list, optional
        Variables whose name is kept in the label of their categories
        (``"variable:category"``); for the rest only the category is shown.
    dichotomies : ``"_all"``, ``"_none"`` or list of str
        Which variables are dichotomies (``value_dicho`` against the rest, e.g.
        0/1 columns) and which are categorical (one bar per category).
    value_dicho
        The value that counts as 1 in the dichotomies.
    weight : str or sequence, optional
        Weight of each case or the column that holds it.
    subsample : bool
        Use only the cases with at least one event.
    sort, decreasing
        Column of the nodes to order the events by.
    nodes : DataFrame, optional
        Your own table of nodes; it restricts and orders the events and adds
        its columns (``label``, ``text``, ``color`` can refer to them).
    name, select, scalebar, note, label, text, color, default_color, cex, language
        As in the page options.
    expected, confidence, significance : bool
        Draw the expected coincidences, their confidence limits (at ``level``)
        and the p-value of the link (``criteria``: ``"Z"`` or ``"hyp"``).
    minimum, maximum : float
        Frequency range of the events.
    percentages : bool
        Show percentages of the cases instead of counts.
    bonferroni, support, min_link, max_link
        Filters of the links, as in :func:`edge_list` (``support`` is the least
        number of coincidences and ``max_link`` bounds the p-value).
    directory : str, optional
        Write the page there.

    Returns
    -------
    BarCoin
    """
    language = language[0] if isinstance(language, (list, tuple)) else language
    name = name or getByLanguage(nameList, language)
    if nodes is not None:
        if name not in nodes.columns:
            warnings.warn(f"name: '{name}' column missing in nodes data frame")
        elif nodes[name].duplicated().any():
            warnings.warn(f"name: '{name}' column values must be unique")
    variables = list(data.columns) if variables is None else list(variables)

    criteria = "hyp" if str(criteria) == "hyp" else "Z"
    procedures = ["Frequencies", "Expected"] + (["Confidence"] if confidence else []) + [criteria]
    criteria_col = {"Z": "p(Z)", "hyp": "p(Fisher)"}[criteria]

    # which variables are dichotomies
    if isinstance(dichotomies, str):
        dichotomies = [dichotomies]
    dichotomies = list(dichotomies)
    if "_all" in dichotomies:
        dichotomies = list(variables)
    if "_none" in dichotomies:
        all_vars, dichotomies = list(variables), []
    else:
        all_vars = list(dict.fromkeys(dichotomies + variables))

    # weights and missing data
    weights = None
    if isinstance(weight, str):
        weights_column = weight
        all_vars = [v for v in all_vars if v != weights_column]
        variables = [v for v in variables if v != weights_column]
        dichotomies = [v for v in dichotomies if v != weights_column]
        complete = data[all_vars].notna().all(axis=1).to_numpy().copy()
        weights = data[weights_column].to_numpy(dtype=float)[complete]
        data = data.loc[complete, all_vars]
    else:
        complete = data[all_vars].notna().all(axis=1).to_numpy().copy()
        if weight is not None:
            weight = np.asarray(weight, dtype=float)
            if len(weight) != len(data):
                raise ValueError("Weights have not the correct dimensions")
            complete &= ~np.isnan(weight)
            weights = weight[complete]
        data = data.loc[complete, all_vars]

    incidences = None
    dichos = None
    if dichotomies:
        dichos = pd.DataFrame({v: (data[v] == value_dicho).astype(int) for v in dichotomies},
                              index=data.index)
        variables = [v for v in variables if v not in dichotomies]
    if variables:
        incidences = dichotomize(data, variables, sep="", min_freq=minimum, length=0, values=None,
                                 sparse=False, add=False, sort=True)
        if dichos is not None:
            incidences = pd.concat([dichos, incidences], axis=1)
    elif dichos is not None:
        incidences = dichos
    if incidences is None:
        raise ValueError("There are no qualitative variables. Try netCorr.")

    # keep only the events of the nodes table
    if nodes is not None:
        wanted = [str(v) for v in nodes[name]]
        missing = [v for v in wanted if v not in incidences.columns]
        if missing:
            warnings.warn(", ".join(missing) + " is/are not present among incidences.")
        incidences = incidences[[v for v in wanted if v in incidences.columns]]

    incidences = incidences.dropna()
    cooc = coin(incidences, minimum, maximum, sort=True, decreasing=True, weight=weights,
                subsample=subsample)

    # the nodes
    events = as_nodes(cooc, frequency=not percentages, percentages=percentages, language=language)
    events = events.rename(columns={events.columns[0]: name})
    events = events.rename(columns={events.columns[1]: "incidences"})
    if nodes is not None:
        user = nodes.copy()
        user[name] = user[name].astype(str)
        extra = user[[c for c in user.columns if c not in frequencyList]]
        events = events.merge(extra, on=name, how="left", sort=False)
        events = events.iloc[sorted(range(len(events)), key=lambda i: _collation_key(events[name].iloc[i]))
                             ].reset_index(drop=True)
    elif commonlabel is not None:
        label = getByLanguage(labelList, language)
        keep = {commonlabel} if isinstance(commonlabel, str) else set(commonlabel)
        full = events[name].astype(str)

        def short(text):
            colon = text.find(":")
            return text if colon < 0 else text[colon + 1:]

        prefix = full.map(lambda t: t[:t.find(":")] if ":" in t else "")
        events[label] = [t if p in keep else short(t) for t, p in zip(full, prefix)]

    # the links
    level = _check_level(level)
    links = edge_list(cooc, procedures, criteria, level, bonferroni, min_link, max_link, support,
                      directed=False, diagonal=False, sort=None)

    expected_option = "expected" if (expected or confidence) else None
    confidence_option = None
    if confidence:
        if significance:
            confidence_option = "confidence"
            links = links.drop(columns=["conf.L", "conf.U"])
        else:
            confidence_option = ["conf.L", "conf.U"]
            links = links.drop(columns=["confidence"])
    significance_option = None
    if significance:
        significance_option = criteria_col
    else:
        links = links.drop(columns=[criteria_col])

    if percentages:
        scaled = [c for c in ("coincidences", "expected", "confidence", "conf.L", "conf.U")
                  if c in links.columns]
        links[scaled] = links[scaled] / cooc.n * 100

    return BarCoin(
        events, links, name=name, select=select, source="Source", target="Target", label=label,
        text=text, color=color, incidences="incidences", coincidences="coincidences",
        expected=expected_option, confidence=confidence_option, level=level,
        significance=significance_option, sort=sort, decreasing=decreasing, scalebar=scalebar,
        default_color=default_color, note=note, cex=cex, language=language, directory=directory,
    )
