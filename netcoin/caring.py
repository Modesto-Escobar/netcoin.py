"""``caring``: a page with several networks of one survey (port of the ``caring`` of netCoin).

``caring`` reads a data file and, following a dictionary of arguments (the same the graphical
interface of netCoin sends to R), builds the networks asked for - coincidences of the categories
(:func:`sur_coin`), scattergram (:func:`sur_scat`), log-linear models (:func:`log_coin`), incidences
(:func:`all_net`) and regressions (:func:`glm_coin`) - and puts them in one page with
:func:`multigraphCreate`.
"""
from __future__ import annotations

import inspect
import os
import re
import tempfile
import warnings
from typing import Mapping

import numpy as np
import pandas as pd

# names of the arguments in R -> in python
_RENAMES = {
    "valueDicho": "value_dicho", "minL": "min_link", "maxL": "max_link", "sortL": "sort_links",
    "decreasingL": "decreasing_links", "Bonferroni": "bonferroni", "dir": "directory",
    "vPatterns": "v_patterns", "maxN": "max_n", "scaleAxes": "scale_axes", "critFactors": "crit_factors",
    "alignTo": "align_to", "clusterOn": "cluster_on", "scaleClusters": "scale_clusters",
    "clusterIndex": "cluster_index", "sortClusters": "sort_clusters", "noFirstCat": "no_first_cat",
    "percentage": "percentage", "percentages": "percentages",
}
_SINGLE = ("sur_coin", "sur_scat", "log_coin", "all_net", "glm_coin", "net_corr")


def _python_name(key: str) -> str:
    from .core import _LEGACY

    if key in _RENAMES:
        return _RENAMES[key]
    if key in _LEGACY:
        return _LEGACY[key]
    return re.sub(r"(?<=[a-z0-9])([A-Z])", lambda m: "_" + m.group(1).lower(), key)


def _accepted(function) -> set[str]:
    """The names a function takes, its ``**netcoin_args`` meaning those of :class:`netCoin`."""
    from .core import netCoin

    names = set(inspect.signature(function).parameters)
    if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in inspect.signature(function).parameters.values()):
        names |= set(inspect.signature(netCoin.__init__).parameters)
    names.discard("self")
    return names


def _keep(arguments: Mapping, function) -> dict:
    """The arguments (named as python names them) that ``function`` takes."""
    ok = _accepted(function)
    return {k: v for k, v in arguments.items() if k in ok and k not in ("data", "kwargs")}


# ----------------------------------------------------------------- reading ----
def caring_read_file(filepath: str) -> pd.DataFrame:
    """Read a ``.sav``, ``.dta``, ``.xlsx`` or delimited text file as a table of cases.

    Text columns become categories, and those that hold only numbers become numeric. A Stata file
    takes the label of each variable as its name and labelled values become categories.
    """
    ext = os.path.splitext(filepath)[1].lower().lstrip(".")
    if ext == "dta":
        reader = pd.io.stata.StataReader(filepath)
        labels = reader.variable_labels()
        data = reader.read(convert_categoricals=True)
        names, seen = [], {}
        for column in data.columns:
            label = labels.get(column) or column
            seen[label] = seen.get(label, -1) + 1
            names.append(label if seen[label] == 0 else f"{label}_{seen[label]}")
        data.columns = names
    elif ext == "sav":
        try:
            import pyreadstat
        except ImportError as exc:
            raise ImportError("reading .sav files needs pyreadstat: pip install pyreadstat") from exc
        data, meta = pyreadstat.read_sav(filepath, apply_value_formats=True)
        labels = dict(zip(meta.column_names, meta.column_labels))
        names, seen = [], {}
        for column in data.columns:
            label = labels.get(column) or column
            seen[label] = seen.get(label, -1) + 1
            names.append(label if seen[label] == 0 else f"{label}_{seen[label]}")
        data.columns = names
    elif ext in ("xlsx", "xlsm"):
        data = pd.read_excel(filepath)
    else:
        for encoding in ("utf-8", "latin-1"):
            try:
                data = pd.read_csv(filepath, sep=None, engine="python", encoding=encoding)
                break
            except UnicodeDecodeError:
                continue
    data.columns = [str(c).replace('"', "'") for c in data.columns]
    for column in data.columns:
        series = data[column]
        if series.dtype == object or isinstance(series.dtype, pd.CategoricalDtype):
            numbers = pd.to_numeric(series.astype(object), errors="coerce")
            if numbers.notna().sum() == series.notna().sum() and series.notna().any():
                data[column] = numbers
            else:
                data[column] = series.astype("category")
    return data


# ---------------------------------------------------------------- the graphs ----
def _values(arguments: Mapping, key: str):
    value = arguments.get(key)
    return value if value is None or isinstance(value, (list, tuple, np.ndarray)) else [value]


def caring_create_graphs(data: pd.DataFrame, arguments: Mapping, show: bool = False) -> str:
    """The networks of ``arguments["plot"]`` in one page; returns the folder.

    ``arguments`` holds the names of R (``variables``, ``weight``, ``maxL``, ``criteria``,
    ``dichotomies``, ``valueDicho``, ``metric``, ``exogenous``, ``chaine``, ``family``, ``nclustersMin``...) or
    python's; ``plot`` lists, one after the other, the kind of network (``"surCoin"``, ``"surScat"``,
    ``"logCoin2"`` to ``"logCoin5"``, ``"allNet"``, ``"glmCoin"``) and the title it gets in the page;
    ``dir`` (or ``directory``) is the folder to write.
    """
    from .allnet import all_net
    from .glmcoin import glm_coin
    from .logcoin import log_coin
    from .multigraph import multigraphCreate
    from .netcorr import net_corr
    from .surcoin import sur_coin
    from .surscat import sur_scat

    args = {k: v for k, v in dict(arguments).items()}
    initial = list(args.get("variables") or data.columns)
    args["variables"] = list(initial)
    plots = list(_values(args, "plot") or [])
    directory = args.pop("dir", None) or args.pop("directory", None)

    weight = _values(args, "weight")
    if weight:
        weight = weight[0]
        if weight not in data.columns or not pd.api.types.is_numeric_dtype(data[weight]):
            weight = None
        args["weight"] = weight
    else:
        args.pop("weight", None)

    pmax = None
    max_l = args.get("maxL", args.get("max_link"))
    if max_l is None or not np.isfinite(max_l):
        if args.get("criteria") in (None, "Z", "hyp"):
            args["maxL"] = 0.5
            args.pop("max_link", None)
    else:
        pmax = max_l

    dichotomies = _values(args, "dichotomies")
    if dichotomies is None:
        found = [c for c in args["variables"]
                 if pd.api.types.is_numeric_dtype(data[c]) and set(data[c].dropna().unique()) <= {0, 1}]
        if found:
            args["dichotomies"], args["valueDicho"] = found, 1
    if args.get("dichotomies") is not None and args.get("metric") is not None:
        args["metric"] = [m for m in _values(args, "metric") if m not in _values(args, "dichotomies")]

    args["frequency"], args["percentages"], args["percentage"] = True, True, True
    args["color"], args["data"] = "variable", data
    names = {}
    translated = {_python_name(k): v for k, v in args.items()}
    translated.pop("plot", None)

    def title_of(kind):
        return plots[plots.index(kind) + 1] if kind in plots and plots.index(kind) + 1 < len(plots) else kind

    graphs = {}
    metric = _values(args, "metric")

    if "surCoin" in plots:
        if metric and not set(translated["variables"]) - set(metric):
            corr_args = _keep(translated, net_corr)
            corr_args.pop("frequency", None)
            corr_args.pop("variables", None)
            if corr_args.get("weight") is not None:
                corr_args["weight"] = data[corr_args["weight"]].to_numpy(dtype=float)
            if str(corr_args.get("criteria", "")).lower() in ("z", "hyp"):
                corr_args["criteria"] = "p"
            graphs[title_of("surCoin")] = net_corr(data[metric], **corr_args)
        else:
            coin_args = _keep(translated, sur_coin)
            coin_args.pop("data", None)
            graphs[title_of("surCoin")] = sur_coin(data, **coin_args)

    if "surScat" in plots:
        scat = _keep(translated, sur_scat)
        scat.pop("data", None)
        if args.get("alignedAxes"):
            scat["suffix"] = "*"
        if scat.get("patterns") is not None:
            scat["v_patterns"] = scat["patterns"] if not isinstance(scat["patterns"], bool) else scat.get("v_patterns", initial)
            scat["patterns"] = True
        low, high, default = args.get("nclustersMin"), args.get("nclustersMax"), args.get("nclustersDef")
        auto = bool(args.get("nclustersAuto")) or str(default).lower() == "auto"
        if any(x is not None for x in (low, high, default)) or auto:
            numbers = None
            if low is not None and high is not None:
                numbers = list(range(int(low), int(high) + 1))
            elif low is not None:
                numbers = [int(low)]
            elif high is not None:
                numbers = [int(high)]
            if auto:
                numbers = ["auto"] + (numbers or [])
            elif default is not None:
                numbers = ([n for n in numbers if n != int(default)] + [int(default)]) if numbers else [int(default)]
            scat["nclusters"] = numbers
        if scat.get("type") == "":
            scat.pop("type")
        for key in ("degree_filter", "layout"):
            scat.pop(key, None)
        for key in ("clusters", "layouts"):
            if scat.get(key) and "all" in _values(scat, key):
                scat[key] = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*not available yet.*")
            graphs[title_of("surScat")] = sur_scat(data, **scat)

    for kind in [p for p in ("logCoin2", "logCoin3", "logCoin4", "logCoin5") if p in plots]:
        log_args = _keep(translated, log_coin)
        log_args.pop("data", None)
        frame = data.copy()
        dich = _values(args, "dichotomies")
        if dich:
            values = _values(args, "valueDicho")
            no_first = []
            for i, col in enumerate(dich):
                value = values[i if len(values) > i else 0]
                text = frame[col].astype(str)
                levels = [str(value)] + [v for v in dict.fromkeys(text) if v != str(value)]
                frame[col] = pd.Categorical(text, categories=levels)
                no_first.append(col)
            log_args["no_first_cat"] = no_first
        if pmax is not None:
            log_args["pmax"] = pmax
        log_args["color"] = "var"
        log_args["order"] = int(kind.replace("logCoin", ""))
        if metric:
            log_args["variables"] = [v for v in log_args.get("variables", initial) if v not in metric]
        graphs[title_of(kind)] = log_coin(frame, **log_args)

    if "allNet" in plots:
        inc_args = _keep(translated, all_net)
        incidences = data[initial].copy()
        if inc_args.get("color") not in incidences.columns:
            inc_args.pop("color", None)
        for column in incidences.columns:
            series = incidences[column]
            if not (pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series)):
                raise ValueError(f"allNet needs incidences: '{column}' is text, dichotomize it first or leave it out of variables")
            incidences[column] = (series != 0).astype(float).where(series.notna())
        inc_args.pop("incidences", None)
        graphs[title_of("allNet")] = all_net(incidences, **inc_args)

    if "glmCoin" in plots:
        graphs[title_of("glmCoin")] = _glm_graph(data, args, translated, initial, pmax, glm_coin)

    if not graphs:
        raise ValueError("plot must name at least one network: surCoin, surScat, logCoin2-5, allNet or glmCoin")
    folder = directory or tempfile.mkdtemp(prefix="netcoin_caring_")
    multigraphCreate(graphs, names=list(graphs), directory=folder, show=show)
    return folder


def glm_equations(data: pd.DataFrame, variables: list[str], exogenous: list[str], chaine: list, family: list,
                  metric=None, dichotomies=None, value_dicho=None) -> tuple[pd.DataFrame, str]:
    """The data and the equations of the regressions ``caring`` asks of ``glm_coin``.

    Every dichotomy becomes the 0/1 variable ``"variable:value"``, every other qualitative variable
    that is not exogenous becomes one 0/1 variable per category, and every variable is explained by
    those whose level in ``chaine`` is lower (the exogenous ones have level 0 and explain, never
    are explained).
    """
    metric = list(metric or [])
    dichotomies = list(dichotomies or [])
    value_dicho = list(value_dicho or [])
    frame = data[variables].copy()
    final = list(variables)
    for i, column in enumerate(dichotomies):
        value = value_dicho[i if len(value_dicho) > i else 0]
        new = str(column) + ":" + str(value)
        new = new.replace(" ", "")
        frame[new] = np.where(frame[column].isna(), np.nan, (frame[column].astype(str) == str(value)).astype(float))
        frame = frame.drop(columns=column)
        final[variables.index(column)] = new
    levels = list(chaine)
    families = list(family)
    text_variables = [v for v in variables if v not in set(metric) | set(dichotomies) | set(exogenous)]
    for variable in text_variables:
        position = variables.index(variable)
        for value in dict.fromkeys(frame[variable].tolist()):
            new = (str(variable) + ":" + str(value)).replace(" ", "")
            frame[new] = np.where(frame[variable].isna(), np.nan, (frame[variable] == value).astype(float))
            levels.append(chaine[position])
            families.append(family[position])
            final.append(new)
        keep = [f != variable for f in final]
        levels = [l for l, k in zip(levels, keep) if k]
        families = [f for f, k in zip(families, keep) if k]
        final = [f for f, k in zip(final, keep) if k]
        frame = frame.drop(columns=variable)
    levels = [0 if f in exogenous else l for f, l in zip(final, levels)]
    families = [None if l == 0 else f for f, l in zip(families, levels)]
    order = sorted(range(len(final)), key=lambda i: (levels[i], i), reverse=True)
    lines = []
    for i in order:
        explaining = [final[j] for j in order if levels[j] < levels[i]]
        if explaining:
            lines.append(f"{final[i]} ~ {'+'.join(explaining)} , {families[i]}")
    return frame, "\n".join(lines)


def _glm_graph(data, args, translated, initial, pmax, glm_coin):
    exogenous = _values(args, "exogenous")
    if not exogenous:
        raise ValueError("missing independent variables")
    chaine = list(args["chaine"])
    family = list(args["family"])
    frame, formulas = glm_equations(data, initial, exogenous, chaine, family, _values(args, "metric"),
                                    _values(args, "dichotomies"), _values(args, "valueDicho"))
    from .glmcoin import glm_coin as glm

    options = {k: v for k, v in translated.items() if k in _accepted(glm) and k not in ("data", "formulas")}
    if pmax is not None:
        options["pmax"] = pmax
    options["note"] = formulas.replace("\n", "</br>")
    weight = options.get("weights")
    if "weight" in translated and translated["weight"] is not None:
        options["weights"] = data[translated["weight"]].to_numpy(dtype=float)
    return glm(formulas, frame, **options)


def caring(filepath: str, arguments: Mapping, show: bool = False) -> str:
    """Read ``filepath`` and write the page of networks that ``arguments`` ask for (see :func:`caring_create_graphs`)."""
    return caring_create_graphs(caring_read_file(filepath), arguments, show=show)
