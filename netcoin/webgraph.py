"""Shared pieces of the graphs that are not networks: bar, timeline and pie charts.

These helpers port the checks and the JSON conversion of the R package rD3plot,
whose JavaScript reads exactly that format.
"""
from __future__ import annotations

import json
import math
import os
import re
import warnings
import webbrowser
from tempfile import mkdtemp

import numpy as np
import pandas as pd

from .utils import createHTML, getLanguageScript

# ------------------------------------------------------------------ colors ----
_COLOR_TABLE: dict[str, str] | None = None

# R's category palette (rD3plot's categoryColors)
CATEGORY_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#e377c2", "#d62728", "#bcbd22", "#9467bd", "#8c564b",
    "#7f7f7f", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#f7b6d2", "#ff9896", "#dbdb8d",
    "#c5b0d5", "#c49c94", "#c7c7c7", "#9edae5",
]


def _color_table() -> dict[str, str]:
    global _COLOR_TABLE
    if _COLOR_TABLE is None:
        path = os.path.join(os.path.dirname(__file__), "colors.json")
        with open(path, encoding="utf-8") as f:
            _COLOR_TABLE = json.load(f)
    return _COLOR_TABLE


def col2hex(color: str) -> str:
    """The ``#rrggbb`` of a colour given as in R: ``"#1F77B4"``, ``"#abc"`` or a name.

    Names are R's (``colors()``), case and spaces are ignored; the alpha channel
    of ``#rrggbbaa`` is dropped. Raises ``ValueError`` for anything else.
    """
    text = str(color).strip().lower().replace(" ", "")
    if text.startswith("#"):
        digits = text[1:]
        if re.fullmatch(r"[0-9a-f]{3,4}", digits):
            digits = "".join(c * 2 for c in digits[:3])
        elif re.fullmatch(r"[0-9a-f]{6}([0-9a-f]{2})?", digits):
            digits = digits[:6]
        else:
            raise ValueError(f"invalid color: {color!r}")
        return "#" + digits
    if text == "transparent":
        return "#ffffff"
    try:
        return _color_table()[text]
    except KeyError:
        raise ValueError(f"invalid color: {color!r}") from None


def is_color(color) -> bool:
    """True if every element of ``color`` is a valid colour."""
    items = [color] if isinstance(color, str) else list(color) if color is not None else []
    if not items:
        return False
    try:
        for c in items:
            if not isinstance(c, str):
                return False
            col2hex(c)
    except ValueError:
        return False
    return True


def category_colors(items) -> list[str]:
    """One colour of the category palette per item (numbers are used as indices)."""
    values = pd.Series(items)
    if not pd.api.types.is_numeric_dtype(values):
        values = pd.Series(pd.Categorical(values).codes + 1, index=values.index)
    index = ((values.astype(int) - 1) % len(CATEGORY_COLORS)).to_numpy()
    return [CATEGORY_COLORS[i] for i in index]


# ------------------------------------------------------------ small checks ----
def check_cex(cex) -> float:
    if isinstance(cex, bool) or not isinstance(cex, (int, float)):
        warnings.warn("cex: must be numeric", stacklevel=3)
        return 1
    if cex > 2:
        warnings.warn("cex: the font may be too large for proper display", stacklevel=3)
    return cex


def check_language(language) -> str:
    language = language[0] if isinstance(language, (list, tuple)) else language
    if language not in ("en", "es", "ca"):
        warnings.warn(f"language: '{language}' is not supported", stacklevel=3)
        return "en"
    return language


def check_default_color(color) -> str | None:
    if color is None:
        return None
    if is_color(color):
        return col2hex(color)
    warnings.warn("defaultColor: you must pass a valid color", stacklevel=3)
    return None


def column_option(options: dict, item: str, variable) -> dict:
    """Set ``options[item]`` to a column name, or remove it if it is empty/``False``."""
    if variable is None or variable is False or (isinstance(variable, str) and variable == ""):
        options.pop(item, None)
    else:
        options[item] = variable
    return options


def check_item_value(items: pd.DataFrame, options: dict, prop: str, value, name: str):
    """Colour (or similar) of every row of ``items``: port of rD3plot's ``checkItemValue``.

    ``value`` can be ``None`` (remove), the name of a column of ``items``, one colour,
    a vector with one value per row (colours, or numbers/categories that get a
    colour of the category palette) or a DataFrame with several of those vectors.
    Returns the (possibly extended) ``items`` and ``options``.
    """
    items = items.copy()
    options = dict(options)

    def prepare(values):
        series = pd.Series(list(values))
        if pd.api.types.is_numeric_dtype(series) or isinstance(values, pd.Categorical) or \
                isinstance(getattr(values, "dtype", None), pd.CategoricalDtype):
            return category_colors(values)
        if all(isinstance(v, str) for v in series) and is_color(list(series)):
            return [col2hex(v) for v in series]
        warnings.warn(f"{name}: this value cannot be a {name}", stacklevel=4)
        return None

    if value is None:
        options.pop(prop, None)
    elif isinstance(value, pd.DataFrame):
        if len(value) == len(items):
            for key in value.columns:
                colors = prepare(value[key])
                if colors is not None:
                    legend = items[key].astype(str) if key in items.columns else value[key].astype(str)
                    items[key] = list(legend)
                    items[f"_{prop}_{key}"] = colors
            options[prop] = str(value.columns[0])
        else:
            warnings.warn(f"{name}: number of rows doesn't match with the items", stacklevel=3)
    elif (not isinstance(value, str) and np.ndim(value) > 0 and len(value) > 1) or \
            (isinstance(value, str) and is_color(value) and value not in items.columns):
        values = [value] * len(items) if isinstance(value, str) else list(value)
        if len(values) == len(items):
            legend = [str(v) for v in values]
            colors = prepare(values if not isinstance(value, pd.Series) else value)
            if colors is not None:
                items[f"-{prop}-"] = legend
                items[f"_{prop}_-{prop}-"] = colors
                options[prop] = f"-{prop}-"
        else:
            warnings.warn(f"{name}: value length doesn't match with the number of rows", stacklevel=3)
    else:
        options = column_option(options, prop, value)
    return items, options


# -------------------------------------------------------------------- JSON ----
_NUMBER = re.compile(r"^\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?\s*$")
_HEX = re.compile(r"^\s*[-+]?0[xX][0-9a-fA-F]+\s*$")
_INFINITE = re.compile(r"^\s*[-+]?(inf|infinity)\s*$", re.IGNORECASE)


def _number(x: float):
    """R's ``prepare_number``: integers as they are, others with 4 significant digits."""
    if math.isnan(x) or math.isinf(x):
        return None
    if x == int(x):
        return int(x)
    return float(f"{x:.4g}")


def r_json(obj):
    """Convert ``obj`` as rD3plot's ``toJSON`` would, ready for ``json.dumps``.

    Non-integer numbers keep four significant digits, missing values become
    ``null`` and a text that is a number (``"1990"``) is written as a number,
    which is what the JavaScript of rD3plot receives from R.
    """
    if obj is None:
        return None
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return _number(float(obj))
    if isinstance(obj, str):
        if obj is pd.NA:
            return None
        if _INFINITE.match(obj):
            return None
        if _NUMBER.match(obj):
            return _number(float(obj))
        if _HEX.match(obj):
            return _number(float(int(obj.strip(), 16)))
        return obj
    if obj is pd.NA or obj is pd.NaT:
        return None
    if isinstance(obj, dict):
        return {str(k): r_json(v) for k, v in obj.items()}
    if isinstance(obj, np.ndarray):
        return r_json(obj.tolist())
    if isinstance(obj, (list, tuple)):
        return [r_json(v) for v in obj]
    if isinstance(obj, pd.Series):
        return r_json(obj.tolist())
    if isinstance(obj, pd.DataFrame):
        return r_json(records(obj))
    return r_json(str(obj))


def records(frame: pd.DataFrame) -> list[dict]:
    """A table as a list of rows (``[{"col": value, ...}, ...]``), as R writes data frames.

    Values are kept as they are (missing ones are ``None``); :func:`r_json` rounds them
    when the page is written.
    """
    names = [str(c) for c in frame.columns]
    rows = frame.astype(object).where(frame.notna(), None).to_numpy().tolist()
    return [dict(zip(names, row)) for row in rows]


def columns(frame: pd.DataFrame) -> list[list]:
    """A table as a list of columns (what the network and pie scripts read)."""
    return [frame[c].astype(object).where(frame[c].notna(), None).tolist() for c in frame.columns]


# ------------------------------------------------------------------- pages ----
class WebGraph:
    """Base of the graphs that are written as a page: bar, timeline and pie charts."""

    options: dict
    _engine = ""            # the JavaScript that draws it
    _needs_jspdf = True

    def _data(self) -> dict:
        raise NotImplementedError

    def _create(self, directory, title=None):
        """Write the web page of this graph (also used by multigraphCreate)."""
        language = getLanguageScript(self) or "en.js"
        scripts = ["d3.min.js"] + (["jspdf.min.js"] if self._needs_jspdf else []) + [
            "images.js", "colorScales.js", "functions.js", language, self._engine]
        createHTML(directory, ["reset.css", "styles.css"], scripts,
                   lambda: r_json(self._data()), title)

    def plot(self, directory: str | None = None):
        """Write the page (in a temporary folder by default) and open it in the browser."""
        directory = directory or mkdtemp()
        self._create(directory)
        webbrowser.open(os.path.join(directory, "index.html"))

    def to_json(self) -> str:
        """The JSON that the page receives."""
        return json.dumps(r_json(self._data()), ensure_ascii=False)
