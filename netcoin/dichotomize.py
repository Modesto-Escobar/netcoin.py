"""Dichotomise categorical and textual variables (port of ``netCoin::dichotomize``)."""
from __future__ import annotations

import re
import unicodedata
import warnings
from collections import Counter
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

_TRIM = " \t\r\n"  # what R's trimws() removes


def _collation_key(term: str):
    """Approximate the locale collation R uses to sort terms.

    R orders the levels of ``table()`` with the collation of the session
    locale: case and accents are secondary (``"azul" < "Green" < "rojo"``) and a
    lowercase letter precedes its uppercase twin. Plain code-point order does
    not do that. This key reproduces that behaviour for the usual cases; it is
    not a full implementation of the platform collation.
    """
    decomposed = unicodedata.normalize("NFD", term)
    base = "".join(c for c in decomposed if not unicodedata.combining(c))
    return (base.casefold(), decomposed.casefold(), term.swapcase())


def _as_character(series: pd.Series) -> list[str]:
    """``as.character`` of a column; missing values become ``""`` as in R."""

    def convert(x):
        if x is None or (not isinstance(x, str) and pd.isna(x)):
            return ""
        if isinstance(x, (bool, np.bool_)):
            return "TRUE" if x else "FALSE"
        if isinstance(x, (float, np.floating)) and float(x).is_integer():
            return str(int(x))
        return str(x)

    return [convert(x) for x in series.tolist()]


def _recycle_sep(sep, variables: list[str]) -> dict[str, str]:
    if isinstance(sep, str):
        return {v: sep for v in variables}
    if isinstance(sep, dict):
        return {v: sep.get(v, "") for v in variables}
    if isinstance(sep, Sequence):
        sep = list(sep)
        if not sep:
            raise ValueError("sep must not be empty")
        if len(sep) == len(variables):
            return dict(zip(variables, sep))
        return {v: sep[0] for v in variables}  # R recycles sep[1]
    raise TypeError("sep must be a string, a sequence of strings or a dict")


def _clean_text(parts: list[str], text: str, stopwords: set[str]) -> list[str]:
    """R: lower/upper case, drop everything but letters, spaces and ``_``."""
    out = []
    for x in parts:
        x = x.lower() if text == "lower" else x.upper()
        x = "".join(ch for ch in x if ch.isalpha() or ch.isspace() or ch == "_")
        x = x.strip(_TRIM)
        if x != "" and x not in stopwords:
            out.append(x)
    return out


def _split_column(series: pd.Series, sep: str, regex: bool, text: str,
                  stopwords: set[str]) -> list[list[str]]:
    """The list of terms of every row of one variable."""
    rows = []
    for raw in _as_character(series):
        if raw.strip(_TRIM) == "":
            rows.append([])
            continue
        if sep not in ("", "C"):
            parts = re.split(sep, raw) if regex else raw.split(sep)
            if parts and parts[-1] == "":  # strsplit() drops one trailing empty piece
                parts.pop()
        else:
            parts = [raw]
        if text != "none":
            parts = _clean_text(parts, text, stopwords)
        else:
            # R keeps empty tokens ("a;;b"), which end up as an unnamed column;
            # they carry no information, so they are dropped here.
            parts = [p for p in parts if p != ""]
        rows.append(parts)
    return rows


def _unique(items: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(items))


def _concatenate_columns(df: pd.DataFrame, variables: list[str], nas: str) -> pd.Series:
    """``sep="C"``: join the non-missing values of each row with ``|``."""
    values = [
        [None if (not isinstance(x, str) and pd.isna(x)) else x for x in df[v].tolist()]
        for v in variables
    ]
    as_text = [_as_character(df[v]) for v in variables]
    joined = []
    for i in range(len(df)):
        present = [as_text[j][i] for j in range(len(variables)) if values[j][i] is not None]
        joined.append("|".join(present) if present else str(nas))
    return pd.Series(joined, index=df.index, name="String")


def dichotomize(
    data: pd.DataFrame,
    variables: str | Sequence[str],
    sep: str | Sequence[str] | dict = "",
    regex: bool = False,
    text: str = "none",
    min_freq: float = 1,
    length: int = 0,
    sort: bool = True,
    values: Sequence | None = None,
    stopwords: Iterable[str] | None = None,
    nas: str | None = "None",
    sparse: bool = False,
    add: bool = False,
) -> pd.DataFrame:
    """Convert categorical or textual columns into 0/1 indicator columns.

    Port of ``dichotomize`` from the R package netCoin. The columns are named
    after the categories (or the words) found in the selected variables.

    Parameters
    ----------
    data : DataFrame
        Table with the variables to convert. A cell can hold a single value or
        several events delimited by ``sep``.
    variables : str or list of str
        Columns to convert.
    sep : str, sequence of str or dict
        Separator that divides the events of a cell. ``""`` (default) takes the
        whole cell as one category. ``"C"`` concatenates the variables of each
        row before dichotomising (requires ``nas``). A sequence gives one
        separator per variable (a shorter one recycles its first element).
    regex : bool
        Treat ``sep`` as a regular expression. With several variables, their
        terms are pooled per row into a single group named after the first
        variable.
    text : {"none", "lower", "upper"}
        Case conversion. ``"lower"``/``"upper"`` also remove punctuation and
        numbers.
    min_freq : float
        Keep only terms with at least this frequency. A value in (0, 1) is a
        proportion of the rows. (``min`` in R.)
    length : int
        Maximum number of columns (most frequent terms first); 0 is no limit.
    sort : bool
        Order the columns by decreasing frequency. Otherwise they follow the
        alphabetical order of the terms. Ties are broken alphabetically.
    values : sequence, optional
        Restrict the conversion to these terms.
    stopwords : iterable of str, optional
        Words to remove; only applied when ``text`` is ``"lower"`` or ``"upper"``.
    nas : str or None
        Name of the column that flags the rows without any term; it is called
        ``"<variable>:<nas>"``. ``None`` omits it (always omitted if ``sparse``).
    sparse : bool
        Return pandas sparse columns instead of dense ones.
    add : bool
        Append the new columns to ``data`` (ignored if ``sparse``).

    Returns
    -------
    DataFrame
        0/1 integer columns, with the index of ``data``. With several
        independent variables the columns are prefixed with ``"<variable>:"``.

    Notes
    -----
    The order of terms with the same frequency follows an approximation of the
    locale collation used by R (case and accents are secondary).
    Empty tokens (``"a;;b"``) are discarded, whereas R keeps an unnamed column.
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError("You must pass a data frame!")
    if text not in ("none", "lower", "upper"):
        raise ValueError("text must be 'none', 'lower' or 'upper'")

    variables = [variables] if isinstance(variables, str) else list(variables)
    missing = [v for v in variables if v not in data.columns]
    if missing:
        raise KeyError(f"variables not in data: {missing}")

    original = data
    index = data.index
    n = len(data)
    if 0 < min_freq < 1:
        min_freq = min_freq * n

    sep_by_var = _recycle_sep(sep, variables)

    # sep = "C": concatenate the variables of each row (R: only if nas is given)
    first_sep = sep_by_var[variables[0]] if variables else ""
    if nas is not None and first_sep == "C":
        data = pd.DataFrame({"String": _concatenate_columns(data, variables, nas)})
        variables = ["String"]
        sep_by_var = {"String": "|"}

    stop = set()
    if stopwords is not None and text != "none":
        stop = {str(w).lower() if text == "lower" else str(w).upper() for w in stopwords}

    terms_by_var = {
        v: _split_column(data[v], sep_by_var[v], regex, text, stop) for v in variables
    }

    # terms of each row, per group of variables
    if regex and len(variables) > 1:
        pooled = [_unique(sum(rows, [])) for rows in zip(*(terms_by_var[v] for v in variables))]
        groups = {variables[0]: pooled}
    else:
        groups = {v: [_unique(r) for r in terms_by_var[v]] for v in variables}
    prefix = len(groups) > 1

    blocks: list[pd.DataFrame] = []
    for group, rows in groups.items():
        all_terms = [t for r in rows for t in r]
        if not all_terms:
            continue

        counts = Counter(all_terms)
        if values is None:
            valid = sorted(counts, key=_collation_key)  # order of R's table()
            if sort:
                valid.sort(key=lambda t: -counts[t])  # stable: ties stay alphabetical
            if min_freq > 0:
                valid = [t for t in valid if counts[t] >= min_freq]
            if length > 0 and len(valid) > length:
                valid = valid[:length]
        else:
            valid = [str(v) for v in values]
            if min_freq > 0 or length > 0:
                position = {t: i for i, t in enumerate(valid)}  # first match, as match()
                freq: Counter = Counter()
                for t in all_terms:
                    if t in position:
                        freq[position[t]] += 1
                idx = sorted(freq)
                if min_freq > 0:
                    idx = [i for i in idx if freq[i] >= min_freq]
                if sort:
                    idx.sort(key=lambda i: -freq[i])
                if length > 0 and len(idx) > length:
                    idx = idx[:length]
                valid = [valid[i] for i in idx]

        if not valid:
            continue

        column_of = {t: j for j, t in enumerate(valid)}
        row_idx, col_idx = [], []
        for i, r in enumerate(rows):
            for t in r:
                j = column_of.get(t)
                if j is not None:
                    row_idx.append(i)
                    col_idx.append(j)

        names = [f"{group}:{t}" for t in valid] if prefix else list(valid)

        if sparse:
            from scipy import sparse as sp

            matrix = sp.csr_matrix(
                (np.ones(len(row_idx), dtype=np.int8), (row_idx, col_idx)),
                shape=(n, len(valid)),
            )
            matrix.data[:] = 1
            blocks.append(pd.DataFrame.sparse.from_spmatrix(matrix, index=index, columns=names))
            continue

        dense = np.zeros((n, len(valid)), dtype=np.int64)
        dense[row_idx, col_idx] = 1
        block = pd.DataFrame(dense, index=index, columns=names)
        if nas is not None:
            block[f"{group}:{nas}"] = (dense.sum(axis=1) == 0).astype(np.int64)
        blocks.append(block)

    if not blocks:
        return original.copy()

    result = pd.concat(blocks, axis=1)

    if sparse:
        if add:
            warnings.warn("'add' is ignored when sparse=True", stacklevel=2)
        return result
    if add:
        return pd.concat([original, result], axis=1)
    return result
