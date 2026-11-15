"""Coincidence matrices: the ``Coin`` result type and ``coin`` (port of netCoin's ``coin``)."""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd

from .utils import frequencyNames, getByLanguage, nameList


def _format_count(x: float) -> str:
    return str(int(x)) if float(x).is_integer() else f"{x:g}"


@dataclass
class Coin:
    """Coincidences among a set of events or characteristics.

    This is what :func:`coin` returns; :func:`sim`, :func:`edge_list` and the
    rest of the functions take it as input. It plays the role of the ``coin``
    class of the R package.

    Attributes
    ----------
    matrix : DataFrame
        Square, symmetric table with the number of cases in which each pair of
        variables coincides. The diagonal holds the frequency of each variable.
    n : float
        Number of cases (the sum of the weights, if they were used).
    m : DataFrame, optional
        Pairwise: number of cases with no missing value in each pair.
    x : DataFrame, optional
        Pairwise: cases where the row variable is missing and the column
        variable is present.
    """

    matrix: pd.DataFrame
    n: float
    m: pd.DataFrame | None = None
    x: pd.DataFrame | None = None

    def __post_init__(self):
        mat = self.matrix
        if not isinstance(mat, pd.DataFrame) or mat.shape[0] != mat.shape[1]:
            raise ValueError("matrix must be a square DataFrame")
        if list(mat.index) != list(mat.columns):
            raise ValueError("matrix must have the same labels in rows and columns")
        if (self.m is None) != (self.x is None):
            raise ValueError("m and x must be given together (pairwise coincidences)")

    # ------------------------------------------------------------ building ----
    @classmethod
    def from_matrix(cls, matrix, n: float, labels: Sequence[str] | None = None) -> "Coin":
        """Build a :class:`Coin` from an existing coincidence matrix.

        ``matrix`` can be complete or only its lower triangle (blank or ``NaN``
        above the diagonal); ``n`` is the number of cases.
        """
        df = pd.DataFrame(matrix)
        if labels is None:
            labels = list(df.columns)
        df.columns = [str(c) for c in labels]
        df.index = df.columns
        values = df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        upper = np.triu_indices_from(values, k=1)
        fill = np.isnan(values[upper])
        values[upper] = np.where(fill, values.T[upper], values[upper])
        return cls(pd.DataFrame(values, index=df.index, columns=df.columns), float(n))

    # ---------------------------------------------------------- properties ----
    @property
    def labels(self) -> list[str]:
        return list(self.matrix.columns)

    @property
    def frequencies(self) -> pd.Series:
        """Number of cases of each variable (the diagonal)."""
        return pd.Series(np.diag(self.matrix.to_numpy()), index=self.matrix.index,
                         name="frequency")

    @property
    def pairwise(self) -> bool:
        return self.m is not None

    def __len__(self) -> int:
        return self.matrix.shape[0]

    # ------------------------------------------------------------- display ----
    def lower(self, decimals: int = 0) -> pd.DataFrame:
        """The lower triangle as text, blank above the diagonal (R's ``lower``)."""
        a = self.matrix.to_numpy(dtype=float)
        text = np.empty(a.shape, dtype=object)
        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                text[i, j] = "" if j > i else f"{a[i, j]:.{decimals}f}"
        return pd.DataFrame(text, index=self.matrix.index, columns=self.matrix.columns)

    def summary(self) -> pd.Series:
        """Relative frequency of each variable (``diag / n``), as R's ``summary``."""
        return self.frequencies / self.n

    def __repr__(self) -> str:
        head = f"n= {_format_count(self.n)}"
        if self.pairwise:
            head += " (pairwise)"
        return head + "\n" + self.lower(0).to_string()

    def to_frame(self) -> pd.DataFrame:
        """A copy of the full symmetric matrix."""
        return self.matrix.copy()


def _to_float_array(df: pd.DataFrame) -> np.ndarray:
    """Incidences as floats; missing values are ``NaN``."""
    columns = []
    for col in df.columns:
        s = df[col]
        if isinstance(s.dtype, pd.SparseDtype):
            s = s.sparse.to_dense()
        columns.append(pd.to_numeric(s, errors="coerce").to_numpy(dtype=float, na_value=np.nan))
    return np.column_stack(columns) if columns else np.empty((len(df), 0))


def coin(
    incidences: pd.DataFrame | np.ndarray,
    minimum: float = 1,
    maximum: float | None = None,
    sort: bool = False,
    decreasing: bool = True,
    total: bool = False,
    subsample: bool = False,
    weight: Sequence[float] | None = None,
    pairwise: bool = False,
) -> Coin:
    """Count the coincidences among the variables of an incidence table.

    Port of ``coin`` from the R package netCoin.

    Parameters
    ----------
    incidences : DataFrame or ndarray
        One row per case, one 0/1 column per variable (see :func:`dichotomize`).
        Missing values are allowed.
    minimum, maximum : float
        Keep the variables whose frequency lies in ``[minimum, maximum]``.
        ``maximum`` defaults to the number of cases (the total weight with
        ``weight``).
    sort, decreasing : bool
        Order the variables by frequency (stable for ties).
    total : bool
        Add a first variable ``"Total"`` that is present in every case.
        Not applicable with weights.
    subsample : bool
        Use only the cases with at least one event.
    weight : sequence of float, optional
        Weight of each case.
    pairwise : bool
        Treat missing values pairwise (each pair uses the cases where both
        variables are known) instead of dropping the incomplete cases.

    Returns
    -------
    Coin
    """
    df = incidences if isinstance(incidences, pd.DataFrame) else pd.DataFrame(incidences)
    labels = [str(c) for c in df.columns]
    x = _to_float_array(df)

    known = ~np.isnan(x)
    if not np.all(np.isnan(x) | (x == 0) | (x == 1)):
        raise ValueError("all data in the incidence matrix has to be dichotomous (0/1); "
                         "see dichotomize()")

    w = None
    if weight is not None:
        w = np.asarray(weight, dtype=float)
        if w.shape != (x.shape[0],):
            raise ValueError("weight must have one value per row of incidences")

    if pairwise:
        nomiss = known.astype(float)
        x = np.where(known, x, 0.0)
    else:
        complete = known.all(axis=1)
        x = x[complete]
        w = w[complete] if w is not None else None

    if subsample:
        used = x.sum(axis=1) > 0
        x = x[used]
        w = w[used] if w is not None else None
        if pairwise:
            nomiss = nomiss[used]

    if total:
        if w is not None:
            warnings.warn("total cannot be applied in weighted tables", stacklevel=2)
            total = False
        else:
            x = np.column_stack([np.ones(x.shape[0]), x])
            labels = ["Total"] + labels
            if pairwise:
                nomiss = np.column_stack([np.ones(nomiss.shape[0]), nomiss])

    rows = x.shape[0]
    if w is None:
        f = x.T @ x
        n = float((nomiss.sum(axis=1) > 0).sum()) if pairwise else float(rows)
    else:
        f = x.T @ (w[:, None] * x)
        n = float(w.sum())

    if maximum is None:
        maximum = n if w is not None else rows
    d = np.diag(f)
    order = np.arange(len(d))
    if sort:
        order = np.argsort(-d if decreasing else d, kind="stable")
    keep = [i for i in order if minimum <= d[i] <= maximum or (total and labels[i] == "Total")]
    if total and "Total" in labels:
        keep = [labels.index("Total")] + [i for i in keep if labels[i] != "Total"]
    if not keep:
        raise ValueError("no variables left: check minimum and maximum")

    def frame(values: np.ndarray) -> pd.DataFrame:
        if w is None:
            values = np.rint(values).astype(np.int64)
        names = [labels[i] for i in keep]
        return pd.DataFrame(values[np.ix_(keep, keep)], index=names, columns=names)

    m = xm = None
    if pairwise:
        if w is None:
            m = nomiss.T @ nomiss
            xm = (1 - nomiss).T @ x
        else:
            m = nomiss.T @ (w[:, None] * nomiss)
            xm = (1 - nomiss).T @ (w[:, None] * x)
        m, xm = frame(m), frame(xm)

    return Coin(frame(f), n, m, xm)


# --------------------------------------------------------------- coocur ----
@dataclass
class Cooc:
    """Result of :func:`coocur`: occurrences shared by pairs of variables.

    ``n`` is the total number of occurrences and ``m`` the sum, over the cases,
    of the largest occurrence of each case.
    """

    matrix: pd.DataFrame
    n: float
    m: float

    def __repr__(self) -> str:
        a = self.matrix.to_numpy(dtype=float)
        text = np.where(np.tril(np.ones(a.shape, dtype=bool)), np.vectorize(_format_count)(a), "")
        table = pd.DataFrame(text, index=self.matrix.index, columns=self.matrix.columns)
        return f"n= {_format_count(self.n)}; m= {_format_count(self.m)}\n" + table.to_string()


def coocur(
    ocurrences: pd.DataFrame | np.ndarray,
    minimum: float = 1,
    maximum: float = np.inf,
    sort: bool = False,
    decreasing: bool = True,
) -> Cooc:
    """Co-occurrence of counts: for each pair, the sum over cases of the smaller count.

    Port of ``coocur`` from the R package netCoin. Unlike :func:`coin`, the
    variables need not be 0/1.
    """
    df = ocurrences if isinstance(ocurrences, pd.DataFrame) else pd.DataFrame(ocurrences)
    x = df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    labels = [str(c) for c in df.columns]
    k = x.shape[1]
    result = np.empty((k, k))
    for i in range(k):
        for j in range(i, k):
            result[i, j] = result[j, i] = np.sum(np.minimum(x[:, i], x[:, j]))

    d = np.diag(result)
    order = np.arange(k)
    if sort:
        order = np.argsort(-d if decreasing else d, kind="stable")
    keep = [i for i in order if minimum <= d[i] <= maximum]
    if not keep:
        raise ValueError("no variables left: check minimum and maximum")
    names = [labels[i] for i in keep]
    sub = x[:, keep]
    return Cooc(pd.DataFrame(result[np.ix_(keep, keep)], index=names, columns=names),
                float(np.sum(sub)), float(np.sum(np.max(sub, axis=1))))


# ------------------------------------------------------------- prop_coin ----
def prop_coin(x: Coin, margin: int = 0, decimals: int = 1) -> pd.DataFrame:
    """Coincidences as percentages.

    Port of ``propCoin`` from the R package netCoin.

    Parameters
    ----------
    x : Coin
    margin : {0, 1, 2}
        0: percentage of the cases (of the cases with data in both variables, if
        pairwise); 1: of the frequency of the row variable; 2: of the column's.
    decimals : int
        Rounding.
    """
    if not isinstance(x, Coin):
        raise TypeError("x must be a Coin object (see coin())")
    a = x.matrix.to_numpy(dtype=float)
    diag = np.diag(a)
    if margin == 0:
        n = x.m.to_numpy(dtype=float) if x.pairwise else x.n
        values = 100 * a / n
    elif margin == 1:
        values = 100 * a / diag[:, None]
    elif margin == 2:
        values = 100 * a / diag[None, :]
    else:
        raise ValueError("margin must be 0, 1 or 2")
    return pd.DataFrame(np.round(values, decimals), index=x.matrix.index, columns=x.matrix.columns)


# -------------------------------------------------------------- as_nodes ----
def as_nodes(
    data: Coin | pd.DataFrame,
    frequency: bool = True,
    percentages: bool = False,
    language: str = "en",
) -> pd.DataFrame:
    """Table of nodes (name, frequency and/or percentage) from a :class:`Coin`.

    Port of ``asNodes`` from the R package netCoin. If ``data`` is an edge list
    (with ``Source`` and ``Target`` columns) the nodes are the sorted names found
    in it.
    """
    name = getByLanguage(nameList, language)
    if isinstance(data, Coin):
        out = pd.DataFrame({name: data.labels})
        divider = np.diag(data.m.to_numpy(dtype=float)) if data.pairwise else data.n
        if frequency:
            out[getByLanguage(frequencyNames, language)] = data.frequencies.to_numpy()
        if percentages:
            out["%"] = data.frequencies.to_numpy(dtype=float) / divider * 100
        return out
    if isinstance(data, pd.DataFrame) and {"Source", "Target"} <= set(data.columns):
        names = sorted(set(data["Source"]) | set(data["Target"]))
        return pd.DataFrame({name: names})
    raise TypeError("data is neither a Coin object nor an edge DataFrame (Source, Target)")
