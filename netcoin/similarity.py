"""Similarities and edge lists from coincidences (port of netCoin's ``sim``/``edgeList``)."""
from __future__ import annotations

import math
import warnings
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from scipy.stats import t as student_t

from .coincidence import Coin, coin

# ------------------------------------------------------------- procedures ----
# R code -> name of the matrix that ``sim`` returns for it
_CODE_NAMES = {
    "M": "matching", "T": "Rogers", "G": "Gower", "S": "Sneath", "B": "Anderberg",
    "J": "Jaccard", "D": "dice", "A": "antiDice", "O": "Ochiai", "K": "Kulczynski",
    "N": "Hamann", "Y": "Yule", "P": "Pearson", "V": "tetrachoric", "C": "odds",
    "R": "Russell", "E": "expected", "L": "confidence", "H": "Haberman", "Z": "Z",
    "W": "Fisher", "F": "coincidences", "X": "relative", "I": "sConditional",
    "0": "tConditional", "U": "c.conditional", "Q": "c.probable",
}
_ALL_CODES = ["M", "T", "G", "S", "B", "J", "D", "A", "O", "K", "N", "Y", "P", "V", "C",
              "R", "E", "L", "H", "Z", "F", "X", "I", "0", "Q", "U"]
# exact names (case-insensitive) are accepted before R's abbreviation rule, which
# would read e.g. "c.conditional" as "C" (odds) and "sConditional" as "S" (Sneath)
_EXACT = {name.lower(): code for code, name in _CODE_NAMES.items()}
_EXACT.update({"conf.l": "L", "conf.u": "L", "p(z)": "Z", "p(fisher)": "W"})
# R's c_method(): substitutions applied before taking the first letter
_SUBSTITUTIONS = [("FIS", "W"), ("HYP", "W"), ("HAM", "NNH"), ("CC", "UC"), ("CP", "QC"),
                  ("OD", "CD"), ("RO", "TA"), ("AND", "BER"), ("TET", "VET"),
                  ("CON", "LCO"), ("II", "0"), ("COI", "F")]
_PVALUE_CODES = {"Z", "W"}

_CONDITIONAL_LABELS = ["Null", "Mere", "Conditional", "Significant", "Quite significant",
                       "Very significant", "Subtotal", "Suptotal", "Total"]
_PROBABLE_LABELS = ["Null", "Mere", "Probable", "Significant", "Quite significant",
                    "Very significant", "Subtotal", "Suptotal", "Total"]
# names of the columns of an edge list
_EDGE_RENAME = {"Z": "p(Z)", "Fisher": "p(Fisher)"}


def _method_code(procedure: str) -> str:
    """The one-letter code of a procedure, as netCoin's ``c_method`` does.

    Accepts the exact name of the matrix (``"Jaccard"``, ``"c.conditional"``), or
    R's abbreviation rule: the first letter of the name in capitals after a few
    substitutions (``"Ham"`` is Hamann, ``"CC"`` is c.conditional, ``"Hyp"`` is
    Fisher, ``"Con"`` is confidence...).
    """
    key = str(procedure).strip()
    if key.lower() in _EXACT:
        return _EXACT[key.lower()]
    s = key.upper()
    for old, new in _SUBSTITUTIONS:
        s = s.replace(old, new, 1)
    code = s[:1]
    if code not in _CODE_NAMES:
        raise ValueError(
            f"Unknown procedure '{procedure}'. Valid names: {sorted(_CODE_NAMES.values())}"
        )
    return code


def _method_codes(procedures) -> list[str]:
    if isinstance(procedures, str):
        procedures = [procedures]
    procedures = list(procedures)
    if procedures and str(procedures[0]).lower() == "all":
        return list(_ALL_CODES)
    return list(dict.fromkeys(_method_code(p) for p in procedures))


def edge_column_name(procedure: str) -> str:
    """Name of the column of an edge list that holds ``procedure`` (e.g. ``"Z"`` is ``"p(Z)"``)."""
    name = _CODE_NAMES[_method_code(procedure)]
    return _EDGE_RENAME.get(name, name)


def _result_names(code: str) -> list[str]:
    return ["conf.L", "confidence", "conf.U"] if code == "L" else [_CODE_NAMES[code]]


def _check_level(level: float) -> float:
    if 1 <= level < 100:
        level = level / 100
    if level <= 0 or level >= 1:
        warnings.warn("Not valid level; 0.95 is used", stacklevel=3)
        return 0.95
    return float(level)


# ----------------------------------------------------------------- counts ----
def _safe_div(num, den):
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full(np.broadcast(num, den).shape, np.nan)
    np.divide(num, den, out=out, where=den != 0)
    return out


def _prepare_counts(data, minimum, maximum, sort, decreasing, weight, pairwise):
    """The 2x2 tables of every pair: a (both), b, c (only one), d (neither), N."""
    if isinstance(data, Coin):
        base = data
        if weight is not None:
            warnings.warn("weight is ignored: a Coin object already carries its n", stacklevel=3)
        if pairwise and base.m is None:
            raise ValueError("pairwise=True needs a Coin built with coin(..., pairwise=True)")
    elif isinstance(data, (pd.DataFrame, np.ndarray)):
        base = coin(data, minimum, maximum, sort, decreasing, weight=weight, pairwise=pairwise)
    else:
        raise TypeError("data must be a Coin object or an incidence DataFrame (see coin())")

    a_df = base.matrix
    diag = np.diag(a_df.to_numpy(dtype=float))
    idx = np.flatnonzero((diag >= minimum) & (diag <= maximum))
    if len(idx) == 0:
        raise ValueError("no variables left: check minimum and maximum")
    if sort or decreasing:
        d_sel = diag[idx]
        idx = idx[np.argsort(-d_sel if decreasing else d_sel, kind="stable")]

    a = a_df.to_numpy(dtype=float)[np.ix_(idx, idx)]
    labels = list(a_df.index[idx])
    if pairwise:
        N = base.m.to_numpy(dtype=float)[np.ix_(idx, idx)]
        X = base.x.to_numpy(dtype=float)[np.ix_(idx, idx)]
        nn = float(np.max(a_df.to_numpy(dtype=float)))
    else:
        N = np.full(a.shape, float(base.n))
        X = np.zeros(a.shape)
        nn = float(base.n)

    b = np.diag(a)[:, None] - a - X.T
    c = b.T
    d = N - a - b - c
    return {"a": a, "b": b, "c": c, "d": d, "N": N, "n": float(base.n), "nn": nn,
            "labels": labels, "pairwise": pairwise}


# -------------------------------------------------------------------- sim ----
def sim(
    data: Coin | pd.DataFrame | np.ndarray,
    procedures: str | Sequence[str] = "Jaccard",
    level: float = 0.95,
    distance: bool = False,
    minimum: float = 1,
    maximum: float = np.inf,
    sort: bool = False,
    decreasing: bool = False,
    weight: Sequence[float] | None = None,
    pairwise: bool = False,
) -> pd.DataFrame | dict[str, pd.DataFrame]:
    """Similarities, associations and statistics for every pair of variables.

    Port of ``sim`` from the R package netCoin.

    Parameters
    ----------
    data : Coin or DataFrame
        A :class:`Coin` object or a table of incidences (see :func:`coin`).
    procedures : str or list of str
        What to compute: ``"matching"``, ``"Rogers"``, ``"Gower"``, ``"Sneath"``,
        ``"Anderberg"``, ``"Jaccard"``, ``"dice"``, ``"antiDice"``, ``"Ochiai"``,
        ``"Kulczynski"``, ``"Hamann"``, ``"Yule"``, ``"Pearson"``,
        ``"tetrachoric"``, ``"odds"``, ``"Russell"``, ``"expected"``,
        ``"confidence"`` (also ``conf.L`` and ``conf.U``), ``"Haberman"``, ``"Z"``,
        ``"Fisher"``, ``"coincidences"``, ``"relative"``, ``"sConditional"``,
        ``"tConditional"``, ``"c.conditional"`` and ``"c.probable"``; or ``"all"``.
        R's short forms work too (``"J"``, ``"Ham"``, ``"CC"``, ``"Hyp"``...).
    level : float
        Confidence level for ``"confidence"``.
    distance : bool
        Return distances (``1 - similarity``) for the similarity coefficients.
    minimum, maximum : float
        Keep the variables whose frequency lies in that range.
    sort, decreasing : bool
        Order the variables by frequency.
    weight, pairwise
        As in :func:`coin` (used when ``data`` is a table of incidences;
        ``pairwise=True`` with a :class:`Coin` needs one built pairwise).

    Returns
    -------
    DataFrame or dict of DataFrame
        One square matrix per procedure: a DataFrame if only one was asked for,
        otherwise a dict keyed by the name of the matrix.
    """
    level = _check_level(level)
    codes = _method_codes(procedures)
    k = _prepare_counts(data, minimum, maximum, sort, decreasing, weight, pairwise)
    a, b, c, d, N = k["a"], k["b"], k["c"], k["d"], k["N"]
    labels = k["labels"]

    def frame(values):
        return pd.DataFrame(values, index=labels, columns=labels)

    def dist(s):
        return 1.0 - s if distance else s

    out: dict[str, pd.DataFrame] = {}
    with np.errstate(all="ignore"):
        m = np.where(a + d == N, 1.0, np.where(b + c == N, -1.0, 0.0))
        ab, ac = a + b, a + c
        den = np.sqrt(ab * ac * (b + d) * (d + c))      # denominator of phi
        phi = _safe_div(a * d - b * c, den)
        z_stat = np.sqrt(N) * phi
        df_n = np.where(N > 0, N, np.nan)

        def fill_special(s, special):
            return np.where(np.isnan(s), special, s)

        for code in codes:
            if code == "M":
                out["matching"] = frame(dist(_safe_div(a + d, a + b + c + d)))
            elif code == "T":
                out["Rogers"] = frame(dist(_safe_div(a + d, a + 2 * (b + c) + d)))
            elif code == "G":
                s = _safe_div(a * d, np.sqrt(ab * ac * (d + b) * (d + c)))
                out["Gower"] = frame(fill_special(dist(s), dist(np.maximum(m, 0))))
            elif code == "S":
                out["Sneath"] = frame(dist(_safe_div(2 * (a + d), 2 * (a + d) + b + c)))
            elif code == "B":
                s = (_safe_div(a, ab) + _safe_div(a, ac) + _safe_div(d, c + d)
                     + _safe_div(d, b + d)) / 4
                out["Anderberg"] = frame(fill_special(dist(s), dist(np.maximum(m, 0))))
            elif code == "J":
                out["Jaccard"] = frame(dist(_safe_div(a, a + b + c)))
            elif code == "D":
                out["dice"] = frame(dist(_safe_div(2 * a, 2 * a + b + c)))
            elif code == "A":
                out["antiDice"] = frame(dist(_safe_div(a, a + 2 * (b + c))))
            elif code == "O":
                s = _safe_div(a, np.sqrt(ab * ac))
                out["Ochiai"] = frame(fill_special(dist(s), dist(np.maximum(m, 0))))
            elif code == "K":
                s = (_safe_div(a, ab) + _safe_div(a, ac)) / 2
                out["Kulczynski"] = frame(fill_special(dist(s), dist(np.maximum(m, 0))))
            elif code == "N":
                out["Hamann"] = frame(dist(_safe_div(a - (b + c) + d, a + b + c + d)))
            elif code == "Y":
                s = _safe_div(a * d - b * c, a * d + b * c)  # R does not apply distance here
                out["Yule"] = frame(fill_special(s, m))
            elif code == "P":
                out["Pearson"] = frame(fill_special(dist(phi), dist(m)))
            elif code == "V":
                ratio = (a * d / (b * c)) ** (math.pi / 4)
                s = (ratio - 1) / (ratio + 1)
                out["tetrachoric"] = frame(fill_special(s, dist(m)))
            elif code == "C":
                s = (np.maximum(a, 0.5) * np.maximum(d, 0.5)) / (
                    np.maximum(b, 0.5) * np.maximum(c, 0.5))
                if distance:
                    s = -s
                np.fill_diagonal(s, -np.inf if distance else np.inf)
                out["odds"] = frame(s)
            elif code == "R":
                s = dist(_safe_div(a, a + b + c + d))
                if not distance:
                    np.fill_diagonal(s, 1.0)
                out["Russell"] = frame(s)
            elif code == "E":
                out["expected"] = frame(ab * ac / N)
            elif code == "L":
                df_t = np.where(N - 1 > 0, N - 1, np.nan)
                exp_ = ab * ac / N
                se = np.sqrt(exp_ * ((1 - ab / N) * (1 - ac / N)))
                t_two = student_t.ppf(level + (1 - level) / 2, df_t)
                sign = 2.0 * (exp_ < a) - 1
                conf = np.maximum(exp_ + sign * student_t.ppf(level, df_t) * se, 0)
                np.fill_diagonal(conf, np.diag(a))
                out["conf.L"] = frame(np.maximum(a - t_two * se, 0))
                out["confidence"] = frame(conf)
                out["conf.U"] = frame(np.minimum(a + t_two * se, N))
            elif code == "H":
                s = z_stat.copy()
                s = np.where(np.isnan(s), math.sqrt(k["nn"] if k["pairwise"] else k["n"]), s)
                if distance:
                    s = (N + s) / (2 * N)
                out["Haberman"] = frame(s)
            elif code == "Z":
                s = 1 - student_t.cdf(z_stat, df_n)
                out["Z"] = frame(np.where(np.isnan(s), 0.0, s))
            elif code == "W":
                out["Fisher"] = frame(_fisher(a, b, c, N))
            elif code == "F":
                out["coincidences"] = frame(a.copy())
            elif code == "X":
                out["relative"] = frame(a / N * 100)
            elif code == "I":
                out["sConditional"] = frame(_safe_div(a, ac) * 100)
            elif code == "0":
                out["tConditional"] = frame(_safe_div(a, ab) * 100)
            elif code == "U":
                z = 1 - student_t.cdf(z_stat, df_n)
                val = _categories(a, b, c, z, with_b=True)
                out["c.conditional"] = frame(val)
            elif code == "Q":
                t_val = _safe_div(_safe_div(a, ac) - 0.5, 1.0 / (2.0 * np.sqrt(ac)))
                z = 1 - student_t.cdf(t_val, np.where(ac > 0, ac, np.nan))
                val = _categories(a, b, c, z, with_b=False)
                out["c.probable"] = frame(val)

    return next(iter(out.values())) if len(out) == 1 else out


def _categories(a, b, c, z, with_b):
    """Strength of the conditional relationship, 0 (null) to 8 (total)."""
    conditions = [b + c == 0, c == 0]
    values = [8, 7]
    if with_b:
        conditions.append(b == 0)
        values.append(6)
    early = np.logical_or.reduce(conditions)
    conditions += [z < 0.001, z < 0.01, z < 0.05, z < 0.5, a > 0]
    values += [5, 4, 3, 2, 1]
    out = np.select(conditions, values, default=0).astype(float)
    return np.where(np.isnan(z) & ~early, np.nan, out)  # R: NA when Z is undefined


def _fisher(a, b, c, N):
    """One-sided Fisher exact p-value, ``1 - phyper(a-1, ...)`` of R.

    With weights the counts need not be integers; as R's ``phyper`` does, the
    quantile is floored and the numbers of white and black balls are rounded.
    """
    n_min = np.minimum(a + b, a + c)
    n_max = np.maximum(a + b, a + c)
    s = np.full(a.shape, np.nan)
    ok = np.isfinite(N) & (N > 0) & np.isfinite(n_min) & np.isfinite(a)
    white = np.rint(n_min[ok])
    black = np.rint(N[ok] - n_min[ok])
    s[ok] = hypergeom.sf(np.floor(a[ok] - 1 + 1e-7), white + black, white, np.rint(n_max[ok]))
    return s


# ------------------------------------------------------------- edge lists ----
def _edge_mask(nrow: int, ncol: int, row_labels, col_labels, directed, diagonal):
    """Which cells of the matrix become edges (R's ``mats2edges``)."""
    if nrow < ncol:  # asymmetric matrix: every column that is not a row, plus the lower part
        mask = np.ones((nrow, ncol), dtype=bool)
        position = {name: j for j, name in enumerate(col_labels)}
        cols = [position[name] for name in row_labels if name in position]
        if len(cols) == nrow:
            mask[:, cols] = np.tril(np.ones((nrow, nrow), dtype=bool), k=0 if diagonal else -1)
        return mask
    lower = np.tril(np.ones((nrow, ncol), dtype=bool), k=0 if diagonal else -1)
    if directed:
        return lower | np.triu(np.ones((nrow, ncol), dtype=bool), k=1)
    return lower


def _empty_edges(columns=("Source", "Target")) -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype=object) for c in columns})


def _matrices_to_edges(matrices: dict[str, pd.DataFrame], counts: np.ndarray,
                       criteria: str, min_val, max_val, support, directed, diagonal):
    first = next(iter(matrices.values()))
    row_labels, col_labels = list(first.index), list(first.columns)
    mask = _edge_mask(len(row_labels), len(col_labels), row_labels, col_labels,
                      directed, diagonal)

    crit = matrices[criteria].to_numpy(dtype=float)
    keep = mask & (crit <= max_val) & (crit >= min_val) & (counts >= support)
    keep &= ~np.isnan(crit)
    flat = np.flatnonzero(keep.ravel(order="F"))           # column by column, as R
    rows, cols = np.unravel_index(flat, keep.shape, order="F")

    edges = {"Source": [col_labels[j] for j in cols], "Target": [row_labels[i] for i in rows]}
    for name, mat in matrices.items():
        edges[name] = mat.to_numpy()[rows, cols]
    return pd.DataFrame(edges)


def edge_list(
    data: Coin | pd.DataFrame | np.ndarray,
    procedures: str | Sequence[str] = "Haberman",
    criteria: str = "Z",
    level: float = 0.95,
    bonferroni: bool = False,
    min: float = -np.inf,
    max: float = np.inf,
    support: float = -np.inf,
    directed: bool = False,
    diagonal: bool = False,
    sort: str | None = None,
    decreasing: bool = True,
    pairwise: bool = False,
) -> pd.DataFrame:
    """List of links (with their weights) between pairs of variables.

    Port of ``edgeList`` from the R package netCoin.

    Parameters
    ----------
    data : Coin, DataFrame or ndarray
        A :class:`Coin` object or a table of incidences. With
        ``procedures="shape"`` it is a numeric matrix whose cells become the
        weights; with ``"shape"`` (or ``"tree"``) on a table of text, an
        adjacency list: the first column is linked to the others.
    procedures : str or list of str
        Statistics to include (see :func:`sim`), or ``"shape"``/``"tree"``.
    criteria : str
        Statistic that ``min`` and ``max`` refer to. For p-values (``"Z"`` and
        ``"Fisher"``) ``max`` defaults to 0.5.
    level : float
        Confidence level for ``"confidence"``.
    bonferroni : bool
        Divide a p-value ``max`` by the number of pairs.
    min, max : float
        Keep the pairs whose ``criteria`` lies in ``[min, max]``.
    support : float
        Keep the pairs with at least this many coincidences.
    directed, diagonal : bool
        Include both directions / the diagonal.
    sort : str, optional
        Statistic to order the links by (``decreasing`` sets the direction).
    pairwise : bool
        As in :func:`coin`, for a table of incidences.

    Returns
    -------
    DataFrame
        Columns ``Source``, ``Target`` and one per statistic. Empty (with a
        warning) if no pair meets the conditions. As in R, the list goes column
        by column of the lower triangle, so ``Source`` is the variable that
        comes first.
    """
    first = procedures if isinstance(procedures, str) else list(procedures)[0]
    first = str(first).lower()

    if first.startswith(("sh", "tr")):
        return _shape_edges(data, first, min, max, directed, diagonal, sort, decreasing)

    level = _check_level(level)
    procs = [procedures] if isinstance(procedures, str) else list(procedures)
    codes = list(dict.fromkeys(_method_code(p) for p in procs))
    crit_code = _method_code(criteria)
    sort_code = _method_code(sort) if sort else None

    if crit_code in _PVALUE_CODES and max == np.inf:
        max = 0.50
    if bonferroni and crit_code in _PVALUE_CODES:
        n_vars = len(data) if isinstance(data, Coin) else np.shape(data)[1]
        if n_vars > 1:
            max = max / math.comb(n_vars, 2)

    if isinstance(data, Coin):
        base = data
    else:
        base = coin(data, minimum=1, pairwise=pairwise)

    wanted = list(dict.fromkeys(([sort_code] if sort_code else []) + codes + [crit_code]))
    names = [n for code in wanted for n in _result_names(code)]
    result = sim(base, wanted + ["F"] if "F" not in wanted else wanted, level=level,
                 pairwise=pairwise and base.pairwise)
    matrices = result if isinstance(result, dict) else {_CODE_NAMES[wanted[0]]: result}
    counts = matrices["coincidences"].to_numpy(dtype=float)
    ordered = {n: matrices[n] for n in names}
    crit_name = _CODE_NAMES[crit_code]

    edges = _matrices_to_edges(ordered, counts, crit_name, min, max, support, directed, diagonal)
    if edges.empty:
        warnings.warn("Check max and min values", stacklevel=2)
        return _empty_edges(["Source", "Target"] + names)

    if "c.conditional" in edges:
        edges["c.conditional"] = _as_category(edges["c.conditional"], _CONDITIONAL_LABELS)
    if "c.probable" in edges:
        edges["c.probable"] = _as_category(edges["c.probable"], _PROBABLE_LABELS)
    if sort_code:
        edges = edges.sort_values(_CODE_NAMES[sort_code], ascending=not decreasing,
                                  kind="stable").reset_index(drop=True)
    return edges.rename(columns=_EDGE_RENAME)


def _as_category(values: pd.Series, labels: list[str]) -> pd.Categorical:
    text = values.map(lambda x: labels[int(x)] if pd.notna(x) and 0 <= x <= 8 else None)
    return pd.Categorical(text, categories=labels, ordered=True)


def _shape_edges(data, mode, min_val, max_val, directed, diagonal, sort, decreasing):
    """``procedures="shape"``: a matrix becomes a list of links; text, an adjacency list."""
    if isinstance(data, Coin):
        data = data.matrix
    if not isinstance(data, (pd.DataFrame, np.ndarray)):
        raise TypeError("procedures='shape' needs a matrix or a DataFrame")
    df = pd.DataFrame(data)

    numeric = all(pd.api.types.is_numeric_dtype(t) for t in df.dtypes)
    if mode.startswith("tr") or not numeric:
        rows = []
        for line in df.itertuples(index=False):
            vals = [str(x) for x in line if pd.notna(x) and str(x) != ""]
            rows.extend({"Source": vals[0], "Target": t} for t in vals[1:])
        return pd.DataFrame(rows, columns=["Source", "Target"])

    if min_val == -np.inf:
        min_val = 1
    mat = df.to_numpy(dtype=float)
    row_labels = [str(i) for i in (df.index if isinstance(data, pd.DataFrame) and
                                   not isinstance(df.index, pd.RangeIndex) else range(1, mat.shape[0] + 1))]
    col_labels = [str(c) for c in (df.columns if isinstance(data, pd.DataFrame) and
                                   not isinstance(df.columns, pd.RangeIndex) else range(1, mat.shape[1] + 1))]
    mask = _edge_mask(mat.shape[0], mat.shape[1], row_labels, col_labels, directed, diagonal)
    keep = mask & (mat <= max_val) & (mat >= min_val)
    flat = np.flatnonzero(keep.ravel(order="F"))
    rows, cols = np.unravel_index(flat, keep.shape, order="F")
    out = pd.DataFrame({"Source": [col_labels[j] for j in cols],
                        "Target": [row_labels[i] for i in rows],
                        "value": mat[rows, cols]})
    if out.empty:
        warnings.warn("Check max and min values", stacklevel=3)
        return out
    if sort:
        out = out.sort_values("value", ascending=not decreasing, kind="stable").reset_index(drop=True)
    return out


# -------------------------------------------------------- expected_list ----
def expected_list(
    data: Coin,
    names: Sequence[str] | None = None,
    min: float = 1,
    confidence: bool = False,
) -> pd.DataFrame:
    """Observed and expected coincidences for every pair.

    Port of ``expectedList`` from the R package netCoin. With
    ``confidence=True`` it adds the upper/lower limit of the expectation
    (``confidence``; R uses a fixed 1.64, i.e. a one-sided 95 % level).

    Parameters
    ----------
    data : Coin
        The coincidences.
    names : list of str, optional
        Restrict the pairs to these variables.
    min : float
        Smallest number of coincidences of a pair to be listed.
    confidence : bool
        Add the limits of the expectation.

    Returns
    -------
    DataFrame
        One row per pair (``Source``, ``Target``, ``coincidences``, ``expected``, and the limits with ``confidence``).
    """
    if not isinstance(data, Coin):
        raise TypeError("data must be a Coin object (see coin())")
    labels = list(names) if names is not None else data.labels
    a = data.matrix.to_numpy(dtype=float)
    b = np.diag(a)[:, None] - a
    c = b.T
    n = float(data.n)

    with np.errstate(all="ignore"):
        exp_ = (a + b) * (a + c) / n
        extra = {}
        if confidence:
            sign = 2.0 * (exp_ < a) - 1
            conf = np.maximum(
                exp_ + sign * 1.64 * np.sqrt(exp_ * ((1 - (a + b) / n) * (1 - (a + c) / n))), 0)
            np.fill_diagonal(conf, np.diag(a))
            extra["confidence"] = conf

    mask = _edge_mask(len(labels), len(labels), labels, labels, False, False)
    flat = np.flatnonzero(mask.ravel(order="F"))
    rows, cols = np.unravel_index(flat, mask.shape, order="F")
    out = pd.DataFrame({
        "Source": [labels[j] for j in cols],
        "Target": [labels[i] for i in rows],
        "coincidences": a[rows, cols],
        "expected": exp_[rows, cols],
        **{k: v[rows, cols] for k, v in extra.items()},
    })
    return out[out["coincidences"] >= min].reset_index(drop=True)


# ------------------------------------------------------------ small tools ----
def distant(s, t: bool = False):
    """If ``t`` is true, turn a similarity into a distance (``1 - s``)."""
    return 1.0 - s if t else s


def lower(matrix, decimals: int = 3) -> pd.DataFrame:
    """The lower triangle of a matrix, formatted as text (blank above the diagonal)."""
    if isinstance(matrix, Coin):
        return matrix.lower(decimals)
    if isinstance(matrix, pd.DataFrame):
        labels = list(matrix.index)
        arr = matrix.to_numpy(dtype=float)
    else:
        arr = np.asarray(matrix, dtype=float)
        labels = list(range(arr.shape[0]))
    text = np.empty(arr.shape, dtype=object)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            text[i, j] = "" if j > i else f"{arr[i, j]:.{decimals}f}"
    return pd.DataFrame(text, index=labels, columns=labels)
