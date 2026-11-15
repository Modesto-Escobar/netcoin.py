"""Compare classifications: ``cluster_match`` and ``cluster_agree`` (ports of ``clusterMatch`` and ``clusterAgree``)."""
from __future__ import annotations

import logging
import warnings
from itertools import permutations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from .clusters import as_factor
from .dichotomize import _collation_key

log = logging.getLogger("netcoin")


def optimal_assign(x: np.ndarray) -> np.ndarray:
    """For every row, the column it is given by the assignment that maximizes the total (Hungarian algorithm)."""
    x = np.asarray(x, dtype=float)
    if x.shape[0] != x.shape[1]:
        raise ValueError("optimalAssign: the matrix must be square")
    rows, cols = linear_sum_assignment(x, maximize=True)
    out = np.empty(len(rows), dtype=int)
    out[rows] = cols
    return out


def _tied_assignments(x: np.ndarray, limit: int = 8) -> bool:
    """Whether more than one assignment gives the maximum total (then R's ``clue`` and scipy may differ)."""
    k = x.shape[0]
    if k > limit:
        return False
    best, count = -np.inf, 0
    for p in permutations(range(k)):
        total = x[np.arange(k), list(p)].sum()
        if total > best + 1e-9:
            best, count = total, 1
        elif abs(total - best) <= 1e-9:
            count += 1
    return count > 1


def _factor(values, mask=None) -> tuple[np.ndarray, list[str]]:
    """Codes (0-based) and levels of ``droplevels(as.factor(values))``."""
    values = pd.Series(values)
    if mask is not None:
        values = values[mask]
    factor = as_factor(values)
    used = [c for c in factor.categories if (factor == c).any()]
    codes = np.array([used.index(v) for v in factor.astype(object)], dtype=int)
    return codes, [str(c) for c in used]


def pair_table(pred, ref, weights=None, who: str = "pairTable") -> tuple[pd.DataFrame, int, int]:
    """The contingency table of a pair of classifications (rows: predicted, columns: the model)."""
    pred, ref = pd.Series(pred).reset_index(drop=True), pd.Series(ref).reset_index(drop=True)
    ok = (pred.notna() & ref.notna()).to_numpy().copy()
    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        if len(weights) != len(pred):
            raise ValueError(f"{who}: weights holds {len(weights)} values and the classifications {len(pred)}")
        if (weights[~np.isnan(weights)] < 0).any():
            raise ValueError(f"{who}: there are negative weights")
        ok &= ~np.isnan(weights)
    rows, row_levels = _factor(pred[ok])
    cols, col_levels = _factor(ref[ok])
    table = np.zeros((len(row_levels), len(col_levels)))
    w = np.ones(ok.sum()) if weights is None else weights[ok]
    np.add.at(table, (rows, cols), w)
    if weights is None:
        table = table.astype(np.int64)
    return pd.DataFrame(table, index=row_levels, columns=col_levels), int(ok.sum()), int((~ok).sum())


def class_metrics(m: pd.DataFrame) -> pd.DataFrame:
    """One-against-the-rest measures of a confusion matrix (prediction in rows, reference in columns)."""
    a = m.to_numpy(dtype=float)
    n = a.sum()
    tp = np.diag(a)
    fp = a.sum(axis=1) - tp
    fn = a.sum(axis=0) - tp
    tn = n - tp - fp - fn
    with np.errstate(divide="ignore", invalid="ignore"):
        out = pd.DataFrame({"Sensitivity": tp / (tp + fn), "Specificity": tn / (tn + fp), "Precision": tp / (tp + fp),
                            "F1": 2 * tp / (2 * tp + fp + fn), "Prevalence": a.sum(axis=0) / n, "N": a.sum(axis=0)},
                           index=list(m.columns))
    return out


def confusion_align(pred, ref, weights=None, align: bool = True):
    """Compare one pair of classifications after matching their categories one to one (``confusionAlign``).

    Returns the (realigned) table, the measures by class and the mapping between the categories of the
    model and the predicted ones. R's mapping repeats the names of the model in both columns; here the
    predicted names are the real ones.
    """
    m, rows, dropped = pair_table(pred, ref, weights, who="confusionAlign")
    if m.shape[0] != m.shape[1]:
        raise ValueError(f"confusionAlign: different number of categories ({m.shape[0]} against {m.shape[1]})")
    predicted = list(m.index)
    tied = False
    if align:
        tied = _tied_assignments(m.to_numpy(dtype=float))
        assignment = optimal_assign(m.to_numpy())
        order = np.argsort(assignment, kind="stable")
        m = m.iloc[order]
        predicted = [predicted[i] for i in order]
        m.index = list(m.columns)
    mapping = pd.DataFrame({"Model": list(m.columns), "Predicted": predicted})
    mapping.attrs["ties"] = tied
    return m, class_metrics(m), mapping, float(m.to_numpy().sum()), rows, dropped


def agree_indices(m: pd.DataFrame) -> dict:
    """Adjusted Rand index, normalized mutual information and variation of information (no matching needed)."""
    a = m.to_numpy(dtype=float)
    n = a.sum()
    rs, cs = a.sum(axis=1), a.sum(axis=0)

    def ch2(z):
        return z * (z - 1) / 2

    obs, expected = ch2(a).sum(), ch2(rs).sum() * ch2(cs).sum() / ch2(n)
    maximum = (ch2(rs).sum() + ch2(cs).sum()) / 2
    ari = np.nan if np.isclose(maximum, expected) else (obs - expected) / (maximum - expected)

    def entropy(p):
        p = p[p > 0]
        return -(p * np.log(p)).sum()

    pu, pv, pm = rs / n, cs / n, a / n
    hu, hv = entropy(pu), entropy(pv)
    outer = np.outer(pu, pv)
    keep = pm > 0
    mi = (pm[keep] * np.log(pm[keep] / outer[keep])).sum()
    nmi = np.nan if hu + hv <= 0 else 2 * mi / (hu + hv)
    vi = hu + hv - 2 * mi
    return {"ARI": ari, "NMI": nmi, "VI": vi, "VI.n": vi / np.log(n) if n > 1 else np.nan}


def _is_text(series: pd.Series) -> bool:
    return series.dtype == object or isinstance(series.dtype, pd.CategoricalDtype) or pd.api.types.is_string_dtype(series)


def match_setup(x, model, model_name, columns, weights, modal, same_k, who):
    nodes = x if isinstance(x, pd.DataFrame) else getattr(x, "nodes", None)
    if not isinstance(nodes, pd.DataFrame):
        raise TypeError(f"{who}: x must be a netCoin object or a DataFrame")
    scat = getattr(x, "scat", {}) if not isinstance(x, pd.DataFrame) else {}
    n_nodes = len(nodes)
    idx = scat.get("case_to_pattern")
    n_cases = None if idx is None else len(idx)
    spanning: dict = {}

    def modalize(z, name):
        if idx is None or not _is_text(z):
            return z
        text = z.astype(object)
        count = int(text.map(lambda s: isinstance(s, str) and "|" in s).sum())
        if count:
            spanning[name] = count
            if modal:
                return text.map(lambda s: s.split("|")[0] if isinstance(s, str) else s)
        return z

    if isinstance(model, str) and model in nodes.columns:
        model_name, model_vec, unit = model, modalize(nodes[model], model), "nodes"
    else:
        if isinstance(model, str):
            candidates = list(scat.get("cluster_columns") or [c for c in nodes.columns if _is_text(nodes[c])])
            raise ValueError(f'{who}: there is no node column called "{model}"'
                             + (". Classifications available: " + ", ".join(candidates[:10])
                                + (f" and {len(candidates) - 10} more" if len(candidates) > 10 else "") if candidates else ""))
        model_vec = pd.Series(np.asarray(model))
        if len(model_vec) == n_nodes:
            unit = "nodes"
        elif n_cases is not None and len(model_vec) == n_cases:
            unit = "cases"
        else:
            raise ValueError(f"{who}: the model holds {len(model_vec)} elements; one per node ({n_nodes})"
                             + (f" or one per case ({n_cases})" if n_cases is not None else "") + " was expected.")
    model_codes, model_levels = _factor(pd.Series(model_vec).reset_index(drop=True), pd.Series(model_vec).notna().to_numpy())
    model_series = pd.Series(model_vec).reset_index(drop=True)
    k = len(model_levels)
    if k < 2:
        raise ValueError(f"{who}: the model must hold at least 2 categories")

    case_w = scat.get("case_weight")
    if case_w is not None and n_cases is not None and len(case_w) != n_cases:
        case_w = None
    w_source = "unit"
    if weights is None:
        if unit == "cases":
            w, w_source = (np.asarray(case_w, dtype=float), "caseWeight") if case_w is not None else (None, w_source)
        elif idx is not None:
            cw = np.ones(n_cases) if case_w is None else np.asarray(case_w, dtype=float)
            w = np.bincount(np.asarray(idx) - 1, weights=cw, minlength=n_nodes)
            w_source = "cases per pattern" if case_w is None else "cases per pattern (caseWeight)"
        else:
            w = None
    else:
        if isinstance(weights, str):
            if weights not in nodes.columns:
                raise ValueError(f'{who}: there is no node column called "{weights}"')
            w, w_source = nodes[weights].to_numpy(dtype=float), f'column "{weights}"'
        else:
            w, w_source = np.asarray(weights, dtype=float), "weights"
        n = n_cases if unit == "cases" else n_nodes
        if w.size == 1:
            w = np.full(n, float(w))
        if len(w) != n:
            if unit == "nodes" and n_cases is not None and len(w) == n_cases:
                w = np.bincount(np.asarray(idx) - 1, weights=w, minlength=n_nodes)
                w_source += ", added up by pattern"
            else:
                raise ValueError(f"{who}: weights holds {len(w)} values and the comparison goes by {unit} ({n})")
    if w is not None and (w[~np.isnan(w)] < 0).any():
        raise ValueError(f"{who}: there are negative weights")

    if columns is None:
        columns = scat.get("cluster_columns")
        if columns is None:
            columns = [c for c in nodes.columns if _is_text(nodes[c])]
    columns = [c for c in columns if c in nodes.columns and c != model_name]

    def get_col(name):
        z = modalize(nodes[name], name)
        z = z.reset_index(drop=True)
        return z.iloc[np.asarray(idx) - 1].reset_index(drop=True) if unit == "cases" else z

    n_cat = {c: len(_factor(get_col(c), get_col(c).notna().to_numpy())[1]) for c in columns}
    if same_k:
        columns = [c for c in columns if n_cat[c] == k]
        if not columns:
            raise ValueError(f"{who}: no classification holds {k} categories."
                             + (f" There are {len(spanning)} columns whose patterns span several groups; with "
                                "modal=True the modal group of each pattern is compared." if spanning and not modal else ""))
    elif not columns:
        raise ValueError(f"{who}: there is no classification to compare with")
    return {"nodes": nodes, "model": model_series, "model_name": model_name, "k": k, "unit": unit, "w": w,
            "w_source": w_source, "columns": columns, "n_cat": {c: n_cat[c] for c in columns}, "get_col": get_col,
            "spanning": spanning}


def sort_result(res: pd.DataFrame, sort, who: str) -> pd.DataFrame:
    """``sort``: a column (from the largest down), several (the rest break ties); ``"-col"`` sorts that one up."""
    if sort is None:
        return res
    sort = [sort] if isinstance(sort, str) else list(sort)
    descending = [not s.startswith("-") for s in sort]
    names = [s.lstrip("-") if s.startswith("-") else s for s in sort]
    bad = [n for n in names if n not in res.columns]
    if bad:
        raise ValueError(f"{who}: there is no column " + ", ".join(f'"{b}"' for b in bad) + " to sort by. Available: "
                         + ", ".join(res.columns) + '. The prefix "-" sorts from the smallest up.')
    keys = []
    for name, desc in zip(names, descending):
        col = res[name]
        z = col.to_numpy(dtype=float) if pd.api.types.is_numeric_dtype(col) else np.array(
            [sorted(set(col), key=_collation_key).index(v) for v in col], dtype=float)
        keys.append(-z if desc else z)
    order = np.lexsort(tuple(reversed(keys)))
    return res.iloc[order].reset_index(drop=True)


def report_setup(s, modal: bool, who: str):
    if s["w"] is not None:
        log.info("%s: comparison by %s, weighted by %s (N = %s)", who, s["unit"], s["w_source"], f"{np.nansum(s['w']):g}")
    if s["spanning"]:
        names = ", ".join(f"{k} ({v} patterns)" for k, v in s["spanning"].items())
        warnings.warn(f"{who}: {names} " + ("span" if len(s["spanning"]) > 1 else "spans") + " more than one group. "
                      + ("The modal group of each pattern has been compared, which is all the object keeps: the cases "
                         "of the pattern going to another group count as hits." if modal
                         else "With modal=False they are left out of the comparison."))


def _round_numeric(res: pd.DataFrame, digits: int, keep: set) -> pd.DataFrame:
    for c in res.columns:
        if pd.api.types.is_numeric_dtype(res[c]) and c not in keep:
            res[c] = res[c].round(digits)
    return res


def cluster_match(x, model, columns=None, weights=None, align: bool = True, weighted: bool = True, modal: bool = True,
                  sort=None, digits: int = 3, model_name: str | None = None) -> pd.DataFrame:
    """Compare a model classification with the others of a scattergram, category by category.

    Port of ``clusterMatch``. ``model`` is the name of a node column or a vector (one per node, or per case
    when the nodes are patterns); ``columns`` are the classifications compared (the groups of the object by
    default; only those with as many categories as the model); the categories are matched one to one
    (``align``) so that the diagonal of the table is as large as possible. ``weighted`` adds the means of
    the specificity and F1 weighted by the prevalence; ``modal`` compares, where a pattern spans several
    groups, its modal group. The result has one row per classification; its ``attrs`` hold ``model``,
    ``unit``, ``weights``, ``spanning``, ``tables``, ``byClass``, ``mapping`` and ``ties`` (the
    classifications for which several matchings were equally good: R may pick another of them).

    Parameters
    ----------
    x : netCoin
        A scattergram made by :func:`sur_scat`.
    model : str or sequence
        A node column, or a vector (per node, or per case when the nodes are patterns).
    columns : list of str, optional
        The classifications compared.
    weights : str or sequence, optional
        Weight of each node (a column or a vector).
    align, weighted, modal : bool
        See above.
    sort : str or list of str, optional
        Column(s) to sort by, from the largest down; a ``"-"`` prefix sorts from the smallest up.
    digits : int
        Decimals of the result.
    model_name : str, optional
        Name of the model when it is a vector.

    Returns
    -------
    DataFrame
        One row per classification: ``Classification``, ``k``, ``N``, ``Accuracy``, ``Kappa``, ``Sensitivity``,
        ``Precision``, ``F1`` and ``Specificity`` (macro means; ``.w`` are the weighted ones).
    """
    s = match_setup(x, model, model_name or (model if isinstance(model, str) else "model"), columns, weights, modal,
                    True, "clusterMatch")
    w = s["w"]
    rows, tables, by_class, mapping = [], {}, {}, {}
    for cl in s["columns"]:
        m, bc, mp, n, used, dropped = confusion_align(s["get_col"](cl), s["model"], w, align)
        prevalence = bc["Prevalence"].to_numpy()
        a = m.to_numpy(dtype=float)
        po = np.trace(a) / n
        pe = (a.sum(axis=1) * a.sum(axis=0)).sum() / n ** 2
        vals = {"Classification": cl, "k": s["k"], "N": n, "Accuracy": po, "Kappa": (po - pe) / (1 - pe),
                "Sensitivity": np.nanmean(bc["Sensitivity"]), "Precision": np.nanmean(bc["Precision"]),
                "F1": np.nanmean(bc["F1"]), "Specificity": np.nanmean(bc["Specificity"])}
        if weighted:
            items = list(vals.items())
            out = []
            for key, value in items:
                out.append((key, value))
                if key == "F1":
                    out.append(("F1.w", float(np.nansum(prevalence * bc["F1"].to_numpy()))))
                if key == "Specificity":
                    out.append(("Specificity.w", float(np.nansum(prevalence * bc["Specificity"].to_numpy()))))
            vals = dict(out)
        if w is not None:
            items = list(vals.items())
            vals = dict(items[:2] + [("Rows", used)] + items[2:])
        rows.append(vals)
        tables[cl], by_class[cl], mapping[cl] = m, bc, mp
    res = pd.DataFrame(rows)
    res = sort_result(res, sort, "clusterMatch")
    res = _round_numeric(res, digits, {"k", "Rows"})
    res.attrs.update({"model": s["model_name"], "unit": s["unit"], "weights": w, "spanning": s["spanning"] or None,
                      "tables": tables, "byClass": by_class, "mapping": mapping,
                      "ties": [c for c, mp in mapping.items() if mp.attrs.get("ties")]})
    report_setup(s, modal, "clusterMatch")
    return res


def cluster_agree(x, model, columns=None, weights=None, modal: bool = True, sort=None, digits: int = 3,
                  model_name: str | None = None) -> pd.DataFrame:
    """The same comparison as agreement between partitions (ARI, NMI, VI): no matching, any number of categories.

    Port of ``clusterAgree``. The ``k`` column is that of each classification (the model's is in
    ``attrs['kModel']``).

    Parameters
    ----------
    x : netCoin
        A scattergram made by :func:`sur_scat`.
    model : str or sequence
        A node column, or a vector (per node, or per case when the nodes are patterns).
    columns, weights, modal, sort, digits, model_name
        As in :func:`cluster_match`.

    Returns
    -------
    DataFrame
        One row per classification: ``Classification``, ``k``, ``N``, ``ARI``, ``NMI``, ``VI`` and ``VI.n``
        (the variation of information, raw and normalized).
    """
    s = match_setup(x, model, model_name or (model if isinstance(model, str) else "model"), columns, weights, modal,
                    False, "clusterAgree")
    w = s["w"]
    rows, tables = [], {}
    for cl in s["columns"]:
        m, used, dropped = pair_table(s["get_col"](cl), s["model"], w, who="clusterAgree")
        vals = {"Classification": cl, "k": s["n_cat"][cl], "N": float(m.to_numpy().sum()), **agree_indices(m)}
        if w is not None:
            items = list(vals.items())
            vals = dict(items[:2] + [("Rows", used)] + items[2:])
        rows.append(vals)
        tables[cl] = m
    res = pd.DataFrame(rows)
    res = sort_result(res, sort, "clusterAgree")
    res = _round_numeric(res, digits, {"k", "Rows"})
    res.attrs.update({"model": s["model_name"], "kModel": s["k"], "unit": s["unit"], "weights": w,
                      "spanning": s["spanning"] or None, "tables": tables})
    report_setup(s, modal, "clusterAgree")
    return res
