"""Regression networks: ``glm_coin`` (port of ``glmCoin``), with its own marginal effects.

The links of ``glmCoin`` are the effects of the categories of the explanatory variables on the
response, each against the weighted mean of its variable (and, besides, against the first
category), and the average slopes of the numeric ones. R computes them with the package
``marginaleffects``; here the average predictions, their contrasts and their standard errors
(delta method, model-based or HC1 covariance) are computed directly.
"""
from __future__ import annotations

import re
import unicodedata
import warnings
from dataclasses import dataclass, field
from itertools import product
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .core import _normalize_legacy, name_by_language, netCoin
from .dichotomize import _collation_key
from .utils import getByLanguage, nameList

_FAMILIES = ("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian", "quasi", "quasibinomial", "quasipoisson")
_Z = stats.norm.ppf(0.975)


def _require_models():
    try:
        import patsy                                                      # noqa: F401
        import statsmodels.api as sm
    except ImportError as exc:                                             # pragma: no cover
        raise ImportError("glm_coin needs statsmodels and patsy: pip install statsmodels") from exc
    return sm


# ------------------------------------------------------------------ parsing ----
def split_top_level(text: str, separators: str = ",\n") -> tuple[str, str | None]:
    """Split ``text`` at the first separator that is not inside brackets or quotes."""
    depth, quote, escape = 0, None, False
    for i, ch in enumerate(text):
        if escape:
            escape = False
            continue
        if ch == "\\":
            escape = True
            continue
        if quote:
            if ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
        elif ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth = max(0, depth - 1)
        elif depth == 0 and ch in separators:
            rest = text[i + 1:].strip()
            return text[:i].strip(), rest or None
    return text.strip(), None


def parse_specs(spec: str | Sequence[str]) -> list[tuple[str, str, str | None]]:
    """The equations of a specification: ``(line, formula, family)`` for each line."""
    lines = spec.split("\n") if isinstance(spec, str) else list(spec)
    lines = [l.strip() for l in lines if l.strip() and not l.strip().startswith("#")]
    out = []
    for line in lines:
        formula, family = split_top_level(line, ",")
        if "~" not in formula:
            raise ValueError(f"It could not be converted into a formula: '{formula}'")
        out.append((line, formula, _check_family(family)))
    return out


def _check_family(text: str | None) -> str | None:
    if text is None:
        return None
    text = re.sub(r"^family\s*=\s*", "", text.strip())
    base = re.sub(r"\s*\(.*$", "", text)
    if base == "multinomial":
        return "multinomial"
    if base not in _FAMILIES:
        raise ValueError(f"Family not allowed: '{base}'. Allowed: {', '.join(_FAMILIES)}, multinomial")
    return base


def clean_variables(names: Sequence[str]) -> list[str]:
    """``cleanVariables``: ASCII, spaces and anything but letters, digits and ``_`` become dots."""
    out = []
    for n in names:
        n = unicodedata.normalize("NFKD", str(n)).encode("ascii", "ignore").decode("ascii")
        out.append(re.sub(r"[^a-zA-Z0-9_]", ".", n.replace(" ", ".")))
    return out


# -------------------------------------------------------------------- models ----
@dataclass
class Model:
    """A fitted equation: what the contrasts and the summaries need."""

    name: str
    formula: str
    family: str                # "gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian", "quasi...", "multinomial"
    declared: str              # the family as stated in the equation ("gaussian" when it is lm)
    is_lm: bool
    target: str
    response: pd.Series
    frame: pd.DataFrame        # the cases used, with every variable of the formula
    weights: np.ndarray        # prior weights of those cases
    design_info: object
    X: np.ndarray
    coef: np.ndarray
    names: list[str]
    cov: np.ndarray            # covariance of the coefficients (model-based or robust)
    dispersion: float
    sm_result: object = None
    factors: dict = field(default_factory=dict)       # variable -> levels
    classes: list | None = None                       # multinomial: categories of the response
    y: np.ndarray | None = None
    cov_model: np.ndarray | None = None               # covariance of the model, even when cov is robust
    column_order: list | None = None                  # patsy's columns in the order of R's model.matrix
    terms: list | None = None                         # the terms in the order of R's term.labels

    @property
    def link(self) -> str:
        return {"gaussian": "identity", "binomial": "logit", "quasibinomial": "logit", "poisson": "log",
                "quasipoisson": "log", "Gamma": "inverse", "inverse.gaussian": "1/mu^2",
                "multinomial": "logit"}[self.family]

    @property
    def is_binomial(self) -> bool:
        return "binomial" in self.family

    @property
    def is_multinomial(self) -> bool:
        return self.family == "multinomial"


def _sm_family(family: str):
    sm = _require_models()
    return {"gaussian": sm.families.Gaussian(), "binomial": sm.families.Binomial(),
            "quasibinomial": sm.families.Binomial(), "poisson": sm.families.Poisson(),
            "quasipoisson": sm.families.Poisson(), "Gamma": sm.families.Gamma(),
            "inverse.gaussian": sm.families.InverseGaussian()}[family]


def linkinv(family: str, eta: np.ndarray) -> np.ndarray:
    """Inverse of the link of ``family``."""
    if family in ("binomial", "quasibinomial", "multinomial"):
        return 1 / (1 + np.exp(-eta))
    if family in ("poisson", "quasipoisson"):
        return np.exp(eta)
    if family == "Gamma":
        return 1 / eta
    if family == "inverse.gaussian":
        return 1 / np.sqrt(eta)
    return eta


def _mean_variance(family: str, mu: np.ndarray) -> np.ndarray:
    if family in ("binomial", "quasibinomial"):
        return mu * (1 - mu)
    if family in ("poisson", "quasipoisson"):
        return mu
    if family == "Gamma":
        return mu ** 2
    if family == "inverse.gaussian":
        return mu ** 3
    return np.ones_like(mu)


def _dmu_deta(family: str, eta: np.ndarray, mu: np.ndarray) -> np.ndarray:
    if family in ("binomial", "quasibinomial"):
        return mu * (1 - mu)
    if family in ("poisson", "quasipoisson"):
        return mu
    if family == "Gamma":
        return -(mu ** 2)
    if family == "inverse.gaussian":
        return -0.5 * mu ** 3
    return np.ones_like(mu)


def _plain(name: str) -> str:
    """``Q("a.b")[T.x]`` as ``a.b[T.x]``."""
    return re.sub(r"Q\([\"'](.+?)[\"']\)", lambda m: m.group(1), name)


def _levels(series: pd.Series) -> list[str]:
    if isinstance(series.dtype, pd.CategoricalDtype):
        present = set(series.dropna().astype(str))
        return [str(c) for c in series.cat.categories if str(c) in present]
    from .dichotomize import _as_character
    return sorted(set(_as_character(series.dropna())), key=_collation_key)


def _is_factor(series: pd.Series) -> bool:
    return (isinstance(series.dtype, pd.CategoricalDtype) or series.dtype == object
            or pd.api.types.is_bool_dtype(series) or pd.api.types.is_string_dtype(series))


def fit_model(line: str, formula: str, declared: str | None, data: pd.DataFrame, weights: np.ndarray | None,
              robust: bool, auto_quasi: bool = True, default_family: str = "gaussian") -> Model:
    """Fit one equation (``lm`` or ``glm``) and keep what the contrasts need."""
    sm = _require_models()
    import patsy

    lhs, rhs = [p.strip() for p in formula.split("~", 1)]
    family = declared or default_family
    is_lm = declared is None and default_family == "gaussian"
    if family == "multinomial":
        raise NotImplementedError("multinomial equations are not available yet")
    if family == "quasi":
        raise NotImplementedError("the 'quasi' family needs its link and variance; use quasibinomial or quasipoisson")
    if auto_quasi and family == "binomial":
        family = "quasibinomial"
    names = [c for c in data.columns if re.search(rf"(?<![\w.]){re.escape(c)}(?![\w.])", formula)]
    work = data[names].copy()
    if weights is not None:
        work["(weights)"] = weights
    work = work.dropna()
    for c in names:
        if _is_factor(work[c]):
            levels = _levels(work[c])
            work[c] = pd.Categorical(work[c].astype(str), categories=levels)
    w = np.ones(len(work)) if weights is None else work["(weights)"].to_numpy(dtype=float)
    for c in sorted(names, key=len, reverse=True):          # patsy reads a dot as an attribute: quote such names
        if not c.isidentifier():
            rhs = re.sub(rf"(?<![\w.]){re.escape(c)}(?![\w.])", lambda m, c=c: f'Q("{c}")', rhs)
    X = patsy.dmatrix(rhs, work, return_type="dataframe")
    info = X.design_info
    # R lists the terms by degree and, within a degree, as they appear in the formula
    desc = patsy.ModelDesc.from_formula("y ~ " + rhs)
    term_list = sorted([t for t in desc.rhs_termlist if t.name() != "Intercept"], key=lambda t: len(t.factors))
    slices = info.term_name_slices
    order = (list(range(*slices["Intercept"].indices(X.shape[1]))) if "Intercept" in slices else [])
    for t in term_list:
        order += list(range(*slices[t.name()].indices(X.shape[1])))
    X = X.iloc[:, order]
    X.columns = [_plain(c) for c in X.columns]
    terms = [_plain(t.name()) for t in term_list]
    y_raw = work[lhs] if lhs in work.columns else patsy.dmatrix("0 + " + lhs, work, return_type="dataframe").iloc[:, 0]
    factors = {c: _levels(work[c]) for c in names if isinstance(work[c].dtype, pd.CategoricalDtype)}
    if "binomial" in family:
        if isinstance(y_raw.dtype, pd.CategoricalDtype):
            y = (y_raw.astype(str) != factors[lhs][0]).astype(float).to_numpy()
        elif pd.api.types.is_bool_dtype(y_raw):
            y = y_raw.astype(float).to_numpy()
        else:
            y = y_raw.to_numpy(dtype=float)
    else:
        y = y_raw.to_numpy(dtype=float)
    scale = None if family in ("gaussian", "Gamma", "inverse.gaussian", "quasibinomial", "quasipoisson") else 1.0
    model = sm.GLM(y, X.to_numpy(), family=_sm_family(family), var_weights=w)
    result = model.fit(scale="X2" if scale is None else scale, maxiter=100, tol=1e-12)
    coef = np.asarray(result.params)
    dispersion = float(result.scale)
    if robust:
        cov = hc1(family, X.to_numpy(), y, coef, w)
    else:
        cov = np.asarray(result.cov_params())
    return Model(line, formula, family, declared or "gaussian", is_lm, lhs, y_raw, work.drop(columns="(weights)", errors="ignore"),
                 w, info, X.to_numpy(), coef, list(X.columns), cov, dispersion, result, factors, None, y,
                 np.asarray(result.cov_params()), order, terms)


def hc1(family: str, X: np.ndarray, y: np.ndarray, coef: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Heteroskedasticity-consistent covariance of the coefficients, type HC1 (as ``sandwich::vcovHC``)."""
    eta = X @ coef
    mu = linkinv(family, eta)
    dmu = _dmu_deta(family, eta, mu)
    var = _mean_variance(family, mu)
    working_weight = w * dmu ** 2 / var
    bread = np.linalg.inv(X.T @ (X * working_weight[:, None]))
    score = X * (w * (y - mu) * dmu / var)[:, None]
    n, k = X.shape
    return n / (n - k) * bread @ (score.T @ score) @ bread


# --------------------------------------------------------- marginal effects ----
def _design_for(model: Model, frame: pd.DataFrame) -> np.ndarray:
    import patsy
    built = patsy.build_design_matrices([model.design_info], frame, return_type="dataframe")[0]
    return np.asarray(built)[:, model.column_order]


def _jacobian(function, beta: np.ndarray) -> np.ndarray:
    """Central-difference Jacobian of ``function(beta)`` (a vector) with respect to ``beta``."""
    base = np.atleast_1d(function(beta))
    jac = np.empty((len(base), len(beta)))
    for j in range(len(beta)):
        h = 1e-6 * max(1.0, abs(beta[j]))
        up, down = beta.copy(), beta.copy()
        up[j] += h
        down[j] -= h
        jac[:, j] = (np.atleast_1d(function(up)) - np.atleast_1d(function(down))) / (2 * h)
    return jac


def _inference(estimate: np.ndarray, se: np.ndarray) -> dict:
    with np.errstate(divide="ignore", invalid="ignore"):
        z = estimate / se
    p = 2 * stats.norm.sf(np.abs(z))
    with np.errstate(divide="ignore"):
        s = -np.log2(p)
    return {"estimate": estimate, "std.error": se, "statistic": z, "p.value": p, "s.value": s,
            "conf.low": estimate - _Z * se, "conf.high": estimate + _Z * se}


def contrasts_gw(model: Model, group_vars: list[str]) -> pd.DataFrame:
    """Contrasts of the average predictions of a (combination of) factor(s) against their weighted mean.

    Port of ``contrastes_gw``. The average prediction of every combination of categories is that of
    all the cases, whatever their own category, with it set to the combination; the contrast
    subtracts the mean of those predictions weighted by the share of every combination in the
    sample. The columns ``*.base`` contrast them with the first combination.
    """
    frame = model.frame
    for v in group_vars:
        if v not in model.factors:
            raise ValueError(f"The variable '{v}' is not a factor.")
        if len(model.factors[v]) < 2:
            raise ValueError(f"The variable '{v}' has only one observed level.")
    combos = list(product(*[model.factors[v] for v in group_vars]))
    w = model.weights

    def averages(beta):
        out = []
        for combo in combos:
            counterfactual = frame.copy()
            for v, level in zip(group_vars, combo):
                counterfactual[v] = pd.Categorical([level] * len(frame), categories=model.factors[v])
            eta = _design_for(model, counterfactual) @ beta
            out.append((w * linkinv(model.family, eta)).sum() / w.sum())
        return np.array(out)

    estimates = averages(model.coef)
    jac = _jacobian(averages, model.coef)
    sigma = jac @ model.cov @ jac.T

    shares = frame[group_vars].astype(str).assign(w=w).groupby(group_vars, sort=False)["w"].sum()
    lookup = {k if isinstance(k, tuple) else (k,): v for k, v in shares.items()}
    try:
        prop = np.array([lookup[c] for c in combos], dtype=float)
    except KeyError:
        raise ValueError("Weights could not be aligned with 'margin' levels")
    prop = prop / lookup_sum(shares)
    k = len(combos)
    transform = np.eye(k) - np.outer(np.ones(k), prop / prop.sum())
    result = _inference(transform @ estimates, np.sqrt(np.diag(transform @ sigma @ transform.T)))

    is_base = np.array([all(level == model.factors[v][0] for v, level in zip(group_vars, combo)) for combo in combos])
    base = int(np.where(is_base)[0][0]) if is_base.any() else 0
    to_base = np.eye(k) - np.outer(np.ones(k), np.eye(k)[base])
    diff = to_base @ estimates
    se_diff = np.sqrt(np.diag(to_base @ sigma @ to_base.T))
    with np.errstate(divide="ignore", invalid="ignore"):
        p_base = 2 * stats.norm.sf(np.abs(diff / se_diff))
    se_diff = np.where(is_base, np.nan, se_diff)
    p_base = np.where(is_base, np.nan, p_base)
    result.update({"estimate.base": diff, "std.error.base": se_diff, "p.value.base": p_base})
    out = pd.DataFrame(result)
    if model.is_binomial:
        for col in ("estimate", "std.error", "conf.low", "conf.high", "estimate.base", "std.error.base"):
            out[col] = out[col] * 100
    label = "*".join(group_vars)
    out.insert(0, "Source", [f"{label}:{'*'.join(c)}" for c in combos])
    out.insert(1, "Target", model.target)
    return out


def lookup_sum(shares: pd.Series) -> float:
    return float(shares.sum())


def avg_slope(model: Model, variable: str) -> pd.DataFrame:
    """The average slope of the prediction on a numeric variable (``marginaleffects::avg_slopes``)."""
    frame, w = model.frame, model.weights
    x = frame[variable].to_numpy(dtype=float)
    h = 1e-5 * max(1.0, float(np.abs(x).max()))

    def slope(beta):
        up, down = frame.copy(), frame.copy()
        up[variable], down[variable] = x + h, x - h
        p_up = linkinv(model.family, _design_for(model, up) @ beta)
        p_down = linkinv(model.family, _design_for(model, down) @ beta)
        return np.array([(w * (p_up - p_down) / (2 * h)).sum() / w.sum()])

    estimate = slope(model.coef)
    jac = _jacobian(slope, model.coef)
    se = np.sqrt(np.diag(jac @ model.cov @ jac.T))
    out = pd.DataFrame(_inference(estimate, se))
    if model.is_binomial:
        for col in ("estimate", "std.error", "conf.low", "conf.high"):
            out[col] = out[col] * 100
    out["estimate.base"], out["std.error.base"], out["p.value.base"] = out["estimate"], out["std.error"], out["p.value"]
    return out


def _wmean(x, w):
    ok = ~np.isnan(x)
    return float((w[ok] * x[ok]).sum() / w[ok].sum())


def summary_nodes_links(model: Model, drop_empty: bool = True, outliers: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Nodes and links of one equation (``contrastes_gw_one``)."""
    frame, w = model.frame, model.weights
    w_total = float(w.sum())
    target = model.target
    y = model.response
    pred_cat = None
    if model.is_binomial:
        if target in model.factors:
            pred_cat = model.factors[target][1]
        elif pd.api.types.is_bool_dtype(y):
            pred_cat = "TRUE"
        else:
            u = sorted(pd.Series(y).dropna().unique())
            pred_cat = _text_number(u[1]) if len(u) > 1 else None
    target_label = f"{target}:{pred_cat}" if pred_cat is not None else target
    terms = model.terms

    links, nodes = [], []
    for term in terms:
        parts = term.split(":")
        if len(parts) > 1:
            if not all(p in model.factors for p in parts):
                continue
            if drop_empty and any(len(model.factors[p]) < 2 for p in parts):
                continue
            try:
                table = contrasts_gw(model, parts)
            except ValueError as exc:
                warnings.warn(f"GW interaction failed '{term}': {exc}")
                continue
            if model.is_binomial and pred_cat is not None:
                table["Target"] = target_label
            links.append(table)
            agg = frame[parts].astype(str).assign(Frequency=w).groupby(parts, sort=False)["Frequency"].sum().reset_index()
            agg = _aggregate_order(agg, parts, model)
            label = "*".join(parts)
            nodes.append(pd.DataFrame({"name": [f"{label}:" + "*".join(str(r[p]) for p in parts) for _, r in agg.iterrows()],
                                       "variable": label, "Frequency": agg["Frequency"].to_numpy(),
                                       "%": agg["Frequency"].to_numpy() / w_total * 100, "Type": "factor_interaction",
                                       "Mean": np.nan, "Min": np.nan, "Max": np.nan}))
            continue
        if term not in frame.columns:
            continue
        if term in model.factors:
            if drop_empty and len(model.factors[term]) < 2:
                continue
            try:
                table = contrasts_gw(model, [term])
            except ValueError as exc:
                warnings.warn(f"GW factor failed '{term}': {exc}")
                continue
            if model.is_binomial and pred_cat is not None:
                table["Target"] = target_label
            links.append(table)
            sums = frame[term].astype(str).groupby(frame[term].astype(str)).size()
            freq = pd.Series(w, index=frame.index).groupby(frame[term].astype(str)).sum()
            levels = model.factors[term]
            nodes.append(pd.DataFrame({"name": [f"{term}:{l}" for l in levels], "variable": term,
                                       "Frequency": [freq.get(l, 0.0) for l in levels],
                                       "%": [freq.get(l, 0.0) / w_total * 100 for l in levels], "Type": "factor",
                                       "Mean": np.nan, "Min": np.nan, "Max": np.nan}))
        elif pd.api.types.is_numeric_dtype(frame[term]):
            table = avg_slope(model, term)
            table["Source"], table["Target"] = term, target_label
            links.append(table)
            x = frame[term].to_numpy(dtype=float)
            nodes.append(_numeric_node(term, x, w, "numeric"))

    # the response
    if target in frame.columns or target in model.factors or True:
        y_values = model.response
        if target in model.factors:
            levels = model.factors[target]
            freq = pd.Series(w, index=frame.index).groupby(y_values.astype(str)).sum()
            keep = [pred_cat] if model.is_binomial and pred_cat is not None else levels
            node = pd.DataFrame({"name": [f"{target}:{l}" for l in keep], "variable": target,
                                 "Frequency": [freq.get(l, 0.0) for l in keep],
                                 "%": [freq.get(l, 0.0) / w_total * 100 for l in keep], "Type": "response_factor",
                                 "Mean": np.nan, "Min": np.nan, "Max": np.nan})
            nodes.append(node)
        elif pd.api.types.is_numeric_dtype(y_values):
            x = y_values.to_numpy(dtype=float)
            u = sorted(np.unique(x[~np.isnan(x)]))
            if model.is_binomial and len(u) == 2 and set(u) <= {0.0, 1.0}:
                freq1 = float(w[x == 1].sum())
                nodes.append(pd.DataFrame({"name": [f"{target}:1"], "variable": target, "Frequency": [freq1],
                                           "%": [freq1 / w_total * 100], "Type": "response_factor",
                                           "Mean": np.nan, "Min": np.nan, "Max": np.nan}))
            else:
                nodes.append(_numeric_node(target, x, w, "response_numeric"))
        elif pd.api.types.is_bool_dtype(y_values) and model.is_binomial:
            freq1 = float(w[y_values.to_numpy()].sum())
            nodes.append(pd.DataFrame({"name": [f"{target}:TRUE"], "variable": target, "Frequency": [freq1],
                                       "%": [freq1 / w_total * 100], "Type": "response_factor",
                                       "Mean": np.nan, "Min": np.nan, "Max": np.nan}))
    link_table = _stack(links)
    if len(link_table):
        link_table = link_table[["Source", "Target"] + [c for c in link_table.columns if c not in ("Source", "Target")]]
    else:
        link_table = pd.DataFrame({"Source": pd.Series(dtype=str), "Target": pd.Series(dtype=str)})
    node_table = _stack(nodes)
    if len(node_table):
        mm = ["Mean", "Min", "Max"]
        if node_table[mm].isna().all().all():
            node_table = node_table.drop(columns=mm)
    return node_table, link_table


def _text_number(value) -> str:
    value = float(value)
    return str(int(value)) if value.is_integer() else f"{value:g}"


def _numeric_node(name: str, x: np.ndarray, w: np.ndarray, kind: str) -> pd.DataFrame:
    finite = x[np.isfinite(x)]
    mean = _wmean(x, w)
    lo, hi = (float(finite.min()), float(finite.max())) if len(finite) else (np.nan, np.nan)
    norm = np.nan if (not np.isfinite(hi - lo) or hi <= lo) else (mean - lo) / (hi - lo) * 100
    return pd.DataFrame({"name": [name], "variable": name, "Frequency": [float(w.sum())], "%": [norm], "Type": kind,
                         "Mean": [mean], "Min": [lo], "Max": [hi]})


def _aggregate_order(agg: pd.DataFrame, parts: list[str], model: Model) -> pd.DataFrame:
    """The order of R's ``aggregate``: the first variable varies fastest, by the order of the levels."""
    keys = [agg[p].map({l: i for i, l in enumerate(model.factors[p])}).to_numpy() for p in parts]
    return agg.iloc[np.lexsort(tuple(keys))].reset_index(drop=True)


def _stack(frames: list[pd.DataFrame]) -> pd.DataFrame:
    frames = [f for f in frames if len(f)]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


# --------------------------------------------------------------------- glm_coin ----
def contr_gw(formulas, data: pd.DataFrame, weights=None, robust: bool = False, drop_empty: bool = True,
             outliers: bool = True, default_family: str = "gaussian", auto_quasi: bool = True) -> dict:
    """Nodes, links and models of a set of equations (``contr.gw`` of R, for specifications)."""
    specs = parse_specs(formulas)
    w = None if weights is None else np.asarray(weights, dtype=float)
    if w is not None and len(w) != len(data):
        raise ValueError("The length of 'weights' does not match 'nrow(data)'.")
    all_nodes, all_links, models = [], [], {}
    for line, formula, family in specs:
        model = fit_model(line, formula, family, data, w, robust, auto_quasi, default_family)
        nodes, links = summary_nodes_links(model, drop_empty, outliers)
        if len(links):
            links["Model"] = line
            all_links.append(links)
        if len(nodes):
            all_nodes.append(nodes)
        models[line] = model
    links = _stack(all_links)
    nodes = None
    for part in all_nodes:
        if nodes is None:
            nodes = part
        else:
            new = part[~part["name"].isin(nodes["name"])]
            nodes = pd.concat([nodes, new], ignore_index=True, sort=False)
    return {"Nodes": nodes, "Links": links, "models": models}


def _restore(text: str, mapping: dict[str, str]) -> str:
    """Labels are ``variable[*variable][:categories]``: the variables go back to their names."""
    head, sep, rest = text.partition(":")
    tokens = [mapping.get(t, t) for t in head.split("*")]
    return "*".join(tokens) + sep + rest


def glm_coin(
    formulas: str | Sequence[str],
    data: pd.DataFrame,
    weights: str | Sequence[float] | None = None,
    pmax: float = 0.05,
    robust: bool = False,
    show_arrows: bool = True,
    color: str = "variable",
    lwidth: str = "s.value",
    lcolor: str = "estimate",
    language: str = "en",
    igraph: bool = False,
    **netcoin_args,
):
    """A network of the effects in one or several regression equations (port of ``glmCoin``).

    Each line of ``formulas`` is an equation, ``"response ~ term + term, family"``, where the family
    is ``gaussian`` (the default, a linear model), ``binomial`` (fitted as quasibinomial, with
    estimated dispersion), ``poisson``, ``quasipoisson``, ``Gamma``, ``inverse.gaussian`` or
    ``quasibinomial``. The categories of every factor are nodes, and so are the numeric variables
    and the response; the links go from every explanatory category to the response, with the
    effect of the category against the mean of its variable weighted by the sample (``estimate``)
    and against its first category (``estimate.base``), and from every numeric variable with its
    average slope. Effects in binomial equations are in percentage points. The standard errors
    come from the delta method with the covariance of the model, or the robust HC1 one.

    Parameters
    ----------
    formulas : str or list of str
        The equations (variable names with spaces and other signs are accepted).
    data : DataFrame
        Text and categorical columns are the factors, numeric columns the quantitative variables.
    weights : str or sequence, optional
        The weights of the cases, or the column that holds them.
    pmax : float
        Links whose p-value is not below it are hidden.
    robust : bool
        Use the HC1 covariance.
    show_arrows, color, lwidth, lcolor, language
        As in :class:`netCoin`.
    igraph : bool
        Return an ``igraph.Graph``.
    **netcoin_args
        The rest goes to :class:`netCoin`.

    Returns
    -------
    netCoin
        With the attribute ``equations`` (fit of every equation, see :func:`eq_summary`) and ``coeffs``
        (their coefficients, see :func:`eq_coefs`).
    """
    args = _normalize_legacy(netcoin_args, strict=False)
    original = list(map(str, data.columns))
    cleaned = clean_variables(original)
    data = data.copy()
    data.columns = cleaned
    changed = [i for i in range(len(original)) if original[i] != cleaned[i]]
    changed.sort(key=lambda i: -len(original[i]))
    text = formulas if isinstance(formulas, str) else "\n".join(formulas)
    for i in changed:
        text = text.replace(original[i], cleaned[i])
    if isinstance(weights, str):
        weights = clean_variables([weights])[0]
        if weights not in data.columns:
            raise ValueError(f"Weights variable '{weights}' not found in data")
        weights = data[weights].to_numpy(dtype=float)

    built = contr_gw(text, data, weights=weights, robust=robust)
    nodes, links = built["Nodes"], built["Links"]
    models = list(built["models"].values())
    equations = eq_summary(models)
    coeffs = eq_coefs(models)
    args.update(show_arrows=show_arrows, color=color, link_filter=f"`p.value` < {pmax}", lwidth=lwidth,
                lcolor=lcolor, language=language)
    args.setdefault("link_bipolar", True)
    back = {cleaned[i]: original[i] for i in changed}
    if changed:
        if nodes is not None:
            nodes["name"] = [_restore(s, back) for s in nodes["name"]]
            if "variable" in nodes:
                nodes["variable"] = [_restore(s, back) for s in nodes["variable"]]
        if len(links):
            links["Source"] = [_restore(s, back) for s in links["Source"]]
            links["Target"] = [_restore(s, back) for s in links["Target"]]
            links["Model"] = [_restore_formula(s, back) for s in links["Model"]]
        equations["model"] = [_restore_formula(s, back) for s in equations["model"]]
        for table in coeffs:
            if table is not None:
                table["Source"] = [_restore(s, back) for s in table["Source"]]
                table["Target"] = [_restore(s, back) for s in table["Target"]]
                table.attrs["latent.sd"] = {_restore(k, back): v for k, v in table.attrs["latent.sd"].items()}
    if language != "en":
        nodes = nodes.rename(columns={"name": getByLanguage(nameList, language)})
        args["name"] = getByLanguage(nameList, language)
    args["nodes"] = nodes
    args["links"] = links if len(links) else None
    net = netCoin(**args)
    net.equations, net.coeffs = equations, coeffs
    if igraph:
        from .igraph_io import to_igraph

        graph = to_igraph(net)
        graph["equations"], graph["coeffs"] = equations, coeffs
        return graph
    return net


def _restore_formula(text: str, mapping: dict[str, str]) -> str:
    for clean, original in sorted(mapping.items(), key=lambda kv: -len(kv[0])):
        text = re.sub(rf"(?<![a-zA-Z0-9._]){re.escape(clean)}(?![a-zA-Z0-9._])", lambda m: original, text)
    return text


# ------------------------------------------------------------ fit and coefficients ----
def _deviance(family: str, y: np.ndarray, mu: np.ndarray, w: np.ndarray) -> float:
    with np.errstate(divide="ignore", invalid="ignore"):
        if family in ("binomial", "quasibinomial"):
            d = np.where(y > 0, y * np.log(y / mu), 0.0) + np.where(y < 1, (1 - y) * np.log((1 - y) / (1 - mu)), 0.0)
            return float(2 * (w * d).sum())
        if family in ("poisson", "quasipoisson"):
            return float(2 * (w * (np.where(y > 0, y * np.log(y / mu), 0.0) - (y - mu))).sum())
        if family == "Gamma":
            return float(2 * (w * (-np.log(y / mu) + (y - mu) / mu)).sum())
        if family == "inverse.gaussian":
            return float((w * (y - mu) ** 2 / (y * mu ** 2)).sum())
    return float((w * (y - mu) ** 2).sum())


def _family_aic(family: str, y: np.ndarray, mu: np.ndarray, w: np.ndarray, dev: float) -> float:
    """The ``aic`` of the family functions of R (NaN for the quasi ones)."""
    from scipy.special import gammaln
    if family == "gaussian":
        n = len(y)
        return float(n * (np.log(2 * np.pi * dev / n) + 1) + 2 - np.log(w[w > 0]).sum())
    if family == "poisson":
        return float(-2 * (w * (y * np.log(mu) - mu - gammaln(y + 1))).sum())
    if family == "Gamma":
        n = w.sum()
        disp = dev / n
        shape = 1 / disp
        logd = (shape - 1) * np.log(y) - y * shape / mu - shape * np.log(mu * disp) - gammaln(shape)
        return float(-2 * (logd * w).sum() + 2)
    if family == "inverse.gaussian":
        return float(w.sum() * (np.log(dev / w.sum() * 2 * np.pi) + 1) + 3 * (np.log(y) * w).sum() + 2)
    return np.nan


def eq_summary(models: list[Model]) -> pd.DataFrame:
    """Fit measures of every equation: pseudo R-squared, log-likelihood, AIC, BIC, deviances (``eqSummary``)."""
    rows = []
    for m in models:
        w, y = m.weights, m.y
        nobs = int((w > 0).sum())
        n = nobs if np.allclose(w, 1) else float(w.sum())
        eta = m.X @ m.coef
        mu = linkinv(m.family, eta)
        dev = _deviance(m.family, y, mu, w)
        mu0 = np.full_like(y, (w * y).sum() / w.sum())
        dev0 = _deviance(m.family, y, mu0, w)
        rank = m.X.shape[1]
        aic_family = _family_aic(m.family, y, mu, w, dev)
        aic_family0 = _family_aic(m.family, y, mu0, w, dev0)
        extra = 1 if m.family in ("gaussian", "Gamma", "inverse.gaussian") else 0
        if np.isnan(aic_family):
            ll = ll0 = aic = bic = np.nan
        else:
            aic = aic_family + 2 * rank
            ll = rank + extra - aic / 2
            ll0 = 1 + extra - (aic_family0 + 2) / 2
            bic = -2 * ll + np.log(nobs) * (rank + extra)
        r2 = 1 - dev / dev0
        mcfadden = 1 - ll / ll0
        if not np.isnan(ll) and not np.isnan(ll0):
            nagelkerke = (1 - np.exp(2 * (ll0 - ll) / n)) / (1 - np.exp(2 * ll0 / n))
        else:
            nagelkerke = (1 - np.exp((dev - dev0) / n)) / (1 - np.exp(-dev0 / n))
        rows.append({"model": m.name.split(",")[0].strip() if "," in m.name else m.name, "family": m.declared,
                     "n": round(n), "n_obs": nobs, "r2": round(r2, 3), "r2.mcfadden": round(mcfadden, 3),
                     "r2.nagelkerke": round(nagelkerke, 3), "logLik": round(ll, 2), "aic": round(aic, 2),
                     "bic": round(bic, 2), "deviance": round(dev, 2), "null.deviance": round(dev0, 2),
                     "df.residual": nobs - rank})
    return pd.DataFrame(rows)


def _r_column_name(column: str) -> tuple[str, list[tuple[str, str | None]]]:
    """``Gender[T.Male]:edu`` as R names it (``GenderMale:edu``) and its parts ``[(variable, level)]``."""
    parts, depth, current = [], 0, ""
    for ch in column:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        if ch == ":" and depth == 0:
            parts.append(current)
            current = ""
        else:
            current += ch
    parts.append(current)
    pieces, names = [], []
    for part in parts:
        m = re.fullmatch(r"(.+)\[T\.(.*)\]", part)
        if m:
            pieces.append((m.group(1), m.group(2)))
            names.append(m.group(1) + m.group(2))
        else:
            pieces.append((part, None))
            names.append(part)
    return ":".join(names), pieces


def eq_coefs(models: list[Model], level: float = 0.95) -> list:
    """The coefficients of every equation, with the node they belong to (``eqCoefs``).

    One DataFrame per equation: ``Target``, ``term`` (as R names it), ``Source`` (the node),
    ``factor``, ``estimate``, ``std.error``, ``statistic``, ``p.value``, ``conf.low``, ``conf.high``
    and ``sd.x``; its ``attrs`` hold ``family``, ``link``, ``sd.y`` and ``latent.sd``. Standard errors
    are those of the model, never robust.
    """
    out = []
    for m in models:
        try:
            out.append(_model_coefs(m, level))
        except Exception as exc:                                           # pragma: no cover
            warnings.warn(f"coefficients not available for equation '{m.name}': {exc}")
            out.append(None)
    return out


def _model_coefs(m: Model, level: float) -> pd.DataFrame:
    w = m.weights
    estimated = m.family in ("gaussian", "Gamma", "inverse.gaussian", "quasibinomial", "quasipoisson")
    se = np.sqrt(np.diag(m.cov_model))
    dfres = (int((w > 0).sum()) - m.X.shape[1]) if estimated else np.inf
    stat = m.coef / se
    crit = stats.t.ppf(1 - (1 - level) / 2, dfres) if np.isfinite(dfres) else stats.norm.ppf(1 - (1 - level) / 2)
    pvalue = 2 * stats.t.sf(np.abs(stat), dfres) if np.isfinite(dfres) else 2 * stats.norm.sf(np.abs(stat))
    names, sources, factor = [], [], []
    for column in m.names:
        if column == "Intercept":
            names.append("(Intercept)")
            sources.append("(Intercept)")
            factor.append(False)
            continue
        name, pieces = _r_column_name(column)
        names.append(name)
        levels = [(v, l) for v, l in pieces if l is not None]
        if levels:
            sources.append("*".join(v for v, _ in pieces) + ":" + "*".join(l for _, l in levels))
        else:
            sources.append(name)
        factor.append(bool(levels))

    def mean(x):
        return (w * x).sum() / w.sum()

    sd_x = [np.nan if c == "Intercept" else float(np.sqrt((w * (m.X[:, j] - mean(m.X[:, j])) ** 2).sum() / w.sum()))
            for j, c in enumerate(m.names)]
    target = m.target
    if m.is_binomial:
        y = m.response
        if target in m.factors:
            target = f"{target}:{m.factors[target][1]}"
        elif pd.api.types.is_bool_dtype(y):
            target = f"{target}:TRUE"
        else:
            u = sorted(pd.Series(y).dropna().unique())
            target = f"{target}:{_text_number(u[1])}" if len(u) > 1 else target
    table = pd.DataFrame({"Target": target, "term": names, "Source": sources, "factor": factor, "estimate": m.coef,
                          "std.error": se, "statistic": stat, "p.value": pvalue, "conf.low": m.coef - crit * se,
                          "conf.high": m.coef + crit * se, "sd.x": sd_x})
    link = m.link
    sd_y = np.nan
    if link == "identity" and pd.api.types.is_numeric_dtype(m.response):
        sd_y = float(np.sqrt((w * (m.y - mean(m.y)) ** 2).sum() / w.sum()))
    latent = np.nan
    error = {"logit": np.pi ** 2 / 3, "probit": 1.0}.get(link)
    if error is not None:
        xb = m.X @ m.coef
        latent = float(np.sqrt((w * (xb - mean(xb)) ** 2).sum() / w.sum() + error))
    table.attrs.update({"family": m.family, "link": link, "sd.y": sd_y, "latent.sd": {target: latent}})
    return table
