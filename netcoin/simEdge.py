import math
import numpy as np
import pandas as pd
from scipy.stats import t as student_t, hypergeom


_PROCEDURE_ALIASES = {
    "f": "frequencies",
    "x": "relative",
    "i": "sconditional",
    "cc": "coincidences",
    "cp": "tconditional",
    "e": "expected",
    "con": "confidence",
    "m": "matching",
    "t": "rogers",
    "g": "gower",
    "s": "sneath",
    "and": "anderberg",
    "j": "jaccard",
    "d": "dice",
    "a": "antidice",
    "o": "ochiai",
    "k": "kulczynski",
    "ham": "hamann",
    "y": "yule",
    "p": "pearson",
    "od": "odds",
    "r": "russell",
    "h": "haberman",
    "z": "z",
    "hyp": "fisher",
}

_CONDITIONAL_LABELS = [
    "Null", "Mere", "Conditional", "Significant",
    "Quite significant", "Very significant",
    "Subtotal", "Suptotal", "Total",
]

_PROBABLE_LABELS = [
    "Null", "Mere", "Probable", "Significant",
    "Quite significant", "Very significant",
    "Subtotal", "Suptotal", "Total",
]

_DISPLAY_NAMES = {
    "frequencies": "frequencies",
    "coincidences": "coincidences",
    "relative": "relative",
    "sconditional": "c.conditional",
    "tconditional": "c.probable",
    "expected": "expected",
    "conf_l": "conf.L",
    "confidence": "confidence",
    "conf_u": "conf.U",
    "matching": "matching",
    "rogers": "Rogers",
    "gower": "Gower",
    "sneath": "Sneath",
    "anderberg": "Anderberg",
    "jaccard": "Jaccard",
    "dice": "dice",
    "antidice": "antiDice",
    "ochiai": "Ochiai",
    "kulczynski": "Kulczynski",
    "hamann": "Hamann",
    "yule": "Yule",
    "pearson": "Pearson",
    "odds": "odds",
    "russell": "Russell",
    "haberman": "Haberman",
    "z": "p(Z)",
    "fisher": "p(Fisher)",
}

_PVALUE_CRITERIA = {"z", "hyp"}


def _normalize_procedures(procedures):
    if isinstance(procedures, str):
        procedures = [procedures]

    out = []
    for proc in procedures:
        key = str(proc).strip().lower()
        if key not in _PROCEDURE_ALIASES:
            raise ValueError(
                f"Procedimiento desconocido: '{proc}'.\n"
                f"Opciones válidas: {sorted(_PROCEDURE_ALIASES.keys())}"
            )
        out.append(_PROCEDURE_ALIASES[key])

    return out


def _safe_div(num, den):
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full_like(num, np.nan, dtype=float)
    mask = den != 0
    out[mask] = num[mask] / den[mask]
    return out


def _to_numeric_binary_df(df):
    out = df.copy()

    for col in out.columns:
        s = out[col]
        if isinstance(s.dtype, pd.SparseDtype):
            s = s.sparse.to_dense()

        s = pd.to_numeric(s, errors="coerce").fillna(0)
        out[col] = (s > 0).astype("int64")

    return out


def _is_square_dataframe(df):
    return isinstance(df, pd.DataFrame) and df.shape[0] == df.shape[1]


def _mask_by_diag(df, minimum=1, maximum=np.inf):
    diag_vals = np.diag(df.to_numpy())
    keep = (diag_vals >= minimum) & (diag_vals <= maximum)
    idx = df.index[keep]
    return df.loc[idx, idx]


def _order_by_diag(df, decreasing=False):
    diag_vals = np.diag(df.to_numpy())
    order = np.argsort(diag_vals)

    if decreasing:
        order = order[::-1]

    idx = df.index[order]
    return df.loc[idx, idx]


def _to_df(arr, labels):
    return pd.DataFrame(arr, index=labels, columns=labels)


def _distant(s, distance=False):
    return 1.0 - s if distance else s


def _prepare_counts(
    data,
    input_type="binary",
    minimum=1,
    maximum=np.inf,
    sort=False,
    decreasing=False,
):
    if input_type not in {"binary", "cooccurrence", "auto"}:
        raise ValueError("input_type debe ser 'binary', 'cooccurrence' o 'auto'")

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data debe ser un pandas DataFrame")

    inferred_type = input_type
    if input_type == "auto":
        inferred_type = "cooccurrence" if _is_square_dataframe(data) else "binary"

    if inferred_type == "binary":
        x_bin = _to_numeric_binary_df(data)
        n = float(len(x_bin))
        a_df = x_bin.T.dot(x_bin).astype(float)
    else:
        a_df = data.apply(pd.to_numeric, errors="coerce").astype(float)
        a_raw = a_df.to_numpy(dtype=float)
        i_upper = np.triu_indices_from(a_raw, k=1)
        a_raw[i_upper] = a_raw.T[i_upper]
        a_df = pd.DataFrame(a_raw, index=a_df.index, columns=a_df.columns)
        n = float(np.nanmax(np.diag(a_raw)))

    a_df = _mask_by_diag(a_df, minimum=minimum, maximum=maximum)

    if sort or decreasing:
        a_df = _order_by_diag(a_df, decreasing=decreasing)

    labels = a_df.index
    a = a_df.to_numpy(dtype=float)
    diag_a = np.diag(a)

    b = diag_a[:, None] - a
    c = b.T
    d = n - a - b - c

    return {"a": a, "b": b, "c": c, "d": d, "N": n, "labels": labels}


def sim(
    data,
    procedures="j",
    level=0.95,
    distance=False,
    minimum=1,
    maximum=np.inf,
    sort=False,
    decreasing=False,
    input_type="auto",
):
    if not (0 < float(level) < 1):
        raise ValueError("level debe estar entre 0 y 1")

    procs = _normalize_procedures(procedures)

    obj = _prepare_counts(
        data,
        input_type=input_type,
        minimum=minimum,
        maximum=maximum,
        sort=sort,
        decreasing=decreasing,
    )

    a = obj["a"]
    b = obj["b"]
    c = obj["c"]
    d = obj["d"]
    n = obj["N"]
    labels = obj["labels"]

    m = np.where(a + d == n, 1, np.where(b + c == n, -1, 0)).astype(float)
    results = {}

    if "frequencies" in procs:
        results["frequencies"] = pd.Series(np.diag(a), index=labels)

    if "coincidences" in procs:
        results["coincidences"] = _to_df(a.copy(), labels)

    if "relative" in procs:
        results["relative"] = _to_df(a / n * 100, labels)

    if "sconditional" in procs:
        z_hab = np.sqrt(n) * _safe_div(
            a * d - b * c,
            np.sqrt((a + b) * (a + c) * (b + d) * (d + c)),
        )
        z_sc = 1.0 - student_t.cdf(z_hab, df=max(int(n), 1))
        z_sc = np.where(np.isnan(z_sc), 0.0, z_sc)

        val = np.where(
            b + c == 0, 8,
            np.where(
                c == 0, 7,
                np.where(
                    b == 0, 6,
                    np.where(
                        z_sc < 0.001, 5,
                        np.where(
                            z_sc < 0.01, 4,
                            np.where(
                                z_sc < 0.05, 3,
                                np.where(z_sc < 0.5, 2, np.where(a > 0, 1, 0)),
                            ),
                        ),
                    ),
                ),
            ),
        )
        results["sconditional"] = _to_df(val.astype(float), labels)

    if "tconditional" in procs:
        ac = a + c
        ac_safe = np.where(ac == 0, np.nan, ac)

        t_val = _safe_div(
            a / ac_safe - 0.5,
            1.0 / (2.0 * np.sqrt(ac_safe)),
        )

        df_mat = np.where(ac > 0, ac, 1).astype(float)
        z_tc = 1.0 - student_t.cdf(t_val, df=df_mat)
        z_tc = np.where(np.isnan(z_tc), 0.0, z_tc)

        val = np.where(
            b + c == 0, 8,
            np.where(
                c == 0, 7,
                np.where(
                    z_tc < 0.001, 5,
                    np.where(
                        z_tc < 0.01, 4,
                        np.where(
                            z_tc < 0.05, 3,
                            np.where(z_tc < 0.5, 2, np.where(a > 0, 1, 0)),
                        ),
                    ),
                ),
            ),
        )
        results["tconditional"] = _to_df(val.astype(float), labels)

    if "expected" in procs:
        results["expected"] = _to_df((a + b) * (a + c) / n, labels)

    if "confidence" in procs:
        df_t = max(int(n) - 1, 1)
        exp_ = (a + b) * (a + c) / n
        se = np.sqrt(exp_ * (1 - (a + b) / n) * (1 - (a + c) / n))
        t_level = student_t.ppf(level + (1 - level) / 2, df_t)

        conf_l = np.maximum(a - t_level * se, 0)
        signo = 2 * (exp_ < a) - 1
        confidence = np.maximum(exp_ + signo * student_t.ppf(level, df_t) * se, 0)
        conf_u = np.minimum(a + t_level * se, n)

        np.fill_diagonal(confidence, np.diag(a))

        results["conf_l"] = _to_df(conf_l, labels)
        results["confidence"] = _to_df(confidence, labels)
        results["conf_u"] = _to_df(conf_u, labels)

    if "matching" in procs:
        results["matching"] = _to_df(
            _distant(_safe_div(a + d, a + b + c + d), distance), labels
        )

    if "rogers" in procs:
        results["rogers"] = _to_df(
            _distant(_safe_div(a + d, a + 2 * (b + c) + d), distance), labels
        )

    if "gower" in procs:
        s = _safe_div(a * d, np.sqrt((a + b) * (a + c) * (d + b) * (d + c)))
        s = np.where(np.isnan(s), _distant(np.maximum(m, 0), distance), _distant(s, distance))
        results["gower"] = _to_df(s, labels)

    if "sneath" in procs:
        results["sneath"] = _to_df(
            _distant(_safe_div(2 * (a + d), 2 * (a + d) + b + c), distance), labels
        )

    if "anderberg" in procs:
        s = (
            _safe_div(a, a + b)
            + _safe_div(a, a + c)
            + _safe_div(d, c + d)
            + _safe_div(d, b + d)
        ) / 4
        s = np.where(np.isnan(s), _distant(np.maximum(m, 0), distance), _distant(s, distance))
        results["anderberg"] = _to_df(s, labels)

    if "jaccard" in procs:
        results["jaccard"] = _to_df(
            _distant(_safe_div(a, a + b + c), distance), labels
        )

    if "dice" in procs:
        results["dice"] = _to_df(
            _distant(_safe_div(2 * a, 2 * a + b + c), distance), labels
        )

    if "antidice" in procs:
        results["antidice"] = _to_df(
            _distant(_safe_div(a, a + 2 * (b + c)), distance), labels
        )

    if "ochiai" in procs:
        s = _safe_div(a, np.sqrt((a + b) * (a + c)))
        s = np.where(np.isnan(s), _distant(np.maximum(m, 0), distance), _distant(s, distance))
        results["ochiai"] = _to_df(s, labels)

    if "kulczynski" in procs:
        s = (_safe_div(a, a + b) + _safe_div(a, a + c)) / 2
        s = np.where(np.isnan(s), _distant(np.maximum(m, 0), distance), _distant(s, distance))
        results["kulczynski"] = _to_df(s, labels)

    if "hamann" in procs:
        results["hamann"] = _to_df(
            _distant(_safe_div(a - (b + c) + d, a + b + c + d), distance), labels
        )

    if "yule" in procs:
        s = _safe_div(a * d - b * c, a * d + b * c)
        results["yule"] = _to_df(np.where(np.isnan(s), m, s), labels)

    if "pearson" in procs:
        s = _safe_div(
            a * d - b * c,
            np.sqrt((a + b) * (a + c) * (b + d) * (d + c)),
        )
        s = np.where(np.isnan(s), _distant(m, distance), _distant(s, distance))
        results["pearson"] = _to_df(s, labels)

    if "odds" in procs:
        s = (np.maximum(a, 0.5) * np.maximum(d, 0.5)) / (
            np.maximum(b, 0.5) * np.maximum(c, 0.5)
        )
        if distance:
            s = -s
        np.fill_diagonal(s, -np.inf if distance else np.inf)
        results["odds"] = _to_df(s, labels)

    if "russell" in procs:
        s = _distant(_safe_div(a, a + b + c + d), distance)
        if not distance:
            np.fill_diagonal(s, 1.0)
        results["russell"] = _to_df(s, labels)

    if "haberman" in procs:
        s = np.sqrt(n) * _safe_div(
            a * d - b * c,
            np.sqrt((a + b) * (a + c) * (b + d) * (d + c)),
        )
        s = np.where(np.isnan(s), np.sqrt(n), s)

        if distance:
            s = (n + s) / (2 * n)

        results["haberman"] = _to_df(s, labels)

    if "z" in procs:
        z_val = np.sqrt(n) * _safe_div(
            a * d - b * c,
            np.sqrt((a + b) * (a + c) * (b + d) * (d + c)),
        )
        s = 1.0 - student_t.cdf(z_val, df=max(int(n), 1))
        results["z"] = _to_df(np.where(np.isnan(s), 0.0, s), labels)

    if "fisher" in procs:
        n_min = np.minimum(a + b, a + c)
        n_max = np.maximum(a + b, a + c)
        s = np.full_like(a, np.nan, dtype=float)

        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                m_ = int(n)
                n_ = int(n_min[i, j])
                n_draw = int(n_max[i, j])
                k_ = int(a[i, j]) - 1

                if m_ > 0 and n_ >= 0 and n_draw >= 0:
                    s[i, j] = hypergeom.sf(k_, m_, n_, n_draw)

        results["fisher"] = _to_df(s, labels)

    return next(iter(results.values())) if len(results) == 1 else results


def _mats2edges(
    labels,
    a_matrix: np.ndarray,
    matrices: dict[str, pd.DataFrame],
    criteria_key: str,
    min_val: float,
    max_val: float,
    support: float,
    directed: bool,
    diagonal: bool,
) -> pd.DataFrame:
    n = len(labels)
    rows = []

    for i in range(n):
        for j in range(n):
            if i == j and not diagonal:
                continue

            if i != j and not directed and j > i:
                continue

            if support > -np.inf:
                cc_val = float(a_matrix[i, j])
                if np.isnan(cc_val) or cc_val < support:
                    continue

            mat = matrices.get(criteria_key)
            if mat is not None:
                val = float(mat.iloc[i, j])
                if np.isnan(val) or not (min_val <= val <= max_val):
                    continue

            row = {"Source": labels[i], "Target": labels[j]}

            for key, mat in matrices.items():
                row[_DISPLAY_NAMES.get(key, key)] = mat.iloc[i, j]

            rows.append(row)

    return pd.DataFrame(rows)


def edge_list(
    data: pd.DataFrame | np.ndarray,
    procedures: str | list = "h",
    criteria: str = "z",
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
) -> pd.DataFrame | None:

    if isinstance(procedures, str):
        procedures = [procedures]

    first_proc = procedures[0].lower()
    is_shape = first_proc.startswith("sh")

    # -----------------------------------------------------
    # Caso shape: matriz directa
    # -----------------------------------------------------
    if is_shape:
        if isinstance(data, np.ndarray):
            data = pd.DataFrame(data)

        if not isinstance(data, pd.DataFrame):
            raise TypeError("Para procedures='shape', data debe ser matriz o DataFrame")

        if min == -np.inf:
            min = 1

        mat = data.apply(pd.to_numeric, errors="coerce").astype(float)
        labels = list(mat.index)

        rows = []
        n = mat.shape[0]

        for i in range(n):
            for j in range(n):
                if i == j and not diagonal:
                    continue
                if i != j and not directed and j > i:
                    continue

                value = mat.iloc[i, j]

                if np.isnan(value) or not (min <= value <= max):
                    continue

                rows.append({
                    "Source": labels[i],
                    "Target": labels[j],
                    "value": value,
                })

        out = pd.DataFrame(rows)

        if out.empty:
            print("Warning: Check max and min values")
            return None

        if sort:
            out = out.sort_values("value", ascending=not decreasing).reset_index(drop=True)

        return out

    # -----------------------------------------------------
    # Caso tree: DataFrame de listas de adyacencia
    # -----------------------------------------------------
    if first_proc.startswith("tr"):
        if not isinstance(data, pd.DataFrame):
            raise TypeError("Para procedures='tree', data debe ser un DataFrame")

        lines = data.astype(str).replace("", np.nan).to_numpy()
        rows = []

        for line in lines:
            vals = [x for x in line if pd.notna(x) and str(x).lower() != "nan"]
            if len(vals) < 2:
                continue

            source = vals[0]
            for target in vals[1:]:
                rows.append({"Source": source, "Target": target})

        return pd.DataFrame(rows)

    # -----------------------------------------------------
    # Caso normal: coin / matriz de coocurrencias
    # -----------------------------------------------------
    procedures = [p.lower() for p in procedures]
    criteria = criteria.lower()
    sort_norm = sort.lower() if sort else None

    for p in procedures + [criteria] + ([sort_norm] if sort_norm else []):
        if p not in _PROCEDURE_ALIASES:
            raise ValueError(
                f"Abreviación desconocida: '{p}'.\n"
                f"Opciones válidas: {sorted(_PROCEDURE_ALIASES.keys())}"
            )

    if max == np.inf and criteria in _PVALUE_CRITERIA:
        max = 0.50

    if bonferroni:
        n_pairs = math.comb(data.shape[0], 2)
        if n_pairs > 0:
            max /= n_pairs

    all_procs = list(dict.fromkeys(
        ([sort_norm] if sort_norm else []) + procedures + [criteria]
    ))

    if support > -np.inf and "cc" not in all_procs:
        all_procs.append("cc")

    sim_result = sim(data, all_procs, level=level)

    if not isinstance(sim_result, dict):
        internal = _PROCEDURE_ALIASES[all_procs[0]]
        sim_result = {internal: sim_result}

    if "coincidences" in sim_result:
        a_matrix = sim_result["coincidences"].to_numpy(dtype=float)
    else:
        a_raw = data.apply(pd.to_numeric, errors="coerce").astype(float).to_numpy()
        i_upper = np.triu_indices_from(a_raw, k=1)
        a_raw[i_upper] = a_raw.T[i_upper]
        a_matrix = a_raw

    any_mat = next(v for v in sim_result.values() if isinstance(v, pd.DataFrame))
    labels = list(any_mat.index)

    df_edges = _mats2edges(
        labels=labels,
        a_matrix=a_matrix,
        matrices=sim_result,
        criteria_key=_PROCEDURE_ALIASES[criteria],
        min_val=min,
        max_val=max,
        support=support,
        directed=directed,
        diagonal=diagonal,
    )

    if df_edges.empty:
        print("Warning: Check max and min values")
        return None

    if "c.conditional" in df_edges.columns:
        df_edges["c.conditional"] = pd.Categorical(
            df_edges["c.conditional"].astype(int).map(
                lambda x: _CONDITIONAL_LABELS[x] if 0 <= x <= 8 else str(x)
            ),
            categories=_CONDITIONAL_LABELS,
            ordered=True,
        )

    if "c.probable" in df_edges.columns:
        df_edges["c.probable"] = pd.Categorical(
            df_edges["c.probable"].astype(int).map(
                lambda x: _PROBABLE_LABELS[x] if 0 <= x <= 8 else str(x)
            ),
            categories=_PROBABLE_LABELS,
            ordered=True,
        )

    if sort_norm:
        sort_internal = _PROCEDURE_ALIASES[sort_norm]
        sort_display = _DISPLAY_NAMES.get(sort_internal, sort_internal)

        if sort_display in df_edges.columns:
            df_edges = (
                df_edges
                .sort_values(sort_display, ascending=not decreasing)
                .reset_index(drop=True)
            )

    return df_edges

def expected_list(
    data: pd.DataFrame,
    names: list[str] | None = None,
    min: float = 1,
    confidence: bool = False,
) -> pd.DataFrame | None:
    """
    Devuelve una lista de aristas con:
    - coincidences
    - expected
    - confidence, opcional
    """

    if names is not None:
        data = data.copy()
        data.index = names
        data.columns = names

    a_df = data.apply(pd.to_numeric, errors="coerce").astype(float)
    a = a_df.to_numpy(dtype=float)

    i_upper = np.triu_indices_from(a, k=1)
    a[i_upper] = a.T[i_upper]

    labels = list(a_df.index)
    n = float(np.nanmax(np.diag(a)))

    diag_a = np.diag(a)

    b = diag_a[:, None] - a
    c = b.T
    d = n - a - b - c

    expected = (a + b) * (a + c) / (a + b + c + d)

    expected_df = pd.DataFrame(expected, index=labels, columns=labels)
    coincidences_df = pd.DataFrame(a, index=labels, columns=labels)

    f = edge_list(coincidences_df, procedures="shape", min=0, max=np.inf)
    e = edge_list(expected_df, procedures="shape", min=0, max=np.inf)

    if f is None or e is None:
        return None

    out = f.copy()
    out["expected"] = e["value"]
    out = out.rename(columns={"value": "coincidences"})

    out = out[out["coincidences"] >= min].reset_index(drop=True)

    if not confidence:
        return out

    total = a + b + c + d
    signo = 2 * (((a + b) * (a + c) / total) < a) - 1

    conf = np.maximum(
        ((a + b) * (a + c) / total)
        + signo
        * 1.64
        * np.sqrt(
            ((a + b) * (a + c) / total)
            * ((1 - (a + b) / total) * (1 - (a + c) / total))
        ),
        0,
    )

    np.fill_diagonal(conf, np.diag(a))

    conf_df = pd.DataFrame(conf, index=labels, columns=labels)
    l = edge_list(conf_df, procedures="shape", min=-np.inf, max=np.inf)

    if l is not None:
        out["confidence"] = l["value"]

    return out


def distant(s, t: bool = False):
    """
    Si t=True, convierte similitud en distancia: 1 - s.
    """
    return 1.0 - s if t else s


def lower(matrix, decimals: int = 3) -> pd.DataFrame:
    """
    Devuelve solo el triángulo inferior formateado como texto.
    """

    if isinstance(matrix, pd.DataFrame):
        labels = list(matrix.index)
        arr = matrix.to_numpy(dtype=float)
    else:
        arr = np.asarray(matrix, dtype=float)
        labels = list(range(arr.shape[0]))

    formatted = np.empty(arr.shape, dtype=object)

    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if j > i:
                formatted[i, j] = ""
            else:
                formatted[i, j] = f"{arr[i, j]:.{decimals}f}"

    return pd.DataFrame(formatted, index=labels, columns=labels)