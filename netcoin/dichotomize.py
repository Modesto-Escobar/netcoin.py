import pandas as pd
import numpy as np
import re


def _normalize_sep(sep, variables):
    if isinstance(sep, str):
        return {var: sep for var in variables}
    elif isinstance(sep, (list, tuple)):
        if len(sep) == 1:
            return {var: sep[0] for var in variables}
        if len(sep) != len(variables):
            raise ValueError("sep debe tener longitud 1 o la misma longitud que variables")
        return dict(zip(variables, sep))
    elif isinstance(sep, dict):
        return {var: sep.get(var, "") for var in variables}
    else:
        raise TypeError("sep debe ser str, list, tuple o dict")


def _tokenize_series(s: pd.Series, sep: str, regex: bool = False):
    s = s.astype("string")

    if sep == "":
        return s.apply(
            lambda x: [] if pd.isna(x) or str(x).strip() == ""
            else str(x).split()
        )
    else:
        if regex:
            return s.apply(
                lambda x: [] if pd.isna(x) or str(x).strip() == ""
                else [t.strip() for t in re.split(sep, str(x)) if t.strip() != ""]
            )
        else:
            return s.apply(
                lambda x: [] if pd.isna(x) or str(x).strip() == ""
                else [t.strip() for t in str(x).split(sep) if t.strip() != ""]
            )


def dichotomize(
    data: pd.DataFrame,
    variables,
    sep="",
    min_freq=1,
    length=0,
    values=None,
    sparse=False,
    add=False,
    sort=True,
    nas="None",
    stopwords=None,
    regex=False
):
    if not isinstance(data, pd.DataFrame):
        raise TypeError("You must pass a data frame!")

    if isinstance(variables, str):
        variables = [variables]
    else:
        variables = list(variables)

    df = data.copy()
    old_df = df.copy()

    if 0 < min_freq < 1:
        min_freq = min_freq * len(df)
    min_freq = int(np.ceil(min_freq))

    if nas is not None and sep == "C":
        tmp = df[variables].copy()

        def join_row(row):
            vals = []
            for x in row:
                if pd.isna(x):
                    vals.append(np.nan)
                else:
                    vals.append(str(x))
            return "|".join(map(str, vals))

        string_series = tmp.apply(join_row, axis=1)
        string_series = string_series.str.replace(r"\|NA", "", regex=True)
        string_series = string_series.str.replace("NA", str(nas), regex=False)

        df = pd.DataFrame({"String": string_series}, index=df.index)
        variables = ["String"]
        sep = "|"

    sep_map = _normalize_sep(sep, variables)
    result_parts = []
    
    if stopwords is None:
        stopwords = set()
    else:
        stopwords = {str(w).strip() for w in stopwords}

    for c in variables:
        current_sep = sep_map[c]
        tokens = _tokenize_series(df[c], current_sep, regex=regex)
        
        if stopwords:
            tokens = tokens.apply(
                lambda xs: [x for x in xs if str(x).strip() not in stopwords]
        )   

        exploded = tokens.explode()
        exploded = exploded.dropna()
        exploded = exploded[exploded.astype(str).str.strip() != ""]
        exploded = exploded.astype(str).str.strip()

        if exploded.empty:
            W = pd.DataFrame(index=df.index)
        else:
            W = pd.crosstab(exploded.index, exploded).clip(upper=1)
            W = W.reindex(df.index, fill_value=0).astype("Int64")

        if values is not None:
            wanted = [str(v).strip() for v in values]
            W = W.reindex(columns=wanted, fill_value=0).astype("Int64")

        if values is None and min_freq > 0 and W.shape[1] > 0:
            vf = W.sum(axis=0)
            keep = vf[vf >= min_freq].index
            W = W.loc[:, keep]

        if W.shape[1] > 0:
            vf = W.sum(axis=0)

            if sort:
                W = W.loc[:, vf.sort_values(ascending=False).index]
            else:
                W = W.reindex(sorted(W.columns), axis=1)

            if length and length > 0:
                W = W.iloc[:, :min(length, W.shape[1])]

        new_cols = []
        for z in W.columns:
            label = ".ND." if z == "" else str(z)
            if len(variables) > 1:
                label = f"{c}:{label}"
            new_cols.append(label)
        W.columns = new_cols

        if nas is not None:
            na_col = f"{c}:<NA>" if len(variables) > 1 else "<NA>"
            W[na_col] = np.where(W.sum(axis=1) == 0, 1, 0)

        W = W.astype("Int64")
        result_parts.append(W)

    if result_parts:
        R_Data = pd.concat(result_parts, axis=1).astype("Int64")
    else:
        R_Data = pd.DataFrame(index=df.index)

    if nas is not None and not sparse:
        R_Data.columns = [
            col.replace(":<NA>", f":{nas}").replace("<NA>", str(nas)).replace(".ND.", str(nas))
            for col in R_Data.columns
        ]

    if sparse:
        for col in R_Data.columns:
            R_Data[col] = pd.arrays.SparseArray(R_Data[col], fill_value=0)
        return R_Data

    if add:
        return pd.concat([old_df, R_Data], axis=1)

    return R_Data