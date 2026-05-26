import pandas as pd
import numpy as np

# Coin function

def coin(
    df: pd.DataFrame,
    lower: bool = True,
    diag: bool = True,
    blank_upper: bool = True,
    print_n: bool = True
) -> pd.DataFrame:
    """
    Calcula la matriz de coocurrencias a partir de una matriz dicotómica.
    Compatible con columnas Int64, booleanas, sparse y valores NA.
    """

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df debe ser un pandas DataFrame")

    if df.empty:
        return pd.DataFrame()

    df_bin = df.copy()

    for col in df_bin.columns:
        s = df_bin[col]

        if isinstance(s.dtype, pd.SparseDtype):
            s = s.sparse.to_dense()

        s = pd.to_numeric(s, errors="coerce")
        s = s.fillna(0)
        s = (s > 0).astype("int64")

        df_bin[col] = s

    if print_n:
        print(f"n= {len(df_bin)}")

    cooc = df_bin.T.dot(df_bin)

    if lower:
        k = 0 if diag else -1
        mask = np.tril(np.ones(cooc.shape, dtype=bool), k=k)
        cooc = cooc.where(mask)

    cooc = cooc.astype("Int64")

    if blank_upper:
        return cooc.astype(object).where(cooc.notna(), "")

    return cooc

