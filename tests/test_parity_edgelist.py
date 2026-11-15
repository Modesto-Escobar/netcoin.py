"""``edge_list`` and ``expected_list`` against netCoin 2.1.49 (``edgeList``, ``expectedList``)."""
import numpy as np
import pandas as pd
import pytest

from netcoin import coin, edge_list, expected_list
from conftest import known_gap, read_expected, read_input


def _coin(name):
    return coin(read_input(name), print_n=False)


# case -> (python input, kwargs).  Python procedure keys differ from R's names.
CASES = {
    "default_galapagos": ("galapagos", {}),
    "jaccard_dice_min": ("galapagos", dict(procedures=["Jaccard", "dice"], criteria="Jaccard", min=0.3)),
    "dice_sorted": ("dice", dict(procedures=["Haberman", "Z"], criteria="Haberman", min=2, sort="Haberman")),
    "bonferroni": ("galapagos", dict(procedures="Haberman", bonferroni=True)),
    "support": ("galapagos", dict(procedures=["Haberman", "coincidences", "expected", "confidence"], support=3)),
    "directed": ("galapagos", dict(procedures="Haberman", directed=True, diagonal=True)),
    "conditional": ("galapagos", dict(procedures=["CC", "CP"], criteria="Z")),
    "fisher": ("galapagos", dict(procedures="Hyp", criteria="Hyp")),
    "level90": ("dice", dict(procedures="confidence", criteria="Haberman", min=3, level=0.90)),
}

_HABERMAN = "R Haberman bug (sqrt(N)[is.na(x)]) leaves NaN, python fills sqrt(N)"
GAPS = {
    "default_galapagos": _HABERMAN,
    "support": _HABERMAN,
    "bonferroni": _HABERMAN,
    "directed": _HABERMAN,
}


def _py_input(name, kind="coin"):
    """The python input: the incidence table itself or its ``Coin``."""
    return coin(read_input(name)) if kind == "coin" else read_input(name)


def _key(df, directed):
    pairs = zip(df["Source"], df["Target"])
    return [tuple(p) if directed else tuple(sorted(p)) for p in pairs]


def _align(actual, expected, directed):
    actual = actual.copy()
    actual.index = _key(actual, directed)
    expected = expected.copy()
    expected.index = _key(expected, directed)
    return actual, expected


@pytest.mark.parametrize("kind", ["incidences", "coin"])
@pytest.mark.parametrize("case", CASES)
def test_edge_list_values(case, kind, request):
    known_gap(request, GAPS, case)
    known_gap(request, GAPS, f"{case}:{kind}")
    ds, kwargs = CASES[case]
    expected = read_expected(f"edge_{case}")
    actual = edge_list(_py_input(ds, kind), **kwargs)
    directed = kwargs.get("directed", False)
    assert actual is not None
    actual, expected = _align(actual, expected, directed)

    assert set(actual.columns) == set(expected.columns), "different columns"
    assert set(actual.index) == set(expected.index), "different edges"
    actual = actual.loc[expected.index]
    for col in expected.columns:
        if col in ("Source", "Target"):
            continue
        a, e = actual[col], expected[col]
        if e.dtype == object or str(a.dtype) == "category":
            assert a.astype(str).tolist() == e.astype(str).tolist(), col
        else:
            np.testing.assert_allclose(a.astype(float), e.astype(float),
                                       rtol=1e-8, atol=1e-10, err_msg=col)


@pytest.mark.parametrize("case", ["jaccard_dice_min", "bonferroni", "conditional", "fisher",
                                  "dice_sorted", "level90", "support"])
def test_edge_list_order_orientation_and_columns(case, request):
    """Same rows in the same order, Source/Target as in R, columns as requested."""
    ds, kwargs = CASES[case]
    expected = read_expected(f"edge_{case}")
    actual = edge_list(_py_input(ds), **kwargs)
    assert list(actual.columns) == list(expected.columns)
    assert actual["Source"].tolist() == expected["Source"].tolist()
    assert actual["Target"].tolist() == expected["Target"].tolist()


def test_edge_list_shape_matrix():
    expected = read_expected("edge_shape")
    actual = edge_list(coin(read_input("galapagos")), procedures="shape")
    assert actual["Source"].tolist() == expected["Source"].tolist()
    assert actual["Target"].tolist() == expected["Target"].tolist()
    np.testing.assert_allclose(actual["value"], expected["value"])


@pytest.mark.parametrize("procedures", ["tree", "shape"])
def test_edge_list_adjacency(procedures):
    """R asks for an adjacency list with "shape" on a data.frame; "tree" is the same here."""
    expected = read_expected("edge_adjacency", strings=True)
    actual = edge_list(read_input("adjacency", strings=True), procedures=procedures)
    assert actual[["Source", "Target"]].values.tolist() == expected[["Source", "Target"]].values.tolist()


@pytest.mark.parametrize("case,kwargs", [
    ("expected_list", {}),
    ("expected_list_min2", {"min": 2}),
])
def test_expected_list(case, kwargs):
    expected = read_expected(f"edge_{case}")
    actual = expected_list(coin(read_input("galapagos")), **kwargs)
    assert actual["Source"].tolist() == expected["Source"].tolist()
    assert actual["Target"].tolist() == expected["Target"].tolist()
    for col in ("coincidences", "expected"):
        np.testing.assert_allclose(actual[col], expected[col], rtol=1e-9, err_msg=col)


def test_expected_list_with_confidence():
    """R's version fails (see tests/reference/README.md); check the formula instead."""
    c = coin(read_input("galapagos"))
    out = expected_list(c, confidence=True)
    assert {"coincidences", "expected", "confidence"} <= set(out.columns)
    a, b = c.matrix.loc["Geospiza fortis", "Geospiza magnirostris"], 0
    n = c.n
    f1, f2 = c.matrix.loc["Geospiza magnirostris", "Geospiza magnirostris"], c.matrix.loc["Geospiza fortis", "Geospiza fortis"]
    e = f1 * f2 / n
    expected = max(e + (2 * (e < a) - 1) * 1.64 * np.sqrt(e * (1 - f1 / n) * (1 - f2 / n)), 0)
    row = out[(out.Source == "Geospiza magnirostris") & (out.Target == "Geospiza fortis")]
    assert row["confidence"].iloc[0] == pytest.approx(expected)


def test_edge_list_empty_result_warns():
    with pytest.warns(UserWarning, match="Check max and min"):
        out = edge_list(coin(read_input("galapagos")), "Jaccard", criteria="Jaccard", min=2)
    assert out.empty and list(out.columns[:2]) == ["Source", "Target"]
