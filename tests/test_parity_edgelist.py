"""``edge_list`` and ``expected_list`` against netCoin 2.1.49 (``edgeList``, ``expectedList``)."""
import numpy as np
import pandas as pd
import pytest

from netcoin import coin, edge_list
from netcoin.simEdge import expected_list
from conftest import known_gap, read_expected, read_input


def _coin(name):
    return coin(read_input(name), print_n=False)


# case -> (python input, kwargs).  Python procedure keys differ from R's names.
CASES = {
    "default_galapagos": ("galapagos", {}),
    "jaccard_dice_min": ("galapagos", dict(procedures=["j", "d"], criteria="j", min=0.3)),
    "dice_sorted": ("dice", dict(procedures=["h", "z"], criteria="h", min=2, sort="h")),
    "bonferroni": ("galapagos", dict(procedures="h", bonferroni=True)),
    "support": ("galapagos", dict(procedures=["h", "cc", "e", "con"], support=3)),
    "directed": ("galapagos", dict(procedures="h", directed=True, diagonal=True)),
    "conditional": ("galapagos", dict(procedures=["i", "cp"], criteria="z")),
    "fisher": ("galapagos", dict(procedures="hyp", criteria="hyp")),
    "level90": ("dice", dict(procedures="con", criteria="h", min=3, level=0.90)),
}

# gaps are keyed by "<case>" (both inputs) or "<case>:coin" (only when the python
# input is a coin matrix, which does not carry n)
GAPS = {
    "default_galapagos": "R Haberman bug (sqrt(N)[is.na(x)]) leaves NaN, python fills sqrt(N)",
    "support": "R Haberman bug (sqrt(N)[is.na(x)]) leaves NaN, python fills sqrt(N)",
    "bonferroni": "R Haberman bug (sqrt(N)[is.na(x)]) leaves NaN, python fills sqrt(N)",
    "directed": "R Haberman bug (sqrt(N)[is.na(x)]) leaves NaN, python fills sqrt(N)",
    "dice_sorted:coin": "python coin matrix has no n (dice has n=100): wrong Haberman/Z",
    "level90:coin": "python coin matrix has no n (dice has n=100): wrong Haberman/Z",
}


def _py_input(name, kind="coin"):
    """The python input: the incidence table itself or its ``coin`` matrix."""
    return _coin(name) if kind == "coin" else read_input(name)


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


# R builds the list column by column of the lower triangle (Source is the
# variable that comes first, Target the later one); python goes row by row
# with Source = later variable.
@pytest.mark.parametrize("case", ["jaccard_dice_min", "bonferroni", "conditional", "fisher"])
def test_edge_list_order_and_orientation(case, request):
    request.applymarker(pytest.mark.xfail(
        strict=True, reason="python lists Source/Target swapped and row-major; R column-major"))
    ds, kwargs = CASES[case]
    expected = read_expected(f"edge_{case}")
    actual = edge_list(_py_input(ds), **kwargs)
    assert actual["Source"].tolist() == expected["Source"].tolist()
    assert actual["Target"].tolist() == expected["Target"].tolist()


@pytest.mark.parametrize("case", ["jaccard_dice_min", "support"])
def test_edge_list_column_order(case, request):
    known_gap(request, {"support": "column order differs from the order requested"}, case)
    ds, kwargs = CASES[case]
    expected = read_expected(f"edge_{case}")
    actual = edge_list(_py_input(ds), **kwargs)
    assert list(actual.columns) == list(expected.columns)


def test_edge_list_shape_matrix():
    expected = read_expected("edge_shape")
    actual = edge_list(_coin("galapagos"), procedures="shape")
    actual, expected = _align(actual, expected, False)
    assert set(actual.index) == set(expected.index)
    np.testing.assert_allclose(actual.loc[expected.index, "value"], expected["value"])


def test_edge_list_adjacency():
    """R requests an adjacency list with "shape" on a data.frame; python uses "tree"."""
    expected = read_expected("edge_adjacency", strings=True)
    actual = edge_list(read_input("adjacency", strings=True), procedures="tree")
    assert actual[["Source", "Target"]].values.tolist() == expected[["Source", "Target"]].values.tolist()


def test_edge_list_adjacency_with_shape(request):
    request.applymarker(pytest.mark.xfail(
        strict=True, reason='python needs procedures="tree"; R uses "shape" on a data.frame'))
    expected = read_expected("edge_adjacency", strings=True)
    actual = edge_list(read_input("adjacency", strings=True), procedures="shape")
    assert actual[["Source", "Target"]].values.tolist() == expected[["Source", "Target"]].values.tolist()


@pytest.mark.parametrize("case,kwargs", [
    ("expected_list", {}),
    ("expected_list_min2", {"min": 2}),
])
def test_expected_list(case, kwargs):
    expected = read_expected(f"edge_{case}")
    actual = expected_list(_coin("galapagos"), **kwargs)
    actual, expected = _align(actual, expected, False)
    assert set(actual.index) == set(expected.index)
    for col in ("coincidences", "expected"):
        np.testing.assert_allclose(actual.loc[expected.index, col],
                                   expected[col], rtol=1e-9, err_msg=col)
