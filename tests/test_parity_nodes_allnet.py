"""``coocur``, ``prop_coin``, ``as_nodes`` and ``all_net`` against netCoin 2.1.49."""
import json

import numpy as np
import pandas as pd
import pytest

from netcoin import coin
from netcoin.allnet import all_net
from netcoin.coincidence import as_nodes, coocur, prop_coin
from conftest import (EXPECTED, assert_matrix_close, read_expected, read_input, read_matrix)

# ---------------------------------------------------------------- coocur ----


@pytest.mark.parametrize("case,kwargs", [
    ("default", {}),
    ("min3_sort", {"minimum": 3, "sort": True}),
])
def test_coocur(case, kwargs):
    result = coocur(read_input("counts"), **kwargs)
    assert_matrix_close(result.matrix, read_matrix(f"coocur_{case}"), label=case)
    attrs = read_expected(f"coocur_{case}_attrs")
    assert result.n == pytest.approx(attrs["n"][0])
    assert result.m == pytest.approx(attrs["m"][0])


# -------------------------------------------------------------- prop_coin ----

@pytest.mark.parametrize("margin", [0, 1, 2])
def test_prop_coin_margins(margin):
    got = prop_coin(coin(read_input("galapagos")), margin=margin)
    assert_matrix_close(got, read_matrix(f"prop_margin{margin}"), atol=1e-9)


def test_prop_coin_decimals():
    got = prop_coin(coin(read_input("galapagos")), decimals=3)
    assert_matrix_close(got, read_matrix("prop_margin0_dec3"), atol=1e-9)


def test_prop_coin_pairwise():
    got = prop_coin(coin(read_input("galapagos_na"), pairwise=True))
    assert_matrix_close(got, read_matrix("prop_pairwise"), atol=1e-9)


def test_prop_coin_bad_margin():
    with pytest.raises(ValueError):
        prop_coin(coin(read_input("galapagos")), margin=3)


# --------------------------------------------------------------- as_nodes ----

@pytest.mark.parametrize("case,args", [
    ("default", dict()),
    ("both", dict(frequency=True, percentages=True)),
    ("pct_es", dict(frequency=False, percentages=True, language="es")),
    ("freq_ca", dict(frequency=True, percentages=False, language="ca")),
])
def test_as_nodes(case, args):
    got = as_nodes(coin(read_input("galapagos")), **args)
    expected = read_expected(f"nodes_{case}")
    assert list(got.columns) == list(expected.columns)
    pd.testing.assert_frame_equal(got, expected, check_dtype=False, rtol=1e-9)


def test_as_nodes_pairwise():
    got = as_nodes(coin(read_input("galapagos_na"), pairwise=True), True, True)
    pd.testing.assert_frame_equal(got, read_expected("nodes_pairwise"), check_dtype=False, rtol=1e-9)


def test_as_nodes_from_edges():
    from netcoin import edge_list
    edges = edge_list(coin(read_input("galapagos")))
    got = as_nodes(edges)
    pd.testing.assert_frame_equal(got, read_expected("nodes_edges"), check_dtype=False)


# ---------------------------------------------------------------- all_net ----

# python arguments of each R call
ALLNET = {
    "default": ("galapagos", {}),
    "links": ("galapagos", dict(procedures=["Jaccard", "Haberman"], lwidth="Jaccard", lcolor="Z",
                                criteria="Jaccard", min_link=0.3, language="es", frequency=True,
                                percentages=False, layout="pc")),
    "nodes": ("galapagos", dict(nodes="user_nodes")),
    "weight": ("galapagos", dict(weight="weights", support=3, sort=True)),
    "mca": ("dice", dict(layout="mc", procedures="Jaccard", criteria="Jaccard", min_link=0.1)),
    "bonferroni": ("dice", dict(bonferroni=True, max_link=0.2, directed=True,
                                sort_links="Haberman")),
}
# R leaves NaN in Haberman where python gives sqrt(N) (see reference/README.md)
R_NAN_COLUMNS = {"Haberman"}
LAYOUT_COLUMNS = ["fx", "fy"]


def _run(case):
    source, kwargs = ALLNET[case]
    kwargs = dict(kwargs)
    if kwargs.get("nodes") == "user_nodes":
        kwargs["nodes"] = read_input("user_nodes", strings=True)
        kwargs["nodes"]["frequency"] = 999.0
    if kwargs.get("weight") == "weights":
        kwargs["weight"] = read_input("weights")["weight"].tolist()
    return all_net(read_input(source), **kwargs)


def _compare(actual: pd.DataFrame, expected: pd.DataFrame, label: str, up_to_sign=()):
    assert list(actual.columns) == list(expected.columns), label
    assert len(actual) == len(expected), f"{label}: {len(actual)} rows, R has {len(expected)}"
    for col in expected.columns:
        a, e = actual[col].reset_index(drop=True), expected[col].reset_index(drop=True)
        if pd.api.types.is_numeric_dtype(e):
            a = a.astype(float)
            e = e.astype(float)
            if col in up_to_sign:
                if np.nansum(np.abs(a + e)) < np.nansum(np.abs(a - e)):
                    a = -a
            known = e.notna() if col in R_NAN_COLUMNS else pd.Series(True, index=e.index)
            np.testing.assert_allclose(a[known], e[known], rtol=1e-7, atol=1e-9,
                                       equal_nan=True, err_msg=f"{label}/{col}")
        else:
            assert a.astype(str).tolist() == e.astype(str).tolist(), f"{label}/{col}"


@pytest.mark.parametrize("case", ALLNET)
def test_all_net_nodes(case):
    net = _run(case)
    _compare(net.nodes_frame(), read_expected(f"allnet_{case}_nodes"), f"{case}/nodes",
             up_to_sign=LAYOUT_COLUMNS)


@pytest.mark.parametrize("case", ALLNET)
def test_all_net_links(case):
    net = _run(case)
    _compare(net.links_frame(), read_expected(f"allnet_{case}_links"), f"{case}/links")


@pytest.mark.parametrize("case", ALLNET)
def test_all_net_options(case):
    net = _run(case)
    ref = json.loads((EXPECTED / f"allnet_{case}_options.json").read_text(encoding="utf-8"))
    options = net.options
    assert options["nodeName"] == ref["nodeName"]
    assert options["scenarios"] == pytest.approx(ref["scenarios"])
    assert options.get("nodeSize") == ref.get("nodeSize")
    assert options["language"] == ref["language"]
    if "linkWidth" in ref:
        assert options["linkWidth"] == ref["linkWidth"]


def test_all_net_link_color_points_to_the_p_value_column():
    """R sets lcolor="Z" although the column is "p(Z)"; python points to the real column."""
    net = _run("links")
    assert net.options["linkColor"] == "p(Z)"
    assert "p(Z)" in net.links


def test_all_net_pairwise_with_missing_values_works():
    """In R this call fails (colSums is NA); python ignores the missing values."""
    net = all_net(read_input("galapagos_na"), pairwise=True)
    assert net.ncount() == 13


def test_all_net_igraph_not_available_yet():
    with pytest.raises(NotImplementedError):
        all_net(read_input("galapagos"), igraph=True)


# --------------------------------------------------------------- netCoin.add ----

def test_add_changes_only_what_is_given():
    net = _run("default")
    other = net.add(main="A title", color="%")
    assert other.options["main"] == "A title" and other.options["nodeColor"] == "%"
    assert "main" not in net.options and "nodeColor" not in net.options
    assert other.options["nodeSize"] == net.options["nodeSize"]
    assert other.nodes_frame().equals(net.nodes_frame())


def test_add_can_replace_nodes_and_links():
    net = _run("default")
    nodes = net.nodes_frame().assign(group=lambda d: np.arange(len(d)) % 2)
    other = net.add(nodes=nodes, group="group")
    assert "group" in other.nodes and other.options["nodeGroup"] == "group"


def test_add_rejects_unknown_arguments():
    with pytest.raises(TypeError):
        _run("default").add(colour="x")


def test_all_net_writes_valid_json(tmp_path):
    import re
    out = tmp_path / "net"
    all_net(read_input("galapagos"), layout="pc", procedures=["Jaccard", "Haberman"],
            lwidth="Jaccard", directory=str(out))
    html = (out / "index.html").read_text(encoding="utf-8")
    payload = re.search(r'<script type="application/json" id="data">(.*?)</script>', html, re.S)
    data = json.loads(payload.group(1), parse_constant=lambda c: pytest.fail(f"{c} in JSON"))
    assert data["nodenames"][0] == "name" and len(data["nodes"][0]) == 13
    assert data["options"]["linkWidth"] == "Jaccard"
