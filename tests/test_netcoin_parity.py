"""The ``netCoin`` object, layouts, igraph conversion and exports against netCoin 2.1.49."""
import json

import numpy as np
import pandas as pd
import pytest

from netcoin import (cal_centr, layout_circle, layout_grid, netCoin, save_pajek, to_igraph)
from conftest import EXPECTED, read_expected, read_input

igraph = pytest.importorskip("igraph")

NODES = read_input("net_nodes", strings=True).astype({"size": float, "value": float})
LINKS = read_input("net_links", strings=True).astype({"weight": float})
TREE = read_input("net_tree", strings=True)
LAYOUT = np.column_stack([np.arange(10), np.tile([0, 1], 5)]).astype(float)
LAYOUT_DF = pd.DataFrame(LAYOUT, columns=["x", "y"])


def _ref(case):
    return json.loads((EXPECTED / f"net_{case}.json").read_text(encoding="utf-8"))


def _as_list(value):
    return value if isinstance(value, list) else [value]


def _compare_options(actual: dict, expected: dict, case: str):
    assert set(actual) == set(expected), f"{case}: {sorted(set(actual) ^ set(expected))}"
    for key, value in expected.items():
        got = actual[key]
        if isinstance(got, dict) and isinstance(value, list):
            got = list(got.values())          # R's plain() dropped the names of the planes
        if isinstance(value, dict):
            assert set(got) == set(value), f"{case}/{key}"
            for k, v in value.items():
                assert _as_list(got[k]) == _as_list(v), f"{case}/{key}/{k}"
            continue
        got, value = _as_list(got), _as_list(value)
        assert len(got) == len(value), f"{case}/{key}: {got} != {value}"
        for g, v in zip(got, value):
            if isinstance(v, float):
                assert g == pytest.approx(v), f"{case}/{key}"
            else:
                assert g == v, f"{case}/{key}: {g!r} != {v!r}"


def _compare_frame(actual: pd.DataFrame, expected: pd.DataFrame, label: str, rtol=1e-8, up_to_sign=()):
    assert list(actual.columns) == list(expected.columns), f"{label}: columns {list(actual.columns)}"
    assert len(actual) == len(expected), f"{label}: {len(actual)} rows, R has {len(expected)}"
    for col in expected.columns:
        a, e = actual[col].reset_index(drop=True), expected[col].reset_index(drop=True)
        if pd.api.types.is_numeric_dtype(e) and pd.api.types.is_numeric_dtype(a):
            a, e = a.astype(float), e.astype(float)
            if col in up_to_sign and np.nansum(np.abs(a + e)) < np.nansum(np.abs(a - e)):
                a = -a
            np.testing.assert_allclose(a, e, rtol=rtol, atol=1e-7, equal_nan=True, err_msg=f"{label}/{col}")
        else:
            norm = lambda s: s.map(lambda x: "NA" if pd.isna(x) else str(x))
            assert norm(a).tolist() == norm(e).tolist(), f"{label}/{col}"


CASES = {
    "basic": lambda: netCoin(NODES, LINKS, label="label", size="size", color="group", shape="kind",
                             group="group", lwidth="weight", lcolor="type", ltext="note", ntext="text",
                             info="text", main="A network", note="A note", scenarios=10),
    "single": lambda: netCoin(NODES, LINKS, color="firebrick", shape="Square", lcolor="#00f", size=12,
                              label="", legend="kind", border="group", label_size="size"),
    "numeric": lambda: netCoin(NODES, LINKS, color="value", node_scale_limits=[0, 5, 14], lcolor="weight",
                               link_scale_limits=[0, 1], intensity="weight", node_bipolar=True,
                               link_bipolar=True),
    "lists": lambda: netCoin(NODES, LINKS, color={"group": {"g1": "red", "g2": "blue", "g3": "green"}},
                             shape={"kind": {"x": "Circle", "y": "star"}},
                             lcolor={"type": {"p": "orange", "q": "grey", "r": "black"}}),
    "vectors": lambda: netCoin(NODES, LINKS, color=["red", "#0000FF", "light blue", "gray 50", "Forest Green"] * 2,
                               shape=list(range(1, 11)), lcolor=list(range(1, 15))),
    "categories": lambda: netCoin(NODES, LINKS, color=pd.Categorical(NODES["group"]),
                                  shape=pd.Categorical(NODES["kind"]), lcolor=pd.Categorical(LINKS["type"])),
    "filters": lambda: netCoin(NODES, LINKS, node_filter="size > 2", link_filter="weight > 0.3"),
    "degree": lambda: netCoin(NODES, LINKS, degree_filter=[2, 4]),
    "degree_one": lambda: netCoin(NODES, LINKS, degree_filter=3, node_filter="kind == 'x'"),
    "options": lambda: netCoin(NODES, LINKS, controls=[1, 2, 3], show_legend=False, mode="network",
                               repulsion=40, distance=15, zoom=2, limits=[0, 0, 10, 10], background="#ABC",
                               default_color="firebrick", help_on=True, help="Help text",
                               show_coordinates=True, show_arrows=True, show_axes=True,
                               axes_labels=["X", "Y"], nodes_label="N", links_label="L", sort="size",
                               decreasing=True, frequencies=True, statistics=True, rounded_items=True,
                               cex=1.3, group_text=True, group="group"),
    "controls": lambda: netCoin(NODES, LINKS, controls=[1, 5, 8], show_legend=False),
    "controls0": lambda: netCoin(NODES, LINKS, controls=0),
    "spanish": lambda: netCoin(NODES, LINKS, language="es", name="name"),
    "tree": lambda: netCoin(NODES, LINKS, tree=TREE),
    "linksonly": lambda: netCoin(links=LINKS),
    "nolinks": lambda: netCoin(NODES, size="size", name="label"),
    "images": lambda: netCoin(NODES, LINKS, image=["label", "text"], image_names=["group", "kind"]),
    "sparse": lambda: netCoin(NODES.iloc[[0, 1, 3]][["name", "size"]], LINKS),
    "layout_matrix": lambda: netCoin(NODES, LINKS, layout=LAYOUT_DF),
    "layout_nocols": lambda: netCoin(NODES, LINKS, layout=LAYOUT),
    "layout_planes": lambda: netCoin(NODES, LINKS, layout={
        "first": LAYOUT_DF, "second": pd.DataFrame({"u": LAYOUT[:, 1], "v": LAYOUT[:, 0]})}),
    "layout_circle": lambda: netCoin(NODES, LINKS, layout="circle"),
    "layout_grid": lambda: netCoin(NODES, LINKS, layout="grid"),
    "layout_star": lambda: netCoin(NODES, LINKS, layout="star"),
    "layout_reingold_tilford": lambda: netCoin(NODES, LINKS, layout="reingold_tilford"),
    "layout_kamada_kawai": lambda: netCoin(NODES, LINKS, layout="kamada_kawai"),
    "layout_sugiyama": lambda: netCoin(NODES, LINKS, layout="sugiyama"),
    "layout_two_names": lambda: netCoin(NODES, LINKS, layout=["circle", "grid"]),
    "layout_mixed": lambda: netCoin(NODES, LINKS, layout={"mine": LAYOUT_DF, "circle": "circle"}),
    "community_walktrap": lambda: netCoin(NODES, LINKS, community="walktrap"),
    "community_fast_greedy": lambda: netCoin(NODES, LINKS, community="fast_greedy"),
    "community_leading_eigen": lambda: netCoin(NODES, LINKS, community="leading_eigen"),
    "community_edge_betweenness": lambda: netCoin(NODES, LINKS, community="edge_betweenness"),
    "community_color": lambda: netCoin(NODES, LINKS, community="walktrap", color="group"),
    "centrality": lambda: netCoin(NODES, LINKS, centrality=["wdegree", "closeness", "betweenness", "eigen"]),
}

# layouts whose coordinates depend on the algorithm's start, the library or a sign
UP_TO_SIGN = {"layout_kamada_kawai": ["fx", "fy"], "layout_sugiyama": []}
APPROX_LAYOUT = {"layout_kamada_kawai": 1e-4, "layout_star": 1e-6, "layout_reingold_tilford": 1e-6}
_CORE = "the igraph C core differs between R (igraph 2.3.4) and python-igraph 1.0.0"
GAPS = {
    "layout_sugiyama": f"Sugiyama layering differs: {_CORE}",
    "layout_kamada_kawai": f"Kamada-Kawai coordinates differ: {_CORE}",
}


def _same_partition(actual, expected, label):
    pairs = set(zip(actual, expected))
    assert len({a for a, _ in pairs}) == len(pairs) == len({e for _, e in pairs}),         f"{label}: different groups ({sorted(pairs)})"


@pytest.mark.parametrize("case", CASES)
def test_netcoin_matches_r(case, request):
    if case in GAPS:
        request.applymarker(pytest.mark.xfail(strict=True, reason=GAPS[case]))
    ref = _ref(case)
    if case in ("images",):
        net = CASES[case]()
    else:
        net = CASES[case]()
    expected_nodes = read_expected(f"net_{case}_nodes")
    rtol = APPROX_LAYOUT.get(case, 1e-8)
    if case.startswith("community_"):
        # the numbers of the groups come from the algorithm's internals and may be
        # permuted; the partition must be the same
        _same_partition(net.nodes["community"], expected_nodes["community"], case)
        expected_nodes = expected_nodes.assign(community=net.nodes["community"].to_numpy())
    _compare_frame(net.nodes, expected_nodes, f"{case}/nodes", rtol=rtol, up_to_sign=UP_TO_SIGN.get(case, ()))
    assert net.lcount() == ref["nlinks"], f"{case}: links"
    if ref["nlinks"]:
        _compare_frame(net.links, read_expected(f"net_{case}_links"), f"{case}/links")
    options = ref["options"]
    if case == "community_color":
        options = dict(options)
    _compare_options(net.options, options, case)
    if "tree" in ref:
        tree = {k: _as_list(v) for k, v in ref["tree"].items()}
        assert net.tree["Source"].tolist() == tree["Source"] and net.tree["Target"].tolist() == tree["Target"]
    else:
        assert net.tree is None
    if "layouts" in ref:
        assert list(net.layouts) == list(ref["layouts"])
        for key, rows in ref["layouts"].items():
            np.testing.assert_allclose(net.layouts[key].to_numpy(), np.array(rows, dtype=float),
                                       rtol=1e-8, atol=1e-7, err_msg=f"{case}/{key}")
    else:
        assert net.layouts is None


# ------------------------------------------------------------ igraph ----

def test_to_igraph_matches_r():
    net = netCoin(NODES, LINKS, label="label", size="size", color="group", lwidth="weight", lcolor="type",
                  show_arrows=True)
    graph = to_igraph(net)
    meta = json.loads((EXPECTED / "igraph_meta.json").read_text(encoding="utf-8"))
    assert graph.is_directed() == meta["directed"]
    vertices = graph.get_vertex_dataframe().reset_index(drop=True)
    expected_v = read_expected("igraph_vertices")
    for col in expected_v.columns:
        assert col in vertices.columns, col
        got, want = vertices[col].reset_index(drop=True), expected_v[col]
        if pd.api.types.is_numeric_dtype(want):
            np.testing.assert_allclose(got.astype(float), want.astype(float), err_msg=col)
        else:
            assert got.astype(str).tolist() == want.astype(str).tolist(), col
    edges = graph.get_edge_dataframe().reset_index(drop=True)
    expected_e = read_expected("igraph_edges")
    names = graph.vs["name"]
    assert [names[i] for i in edges["source"]] == expected_e["from"].tolist()
    assert [names[i] for i in edges["target"]] == expected_e["to"].tolist()
    for col in expected_e.columns[2:]:
        got, want = edges[col].reset_index(drop=True), expected_e[col]
        if pd.api.types.is_numeric_dtype(want):
            np.testing.assert_allclose(got.astype(float), want.astype(float), err_msg=col)
        else:
            assert got.astype(str).tolist() == want.astype(str).tolist(), col


def test_from_igraph_roundtrip():
    from netcoin import from_igraph
    net = netCoin(NODES, LINKS, size="size", lwidth="weight", show_arrows=True, main="Roundtrip")
    back = from_igraph(to_igraph(net))
    assert back.options["showArrows"] is True and back.options["nodeSize"] == "size"
    assert back.options["linkWidth"] == "width" or back.options["linkWidth"] == "weight"
    assert back.nodes["name"].tolist() == NODES["name"].tolist()
    assert back.links[["Source", "Target"]].values.tolist() == LINKS[["Source", "Target"]].values.tolist()


def test_from_igraph_takes_graph_name_and_attributes():
    from netcoin import from_igraph
    g = igraph.Graph([(0, 1), (1, 2)], directed=False)
    g.vs["name"] = ["x", "y", "z"]
    g.vs["size"] = [1, 2, 3]
    g.es["weight"] = [0.5, 0.7]
    g["name"] = "Small graph"
    net = from_igraph(g)
    assert net.options["main"] == "Small graph" and net.options["nodeSize"] == "size"
    assert net.options["linkWeight"] == "weight" and "showArrows" not in net.options
    assert net.links["weight"].tolist() == [0.5, 0.7]


# ---------------------------------------------------------- centrality ----

def test_cal_centr_matches_r():
    net = netCoin(NODES, LINKS)
    got = cal_centr(net, "all")
    expected = read_expected("centrality_all")
    _compare_frame(got, expected, "centrality_all")


def test_cal_centr_directed_matches_r():
    net = netCoin(NODES, LINKS, show_arrows=True)
    got = cal_centr(net, ["degree", "closeness"])
    _compare_frame(got, read_expected("centrality_directed"), "centrality_directed")


def test_cal_centr_orders_the_result():
    out = cal_centr(netCoin(NODES, LINKS), ["degree", "betweenness"], order="betweenness")
    assert out["betweenness"].is_monotonic_decreasing


# ---------------------------------------------------- circle and grid ----

def _same_layout(actual: pd.DataFrame, name: str):
    expected = read_expected(name).set_index(".row")
    assert list(actual.index) == [str(i) for i in expected.index]
    np.testing.assert_allclose(actual.to_numpy(), expected.to_numpy(dtype=float), atol=1e-9, equal_nan=True)


def test_layout_circle_matches_r():
    _same_layout(layout_circle(NODES, name="name"), "layoutcircle_all")
    _same_layout(layout_circle(NODES, name="name", deg=90), "layoutcircle_deg")
    _same_layout(layout_circle(NODES, items=["c", "a", "e", "zz"], name="name"), "layoutcircle_subset")
    _same_layout(layout_circle(NODES, items=["b"], name="name"), "layoutcircle_one")


@pytest.mark.parametrize("case", ["simple", "gaps", "bottom", "one", "byrow"])
def test_layout_grid_matches_r(case):
    strings = json.loads((EXPECTED / "layoutgrid_strings.json").read_text(encoding="utf-8"))
    _same_layout(layout_grid(NODES, strings[case], name="name", byrow=case == "byrow"), f"layoutgrid_{case}")


def test_layout_circle_by_position_and_default_names():
    out = layout_circle(NODES.iloc[:4])
    assert list(out.index) == ["1", "2", "3", "4"]
    assert layout_circle(NODES.iloc[:4], items=[0, 2]).iloc[:, 0].tolist()[:1] == [1.0]


# -------------------------------------------------------------- exports ----

def _lines(text: str):
    return [line.rstrip("\r") for line in text.splitlines()]


def test_save_pajek_matches_r(tmp_path):
    net = netCoin(NODES, LINKS, main="Pajek test", show_arrows=True)
    path = save_pajek(net, str(tmp_path / "pajek_basic.net"), arcs=["weight"], edges=["weight"])
    got = _lines(open(path, encoding="utf-8").read())
    want = _lines((EXPECTED / "pajek_basic.txt").read_text(encoding="utf-8"))
    assert got == want


def test_save_pajek_project_file(tmp_path):
    net = netCoin(NODES, LINKS, main="Pajek test", show_arrows=True)
    path = save_pajek(net, str(tmp_path / "pajek_project.net"), arcs=["weight"],
                      partitions=["group", "kind"], vectors=["size"])
    assert path.endswith(".paj")
    got = _lines(open(path, encoding="utf-8").read())
    want = _lines((EXPECTED / "pajek_project.txt").read_text(encoding="utf-8"))
    assert got == want


def test_save_pajek_adds_the_extension_and_checks_the_columns(tmp_path):
    net = netCoin(NODES, LINKS)
    assert save_pajek(net, str(tmp_path / "plain")).endswith("plain.net")
    with pytest.raises(ValueError, match="partition"):
        save_pajek(net, str(tmp_path / "x"), partitions=["nope"])
    with pytest.raises(ValueError, match="arc"):
        save_pajek(net, str(tmp_path / "x"), arcs=["nope"])


def test_save_graphml_roundtrip(tmp_path):
    from netcoin import save_graphml
    net = netCoin(NODES, LINKS, size="size", lwidth="weight", color=pd.Categorical(NODES["group"]))
    path = save_graphml(net, str(tmp_path / "net"))
    assert path.endswith(".graphml")
    graph = igraph.Graph.Read_GraphML(path)
    assert graph.vcount() == 10 and graph.ecount() == 14
    assert graph.vs["name"][0] == "a" and graph.es["weight"][0] == pytest.approx(0.9)
    with pytest.raises(TypeError):
        save_graphml("not a network", str(tmp_path / "x"))
