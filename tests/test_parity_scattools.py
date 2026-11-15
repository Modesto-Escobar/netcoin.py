"""``add_axes``, ``add_clusters``, ``replace_clusters`` and ``layouts_align`` against netCoin 2.1.49."""
import json
import warnings

import numpy as np
import pandas as pd
import pytest

from netcoin import add_axes, add_clusters, layouts_align, netCoin, replace_clusters
from conftest import EXPECTED, read_expected, read_input

BASE = read_input("sa_base", strings=True).astype({"v": float, "x": float, "y": float})
AXES = read_input("sa_axes")
CASES_DATA = read_input("sa_cases")
CL = read_input("sa_clusters")
NODES = BASE[["name", "v", "g"]]
XY = BASE[["x", "y"]].to_numpy()
CASE_AXES = CASES_DATA[["C1", "C2"]]
IDX, CW = CASES_DATA["idx"].to_numpy(), CASES_DATA["w"].to_numpy()
CL_CASES = CASES_DATA["cl"].to_numpy()


def base(patterns=False):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        net = netCoin(NODES, layout=XY, size="v", axes_labels=["PC1 (60.0%)", "PC2 (25.0%)"])
    net.scat = {"plane_name": "PCA"}
    net.clusters, net.cluster_means = {}, None
    if patterns:
        net.scat["case_to_pattern"], net.scat["case_weight"] = IDX, CW
    return net


def three():
    return add_axes(base(), AXES, name="p")


CASES = {
    "ax_matrix": lambda: add_axes(base(), AXES),
    "ax_named": lambda: add_axes(base(), AXES, name="pair"),
    "ax_single": lambda: add_axes(base(), AXES.iloc[:, :2], name="one"),
    "ax_twice": lambda: add_axes(add_axes(base(), AXES.iloc[:, :2], name="one"), AXES.iloc[:, 2:], name="one"),
    "ax_cases": lambda: add_axes(base(True), CASE_AXES),
    "ax_cases_w": lambda: add_axes(base(True), CASE_AXES, weight=np.ones(100)),
    "ax_labels": lambda: add_axes(base(), AXES, axes_labels=["X", "Y"]),
    "cl_vec": lambda: add_clusters(base(), CL["a"].to_numpy(), name="grp"),
    "cl_nosort": lambda: add_clusters(base(), CL["a"].to_numpy(), name="grp", sort=False),
    "cl_weight": lambda: add_clusters(base(), CL["a"].to_numpy(), name="grp", weight=NODES["v"].to_numpy()),
    "cl_df": lambda: add_clusters(base(), CL),
    "cl_cases": lambda: add_clusters(base(True), CL_CASES, name="g"),
    "cl_rest": lambda: add_clusters(base(), CL["a"].to_numpy(), name="g", rest=3),
    "cl_twice": lambda: add_clusters(add_clusters(base(), CL["a"].to_numpy(), name="g"), CL["b"].to_numpy(), name="g"),
    "cl_replace": lambda: replace_clusters(add_clusters(base(), CL["a"].to_numpy(), name="g"), CL["b"].to_numpy()),
    "cl_maxgroups": lambda: add_clusters(base(), np.arange(1, 41)),
    "cl_axes_then": lambda: add_clusters(add_axes(base(), AXES, name="p"), CL["a"].to_numpy(), name="g"),
    "cl_after_axes_cases": lambda: add_axes(add_clusters(base(True), CL_CASES, name="g"), CASE_AXES),
    "al_proc": lambda: layouts_align(add_axes(base(), AXES.iloc[:, 2:], name="B")),
    "al_proc3": lambda: layouts_align(three()),
    "al_to2": lambda: layouts_align(three(), to=2),
    "al_suffix": lambda: layouts_align(three(), suffix="*"),
    "al_gpa": lambda: layouts_align(three(), gpa=True),
    "al_gpa_suffix": lambda: layouts_align(three(), gpa=True, suffix="~"),
    "al_cloud": lambda: layouts_align(three(), align="cloud"),
    "al_canvas": lambda: layouts_align(three(), align="canvas"),
    "al_none": lambda: layouts_align(three(), align="none"),
    "al_weights": lambda: layouts_align(three(), weights="v"),
    "al_which": lambda: layouts_align(three(), which=[2, 3]),
    "al_clusters": lambda: layouts_align(add_clusters(three(), CL["a"].to_numpy(), name="g"), suffix="*"),
    "al_cases": lambda: layouts_align(add_axes(base(True), CASE_AXES)),
}


def _info(case):
    return json.loads((EXPECTED / f"sa_{case}_info.json").read_text(encoding="utf-8"))


def _run(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return CASES[case]()


def _matrix(entry):
    return np.array(entry["data"], dtype=float).reshape(-1, len(entry["data"][0]) if entry["data"] and isinstance(entry["data"][0], list) else 1)


@pytest.mark.parametrize("case", CASES)
def test_nodes_match_r(case):
    net = _run(case)
    want = read_expected(f"net_sa_{case}_nodes", strings=True)
    assert list(net.nodes.columns) == list(want.columns), f"{case}: {list(net.nodes.columns)}"
    for col in want.columns:
        got = net.nodes[col]
        if col in ("fx", "fy", "v") or pd.api.types.is_float_dtype(got):
            np.testing.assert_allclose(got.astype(float), want[col].astype(float), rtol=1e-8, atol=1e-9, err_msg=f"{case}/{col}")
        else:
            assert got.astype(str).tolist() == want[col].astype(str).tolist(), f"{case}/{col}"


@pytest.mark.parametrize("case", CASES)
def test_planes_match_r(case):
    net, info = _run(case), _info(case)
    layouts = info["layouts"]
    if not layouts or layouts == {}:
        assert not net.layouts
        return
    assert list(net.layouts) == list(layouts), case
    for name, plane in layouts.items():
        want = np.array(plane["data"], dtype=float)
        np.testing.assert_allclose(net.layouts[name].to_numpy(), want, rtol=1e-6, atol=1e-8, err_msg=f"{case}/{name}")
        assert list(net.layouts[name].columns) == (plane["names"] if isinstance(plane["names"], list) else [plane["names"]])
    labels = info.get("layoutsAxesLabels")
    if labels in (None, {}):
        assert "layoutsAxesLabels" not in net.options
    else:
        assert net.options["layoutsAxesLabels"] == {k: (v if isinstance(v, list) else [v]) for k, v in labels.items()}
    assert net.options.get("axesLabels") == info.get("axesLabels") or info.get("axesLabels") in (None, {})


@pytest.mark.parametrize("case", CASES)
def test_clusters_match_r(case):
    net, info = _run(case), _info(case)
    columns = info.get("clusterColumns")
    columns = [] if columns in (None, {}) else (columns if isinstance(columns, list) else [columns])
    assert net.scat.get("cluster_columns", []) == columns, case
    want = info.get("clusters") or {}
    want = {} if want == [] else want
    assert list(net.clusters) == list(want), case
    for key, entry in want.items():
        got = net.clusters[key]
        assert list(got.index) == (entry["rows"] if isinstance(entry["rows"], list) else [entry["rows"]])
        np.testing.assert_allclose(got.to_numpy(), np.array(entry["data"], dtype=float), rtol=1e-6, atol=1e-8, err_msg=f"{case}/{key}")
    shares = info.get("shares") or {}
    for key, rows in shares.items():
        np.testing.assert_allclose(net.scat["cluster_shares"][key].to_numpy(), np.array(rows, dtype=float), rtol=1e-8, atol=1e-9)
    means = info.get("clusterMeans")
    if means:
        for key, entry in means.items():
            got = net.cluster_means[key]
            assert list(got.columns) == (entry["cols"] if isinstance(entry["cols"], list) else [entry["cols"]])
            np.testing.assert_allclose(got.to_numpy(), np.array(entry["data"], dtype=float).reshape(got.shape), rtol=1e-6, atol=1e-8)
    else:
        assert not net.cluster_means
