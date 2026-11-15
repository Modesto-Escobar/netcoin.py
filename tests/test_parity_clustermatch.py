"""``cluster_match`` and ``cluster_agree`` against netCoin 2.1.49."""
import json
import warnings

import numpy as np
import pandas as pd
import pytest

from netcoin import add_clusters, cluster_agree, cluster_match, netCoin
from conftest import EXPECTED, read_expected, read_input

BASE = read_input("sa_base", strings=True).astype({"v": float, "x": float, "y": float})
NODES = BASE[["name", "v", "g"]]
XY = BASE[["x", "y"]].to_numpy()
NM = read_input("cm_nodes", strings=True)
PM = read_input("cm_pattern_nodes", strings=True)
CASES_DATA = read_input("sa_cases")
IDX, CW = CASES_DATA["idx"].to_numpy(), CASES_DATA["w"].to_numpy()
CASE_MODEL = read_input("cm_case_model")["case_model"].to_numpy()


def _net(patterns=False):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        net = netCoin(NODES, layout=XY, size="v")
    net.scat, net.clusters, net.cluster_means = {}, {}, None
    if patterns:
        net.scat["case_to_pattern"], net.scat["case_weight"] = IDX, CW
    return net


def nodes_net():
    """40 nodes with the clusterizations a, b, c, d (as ``add_clusters`` makes them) and two plain columns."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        net = add_clusters(_net(), NM[["a", "b", "c", "d"]].astype(int))
    net.nodes["m3"] = pd.Categorical(NM["m3"])
    net.nodes["txt"] = NM["txt"].to_numpy()
    return net


def patterns_net():
    """The patterns of 100 cases, with their groups as text (``Group: 1|Group: 2`` where a pattern spans several)."""
    net = _net(True)
    for c in ("g1", "g2", "g4"):
        net.nodes[c] = PM[c].to_numpy()
    net.scat["cluster_columns"] = ["g1", "g2", "g4"]
    return net


VEC = NM["vec"].to_numpy()
CASES = {
    "match_basic": lambda: cluster_match(nodes_net(), "a"),
    "match_noalign": lambda: cluster_match(nodes_net(), "a", align=False),
    "match_unweighted": lambda: cluster_match(nodes_net(), "a", weighted=False),
    "match_sort": lambda: cluster_match(nodes_net(), "a", sort="-Accuracy"),
    "match_sort2": lambda: cluster_match(nodes_net(), "a", sort=["F1", "Kappa"]),
    "match_digits": lambda: cluster_match(nodes_net(), "a", digits=1),
    "match_weights": lambda: cluster_match(nodes_net(), "a", weights="v"),
    "match_columns": lambda: cluster_match(nodes_net(), "a", columns=["b", "d", "m3", "txt"]),
    "match_vector": lambda: cluster_match(nodes_net(), VEC),
    "match_factor_model": lambda: cluster_match(nodes_net(), "m3", columns=["a", "b", "d"]),
    "match_patterns": lambda: cluster_match(patterns_net(), "g1"),
    "match_cases": lambda: cluster_match(patterns_net(), CASE_MODEL),
    "match_cases_w": lambda: cluster_match(patterns_net(), CASE_MODEL, weights=CW),
    "agree_basic": lambda: cluster_agree(nodes_net(), "a"),
    "agree_sort": lambda: cluster_agree(nodes_net(), "a", sort="-VI"),
    "agree_weights": lambda: cluster_agree(nodes_net(), "a", weights="v", digits=4),
    "agree_vector": lambda: cluster_agree(nodes_net(), VEC),
    "agree_patterns": lambda: cluster_agree(patterns_net(), "g1"),
    "agree_patterns_cases": lambda: cluster_agree(patterns_net(), CASE_MODEL),
    "agree_columns": lambda: cluster_agree(nodes_net(), "c", columns=["a", "b"]),
}


def _info(case):
    return json.loads((EXPECTED / f"cm_{case}_info.json").read_text(encoding="utf-8"))


def _run(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return CASES[case]()


# where several matchings of the categories are equally good R's clue and scipy may choose different ones:
# only what does not depend on the choice is compared for those classifications
INVARIANT = {"Classification", "k", "N", "Rows", "Accuracy"}


@pytest.mark.parametrize("case", CASES)
def test_the_table_matches_r(case):
    got = _run(case)
    want = read_expected(f"cm_{case}")
    assert list(got.columns) == list(want.columns), case
    assert len(got) == len(want)
    tied = set(got.attrs.get("ties") or [])
    for col in want.columns:
        if col == "Classification":
            assert got[col].tolist() == want[col].tolist()
            continue
        mask = ~got["Classification"].isin(tied).to_numpy() if col not in INVARIANT else np.ones(len(got), bool)
        np.testing.assert_allclose(got[col].astype(float)[mask], want[col].astype(float)[mask], rtol=1e-6, atol=1e-9,
                                   equal_nan=True, err_msg=f"{case}/{col}")


@pytest.mark.parametrize("case", CASES)
def test_the_details_match_r(case):
    got, info = _run(case), _info(case)
    # R names a vector after the expression it was given, python has no such name
    assert got.attrs["unit"] == info["unit"] and (got.attrs["model"] == info["model"] or "model" == got.attrs["model"])
    if info.get("kModel") not in (None, {}):
        assert got.attrs["kModel"] == info["kModel"]
    spanning = info.get("spanning") or {}
    assert (got.attrs.get("spanning") or {}) == {k: v for k, v in spanning.items()}
    tied = set(got.attrs.get("ties") or [])
    for name, entry in info["tables"].items():
        if name in tied:
            continue
        table = got.attrs["tables"][name]
        assert list(table.index) == (entry["rows"] if isinstance(entry["rows"], list) else [entry["rows"]]) or "agree" in case
        np.testing.assert_allclose(table.to_numpy(dtype=float), np.array(entry["data"], dtype=float).reshape(table.shape))
    for name, entry in (info.get("byClass") or {}).items():
        if name in tied:
            continue
        np.testing.assert_allclose(got.attrs["byClass"][name].to_numpy(), np.array(entry["data"], dtype=float).reshape(-1, 6),
                                   rtol=1e-9, atol=1e-12, equal_nan=True)
