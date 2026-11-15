"""``net_corr``, ``d_net_corr`` and ``layout_fact`` against netCoin 2.1.49."""
import json
import warnings

import numpy as np
import pandas as pd
import pytest

from netcoin import d_net_corr, layout_fact, net_corr
from conftest import (EXPECTED, compare_frame, compare_options, read_expected, read_input)

MT = read_input("corr_mt")
MT_W = read_input("corr_mt_w")
MT_NA = read_input("corr_mt_na")
MT_NAW = MT_NA.assign(w=MT_W["w"])
SMALL = read_input("corr_small")
SMALL_TIES = read_input("corr_small_ties")
SMALL_NA = read_input("corr_small_na")
BLOCKS = read_input("corr_blocks")
NODES = read_input("corr_nodes", strings=True)
NODES_ES = read_input("corr_nodes_es", strings=True)
SEQUENCE = np.round(np.arange(20, 101, 5) / 100, 2)

# python call of each R call (the R calls are in tests/reference/make_reference.R)
CASES = {
    "default": lambda: net_corr(MT),
    "spearman": lambda: net_corr(MT, method="spearman"),
    "kendall": lambda: net_corr(MT, method="kendall"),
    "partial": lambda: net_corr(MT, method="partial"),
    "dcor": lambda: net_corr(MT, method="dcor"),
    "multi": lambda: net_corr(MT, method=["pearson", "spearman", "kendall"], criteria="p(kendall)"),
    "small": lambda: net_corr(SMALL, method=["spearman", "kendall"], criteria="p", max_link=1),
    "small_one": lambda: net_corr(SMALL, method=["spearman", "kendall"], link_bipolar=False, max_link=1),
    "small_ties": lambda: net_corr(SMALL_TIES, method=["spearman", "kendall"], max_link=1),
    "small_na": lambda: net_corr(SMALL_NA, pairwise=True, method=["spearman", "kendall", "pearson"], max_link=1),
    "weight": lambda: net_corr(MT_W, weight="w", method=["pearson", "spearman", "partial"]),
    "weight_vector": lambda: net_corr(MT, weight=MT_W["w"].to_numpy(), method="spearman"),
    "na_listwise": lambda: net_corr(MT_NA),
    "na_pairwise": lambda: net_corr(MT_NA, pairwise=True, method=["pearson", "spearman", "kendall", "dcor"]),
    "na_pairwise_weight": lambda: net_corr(MT_NAW, weight="w", pairwise=True, method=["pearson", "spearman"]),
    "na_partial": lambda: net_corr(MT_NA, pairwise=True, method="partial"),
    "value": lambda: net_corr(MT, criteria="value", min_link=0.5),
    "value_one": lambda: net_corr(MT, criteria="value", min_link=0.5, max_link=0.9, link_bipolar=False),
    "one_sided": lambda: net_corr(MT, link_bipolar=False),
    "bonferroni": lambda: net_corr(MT, bonferroni=True, max_link=0.05),
    "sort": lambda: net_corr(MT_NA, pairwise=True, sort=True, decreasing=False, sort_links="value",
                             decreasing_links=False),
    "sort_decreasing": lambda: net_corr(MT_NA, pairwise=True, sort=True, sort_links="z"),
    "nomeans": lambda: net_corr(MT, means=False, frequency=True),
    "nodes": lambda: net_corr(MT, nodes=NODES_ES, language="es", color="grp"),
    "layout_pc": lambda: net_corr(MT, layout="pc"),
    "layout_pc_nodes": lambda: net_corr(MT, nodes=NODES, layout="pc"),
    "options": lambda: net_corr(MT, lwidth="z", lcolor="p", size="mean", color="std", main="Cars"),
    "d_default": lambda: d_net_corr(BLOCKS, sequence=SEQUENCE),
    "d_pc": lambda: d_net_corr(BLOCKS, factorial="pc", sequence=SEQUENCE, method="spearman"),
    "d_fact": lambda: d_net_corr(BLOCKS, factorial=["pc", "nf", "vf"], sequence=SEQUENCE),
    "d_all": lambda: d_net_corr(BLOCKS, components=False, sequence=np.round(np.arange(3, 10) / 10, 1), means=False),
    "d_back": lambda: d_net_corr(BLOCKS, backcomponents=True, sequence=SEQUENCE, text_filter=(0.5, 0.1)),
}
GAPS = {
    "layout_pc_nodes": "R lays the coordinates in the order of the correlation matrix, not of the "
                       "nodes (merged with the table of nodes): python aligns them by name",
}
UP_TO_SIGN = ["fx", "fy"]
# R's toIgraph looks for the options lwidth, lweight... (they are called linkWidth...) and never
# passes the weights of the links, so its weighted degree equals the degree; python weights it
SKIP = {case: ["wdegree", "eigen"] for case in ("d_default", "d_pc", "d_fact", "d_all", "d_back")}
# factor loadings come from an optimiser: a looser tolerance
LOOSE = {"d_fact": 2e-3}


def _ref(case):
    return json.loads((EXPECTED / f"net_nc_{case}.json").read_text(encoding="utf-8"))


def _run(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return CASES[case]()


@pytest.mark.parametrize("case", CASES)
def test_nodes_and_options_match_r(case, request):
    if case in GAPS:
        request.applymarker(pytest.mark.xfail(strict=True, reason=GAPS[case]))
    net, ref = _run(case), _ref(case)
    nodes = net.nodes
    compare_frame(nodes, read_expected(f"net_nc_{case}_nodes"), f"{case}/nodes", up_to_sign=UP_TO_SIGN,
                  rtol=LOOSE.get(case, 1e-8), skip=SKIP.get(case, ()))
    compare_options(net.options, ref["options"], case)


@pytest.mark.parametrize("case", CASES)
def test_links_match_r(case):
    net, ref = _run(case), _ref(case)
    assert net.lcount() == ref["nlinks"], f"{case}: links"
    if ref["nlinks"]:
        compare_frame(net.links, read_expected(f"net_nc_{case}_links"), f"{case}/links", rtol=1e-7)


@pytest.mark.parametrize("method", ["pc", "nf", "vf", "of"])
def test_layout_fact_matches_r(method):
    got = layout_fact(BLOCKS, method)
    expected = read_expected(f"layoutfact_{method}").set_index(".row")
    assert list(got.index) == list(expected.index)
    a, e = got.to_numpy(), expected.to_numpy(dtype=float)
    for j in range(2):          # the sign of a component is arbitrary
        sign = 1 if np.abs(a[:, j] - e[:, j]).sum() <= np.abs(a[:, j] + e[:, j]).sum() else -1
        np.testing.assert_allclose(sign * a[:, j], e[:, j], atol=2e-3, err_msg=f"{method}/{j}")


def test_layout_fact_several_methods_and_names():
    planes = layout_fact(BLOCKS, ["pc", "nf", "vf", "of"])
    names = json.loads((EXPECTED / "layoutfact_names.json").read_text(encoding="utf-8"))["names"]
    assert list(planes) == names
    assert layout_fact(BLOCKS, "null") is None
    with pytest.warns(UserWarning, match="not implemented"):
        assert layout_fact(BLOCKS, ["pc", "zz"]).shape == (8, 2)
