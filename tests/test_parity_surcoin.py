"""``sur_coin`` and ``cob_coin`` against netCoin 2.1.49."""
import json
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from netcoin import cob_coin, sur_coin
from conftest import (EXPECTED, as_list, assert_matrix_close, compare_frame, compare_options,
                      read_expected, read_input, read_matrix)

ESS = read_input("ess_cweight", strings=True).astype({"cweight": float})
ESS_M = read_input("ess_metric", strings=True).astype({"cweight": float, "income": float, "edu": float})
ESS_CC = read_input("ess_complete", strings=True).astype({"cweight": float})
NODES = read_input("sur_nodes", strings=True)
NODES_M = read_input("sur_nodes_metric", strings=True)
VARS3 = ["Gender", "Age", "Social participation"]
VARS4 = VARS3 + ["Political participation"]
LAYOUT = np.column_stack([np.arange(6), np.arange(5, -1, -1)]).astype(float)

# python call of each R call (the R calls are in tests/reference/make_reference.R)
CASES = {
    "default": lambda: sur_coin(ESS, VARS3),
    "one_dicho": lambda: sur_coin(ESS, ["Gender"], dichotomies="Social participation",
                                  value_dicho="Yes", max_link=1),
    "dichos": lambda: sur_coin(ESS, ["Gender", "Age"], dichotomies=["Social participation", "Political participation"],
                               value_dicho="Yes", max_link=1),
    "dicho_inside": lambda: sur_coin(ESS, ["Social participation", "Gender", "Age"],
                                     dichotomies="Social participation", value_dicho="Yes", max_link=1),
    "links": lambda: sur_coin(ESS, VARS3, frequency=True, percentages=False,
                              procedures=["Haberman", "Jaccard"], lwidth="Jaccard", lcolor="Haberman",
                              criteria="Jaccard", min_link=0.2, language="es", sort=True, decreasing=False,
                              sort_links="Jaccard", main="Survey"),
    "minimum": lambda: sur_coin(ESS, VARS3, minimum=150, max_link=1),
    "weight": lambda: sur_coin(ESS, ["Gender", "Age", "cweight"], weight="cweight", max_link=1),
    "weight_vector": lambda: sur_coin(ESS_CC, VARS3, weight=ESS_CC["cweight"].to_numpy(), max_link=1),
    "subsample": lambda: sur_coin(ESS, VARS3, subsample=True, max_link=1),
    "pairwise": lambda: sur_coin(ESS, VARS3, pairwise=True, max_link=1),
    "nodes": lambda: sur_coin(ESS, VARS3, nodes=NODES, color="grp", max_link=1),
    "layout_mc": lambda: sur_coin(ESS, VARS3, layout="mc", max_link=1),
    "layout_pc": lambda: sur_coin(ESS, VARS3, layout="pc"),
    "layout_two": lambda: sur_coin(ESS, VARS3, layout=["mc", "pc"]),
    "nodes_mc": lambda: sur_coin(ESS, VARS3, nodes=NODES, layout="mc", max_link=1),
    "layout_matrix": lambda: sur_coin(ESS, VARS3, nodes=NODES, layout=LAYOUT, max_link=1),
    "commonlabel": lambda: sur_coin(ESS, ["Gender", "Age"], commonlabel="Gender", max_link=1),
    "exogenous": lambda: sur_coin(ESS, VARS3, exogenous=["Gender", "Age"], max_link=1),
    "exogenous_filter": lambda: sur_coin(ESS, VARS3, exogenous="Gender", link_filter="Haberman > 0",
                                         max_link=1),
    "spanish": lambda: sur_coin(ESS, VARS3, language="es", lwidth="Haberman", ltext="Z"),
    "metric": lambda: sur_coin(ESS_M, ["Gender", "Social participation"], metric=["income", "edu"], max_link=1),
    "metric_pearson": lambda: sur_coin(ESS_M, VARS3, metric="income", procedures=["Pearson", "Haberman"],
                                       criteria="Pearson", min_link=0.05, max_link=np.inf),
    "metric_pairwise": lambda: sur_coin(ESS_M, ["Gender", "Age"], metric=["income", "edu"], pairwise=True,
                                        max_link=1),
    "metric_weight": lambda: sur_coin(ESS_M, ["Gender", "Age", "cweight"], metric="edu", weight="cweight",
                                      max_link=1),
    "metric_nodes": lambda: sur_coin(ESS_M, VARS3, metric=["income", "edu"], nodes=NODES_M, color="grp",
                                     max_link=1),
    "metric_nopct": lambda: sur_coin(ESS_M, VARS3, metric="edu", percentages=False, frequency=True,
                                     max_link=1),
    "metric_exogenous": lambda: sur_coin(ESS_M, VARS3, metric="income", exogenous=["Gender", "income"],
                                         max_link=1),
    "cob_default": lambda: cob_coin(ESS, VARS4),
    "cob_args": lambda: cob_coin(ESS, VARS4, degree=45, significance=0.1, language="es", main="CoWeb"),
}

METRICS = {"metric": ["income", "edu"], "metric_pearson": ["income"], "metric_pairwise": ["income", "edu"],
           "metric_weight": ["edu"], "metric_nodes": ["income", "edu"], "metric_nopct": ["edu"],
           "metric_exogenous": ["income"]}
_METRIC_P = ("R tests the correlations with t = r / sqrt((1 - r) / (n - 2)) instead of "
             "r * sqrt((n - 2) / (1 - r^2)); python uses the usual one (see reference/README.md)")
# R fills the layout by position, so it does not follow the nodes when they are reordered
GAPS = {"nodes_mc": "R puts the coordinates of the correspondence analysis in the order of the "
                    "incidences, not of the (sorted) nodes: python aligns them by name"}


def _ref(case):
    return json.loads((EXPECTED / f"net_sur_{case}.json").read_text(encoding="utf-8"))


def _mute(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return CASES[case]()


def _links_without_metric_p(frame, metric):
    """The p-value of the metric links is left out: python and R differ there."""
    frame = frame.copy()
    if "p(Z)" in frame:
        is_metric = frame["Source"].isin(metric) | frame["Target"].isin(metric)
        frame.loc[is_metric.to_numpy(), "p(Z)"] = np.nan
    return frame


@pytest.mark.parametrize("case", CASES)
def test_nodes_and_options_match_r(case, request):
    if case in GAPS:
        request.applymarker(pytest.mark.xfail(strict=True, reason=GAPS[case]))
    net, ref = _mute(case), _ref(case)
    compare_frame(net.nodes, read_expected(f"net_sur_{case}_nodes"), f"{case}/nodes",
                  up_to_sign=["fx", "fy"] if "layout" in case else ())
    compare_options(net.options, ref["options"], case)
    if "layouts" in ref:
        assert list(net.layouts) == list(ref["layouts"])


@pytest.mark.parametrize("case", CASES)
def test_links_match_r(case):
    net, ref = _mute(case), _ref(case)
    assert net.lcount() == ref["nlinks"], f"{case}: links"
    if not ref["nlinks"]:
        return
    expected = read_expected(f"net_sur_{case}_links")
    actual = net.links
    if case in METRICS:
        actual = _links_without_metric_p(actual, METRICS[case])
        expected = _links_without_metric_p(expected, METRICS[case])
    compare_frame(actual, expected, f"{case}/links")


def _indicator(frame, label):
    """The 0/1 column of a category ``"variable:category"`` (or the metric column itself)."""
    if label in frame.columns:
        return frame[label].astype(float).to_numpy()
    variable, category = label.split(":", 1)
    return (frame[variable] == category).astype(float).to_numpy()


def test_the_metric_links_use_the_usual_test_of_the_correlation():
    """What R gets wrong (see GAPS): the Haberman residual is r * sqrt(n), the p-value that of
    t = r * sqrt((n - 2) / (1 - r^2)) with n - 2 degrees of freedom; n is the number of cases of
    the survey and r is computed with the cases that have every metric variable."""
    net = _mute("metric")
    metric = ["income", "edu"]
    cases = ESS_M.dropna(subset=metric)
    n = len(ESS_M)
    rows = net.links[net.links["Target"].isin(metric)]
    assert len(rows) >= 6
    for _, row in rows.iterrows():
        r = np.corrcoef(_indicator(cases, row["Source"]), _indicator(cases, row["Target"]))[0, 1]
        t = r * np.sqrt((n - 2) / (1 - r ** 2))
        assert row["Haberman"] == pytest.approx(r * np.sqrt(n), rel=1e-9)
        assert row["p(Z)"] == pytest.approx(stats.t.sf(t, n - 2), rel=1e-9, abs=1e-15)


# ------------------------------------------------------------ the Coin object ----

@pytest.mark.parametrize("case,call", [
    ("sur", lambda: sur_coin(ESS, ["Gender", "Age"], coin=True)),
    ("sur_weight", lambda: sur_coin(ESS, ["Gender", "Age", "cweight"], weight="cweight", coin=True)),
])
def test_coin_object_matches_r(case, call):
    result = call()
    assert_matrix_close(result.matrix, read_matrix(f"coin_{case}"), label=case)
    attrs = read_expected(f"coin_{case}_attrs")
    assert result.n == pytest.approx(attrs["n"][0])
