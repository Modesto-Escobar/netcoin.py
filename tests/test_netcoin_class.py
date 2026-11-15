"""The ``netCoin`` class: inputs, link filtering (as in R) and rebuilding."""
import json

import numpy as np
import pandas as pd
import pytest

from netcoin import netCoin
from netcoin.webgraph import r_json


def test_nodes_come_from_the_links_in_order_of_appearance():
    net = netCoin(links={"Source": ["b", "a", "c"], "Target": ["a", "c", "d"]})
    assert net.nodes["name"].tolist() == ["b", "a", "c", "d"]


def test_self_links_and_unknown_endpoints_are_dropped():
    nodes = {"name": ["a", "b", "c"]}
    links = {"Source": ["a", "a", "b", "x"], "Target": ["b", "a", "c", "c"], "w": [1, 2, 3, 4]}
    net = netCoin(nodes, links)
    assert net.links["Source"].tolist() == ["a", "b"] and net.links["w"].tolist() == [1, 3]
    assert net.lcount() == 2


def test_all_links_dropped_warns_and_leaves_a_network_without_links():
    with pytest.warns(UserWarning, match="no row"):
        net = netCoin({"name": ["a", "b"]}, {"Source": ["x"], "Target": ["y"]})
    assert net.lcount() == 0 and net.links is None
    assert "linkSource" not in net.options


def test_dataframes_numpy_layout_and_missing_values():
    nodes = pd.DataFrame({"name": ["A", "B", "C"], "size": [1.0, 2.0, np.nan]})
    links = pd.DataFrame({"Source": ["A", "B"], "Target": ["B", "C"], "w": [0.5, 0.7]})
    net = netCoin(nodes, links, size="size", layout=np.array([[0, 0], [1, 1], [2, 0]]))
    assert net.nodes["fx"].tolist() == [0, 1, 2]
    assert net.nodes["size"].isna().tolist() == [False, False, True]
    from netcoin.network import networkJSON
    json.dumps(r_json(networkJSON(net)), allow_nan=False)      # valid JSON: no NaN anywhere


def test_numeric_links_are_one_based_positions():
    net = netCoin({"name": ["a", "b", "c"]}, {"Source": [1, 2], "Target": [2, 3]})
    assert net.links["Source"].tolist() == ["a", "b"] and net.links["Target"].tolist() == ["b", "c"]


def test_duplicated_names_are_made_unique_with_a_warning():
    with pytest.warns(UserWarning, match="unique"):
        net = netCoin({"name": ["a", "a", "b"]})
    assert net.nodes["name"].tolist() == ["a", "a.1", "b"]


def test_frames_roundtrip():
    net = netCoin({"name": ["a", "b"], "v": [1, 2]}, {"Source": ["a"], "Target": ["b"]})
    assert list(net.nodes_frame().columns) == ["name", "v"]
    assert net.links_frame()[["Source", "Target"]].values.tolist() == [["a", "b"]]
    assert list(netCoin({"name": ["a"]}).links_frame().columns) == ["Source", "Target"]


def test_scenarios_is_an_integer_when_it_is_one():
    net = netCoin({"name": ["a"]}, scenarios=17.0)
    assert net.options["scenarios"] == 17 and isinstance(net.options["scenarios"], int)


def test_str_shows_the_first_rows_of_nodes_and_links():
    nodes = {"name": [f"n{i}" for i in range(8)], "v": list(range(8))}
    links = {"Source": [f"n{i}" for i in range(7)], "Target": [f"n{i + 1}" for i in range(7)]}
    text = str(netCoin(nodes, links, main="Title"))
    assert text.startswith("Title: Title") and "n5" in text and "n7" not in text
    assert text.count("...") == 2


def test_repr_is_a_short_summary():
    net = netCoin({"name": ["a", "b"], "v": [1, 2]}, {"Source": ["a"], "Target": ["b"], "w": [1]})
    assert repr(net) == "<netCoin object>\nNodes: 2  (name, v)\nLinks: 1  (Source, Target, w)"


def test_nodes_or_links_are_needed():
    with pytest.raises(ValueError, match="nodes or links"):
        netCoin()


def test_a_netcoin_cannot_be_the_nodes_of_another():
    with pytest.raises(TypeError, match="add"):
        netCoin(netCoin({"name": ["a"]}))


# ------------------------------------------------------- camelCase names ----

def test_camel_case_names_still_work_with_a_deprecation_warning():
    with pytest.warns(DeprecationWarning, match="showArrows"):
        net = netCoin({"name": ["a"]}, showArrows=True, helpOn=True)
    assert net.options["showArrows"] is True and net.options["helpOn"] is True


def test_order_a_and_order_d_became_sort_and_decreasing():
    with pytest.warns(DeprecationWarning, match="orderD"):
        net = netCoin({"name": ["a"], "v": [1]}, orderD="v")
    assert net.options["nodeOrder"] == "v" and net.options["decreasing"] is True


def test_unknown_and_repeated_arguments_fail():
    with pytest.raises(TypeError, match="unexpected"):
        netCoin({"name": ["a"]}, colour="red")
    with pytest.raises(TypeError, match="multiple"):
        with pytest.warns(DeprecationWarning):
            netCoin({"name": ["a"]}, show_arrows=True, showArrows=True)


def test_from_matrix_and_its_deprecated_name():
    matrix = [[1, 2, 0.5], [2, 3, 0.7]]
    net = netCoin.from_matrix(matrix, nodes={"name": ["a", "b", "c"]}, link_attr_names=["Source", "Target", "w"])
    assert net.links["Source"].tolist() == ["a", "b"] and net.links["w"].tolist() == [0.5, 0.7]
    with pytest.warns(DeprecationWarning, match="fromMatrix"):
        old = netCoin.fromMatrix(matrix, nodes={"name": ["a", "b", "c"]}, linkAttrNames=["Source", "Target", "w"])
    assert old.links.equals(net.links)
    with pytest.raises(ValueError, match="same length"):
        netCoin.from_matrix([[1, 2], [1, 2, 3]], nodes={"name": ["a", "b", "c"]})
