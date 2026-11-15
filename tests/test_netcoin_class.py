"""The ``netCoin`` class: inputs, link filtering (as in R) and rebuilding."""
import json

import numpy as np
import pandas as pd
import pytest

from netcoin import netCoin


def test_nodes_come_from_the_links_in_order_of_appearance():
    net = netCoin(links={"Source": ["b", "a", "c"], "Target": ["a", "c", "d"]})
    assert net.nodes["name"] == ["b", "a", "c", "d"]


def test_self_links_and_unknown_endpoints_are_dropped():
    nodes = {"name": ["a", "b", "c"]}
    links = {"Source": ["a", "a", "b", "x"], "Target": ["b", "a", "c", "c"], "w": [1, 2, 3, 4]}
    net = netCoin(nodes, links)
    assert net.links["Source"] == ["a", "b"] and net.links["w"] == [1, 3]
    assert net.lcount() == 2


def test_all_links_dropped_warns_and_leaves_a_network_without_links():
    with pytest.warns(UserWarning, match="no row"):
        net = netCoin({"name": ["a", "b"]}, {"Source": ["x"], "Target": ["y"]})
    assert net.lcount() == 0 and net.links == {}


def test_dataframes_numpy_layout_and_missing_values():
    nodes = pd.DataFrame({"name": ["A", "B", "C"], "size": [1.0, 2.0, np.nan]})
    links = pd.DataFrame({"Source": ["A", "B"], "Target": ["B", "C"], "w": [0.5, 0.7]})
    net = netCoin(nodes, links, size="size", layout=np.array([[0, 0], [1, 1], [2, 0]]))
    assert net.nodes["size"] == [1.0, 2.0, None]
    assert net.nodes["fx"] == [0, 1, 2]
    json.dumps(net.getNet(), allow_nan=False)      # valid JSON: no NaN anywhere


def test_frames_roundtrip():
    net = netCoin({"name": ["a", "b"], "v": [1, 2]}, {"Source": ["a"], "Target": ["b"]})
    assert list(net.nodes_frame().columns) == ["name", "v"]
    assert net.links_frame()[["Source", "Target"]].values.tolist() == [["a", "b"]]


def test_scenarios_is_an_integer_when_it_is_one():
    net = netCoin({"name": ["a"]}, scenarios=17.0)
    assert net.options["scenarios"] == 17 and isinstance(net.options["scenarios"], int)
