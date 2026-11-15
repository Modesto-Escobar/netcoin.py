"""Regression tests for the bugs fixed in phase 0."""
import numpy as np
import pandas as pd
import pytest

from netcoin import coin, edge_list, sim
from netcoin.simEdge import _looks_like_cooccurrence


@pytest.fixture
def incidences():
    return pd.DataFrame({
        "feature_1": [1, 0, 1, 1],
        "feature_2": [1, 1, 0, 1],
        "feature_3": [0, 1, 1, 0],
    })


def test_edge_list_accepts_incidence_data(incidences):
    edges = edge_list(incidences, procedures=["j", "d", "h"], criteria="j", min=0)
    assert list(edges.columns) == ["Source", "Target", "Jaccard", "dice", "Haberman"]
    assert len(edges) == 3


def test_edge_list_incidences_equal_coin_matrix(incidences):
    from_inc = edge_list(incidences, procedures="j", criteria="j", min=0)
    from_coin = edge_list(coin(incidences, print_n=False), procedures="j", criteria="j", min=0)
    pd.testing.assert_frame_equal(from_inc, from_coin)


def test_edge_list_does_not_add_coincidences_column(incidences):
    edges = edge_list(incidences, procedures="j", criteria="j", min=0, support=2)
    assert "coincidences" not in edges.columns
    assert len(edges) == 1  # only feature_1/feature_2 co-occur at least twice


def test_square_incidence_table_is_not_a_cooccurrence_matrix():
    rng = np.random.default_rng(1)
    square = pd.DataFrame(rng.integers(0, 2, (4, 4)), columns=list("abcd"))
    assert not _looks_like_cooccurrence(square)
    # and sim treats it as incidences: n = 4 cases
    assert sim(square, "f").sum() == square.to_numpy().sum()


def test_coin_matrix_is_detected(incidences):
    assert _looks_like_cooccurrence(coin(incidences, print_n=False))
