"""``Coin`` and ``coin`` behaviour that is specific to the python port."""
import numpy as np
import pandas as pd
import pytest

from netcoin import Coin, coin, edge_list, expected_list, sim


@pytest.fixture
def incidences():
    return pd.DataFrame({
        "feature_1": [1, 0, 1, 1],
        "feature_2": [1, 1, 0, 1],
        "feature_3": [0, 1, 1, 0],
    })


def test_coin_returns_a_coin_with_n(incidences):
    c = coin(incidences)
    assert isinstance(c, Coin)
    assert c.n == 4
    assert c.labels == ["feature_1", "feature_2", "feature_3"]
    assert c.frequencies.tolist() == [3, 3, 2]
    assert c.summary().tolist() == [0.75, 0.75, 0.5]
    assert len(c) == 3


def test_coin_repr_is_like_r(incidences):
    text = repr(coin(incidences))
    assert text.splitlines()[0] == "n= 4"
    assert "feature_3" in text


def test_coin_rejects_non_dichotomous_data():
    with pytest.raises(ValueError, match="dichotomous"):
        coin(pd.DataFrame({"a": [0, 2, 1]}))


def test_coin_no_variables_left(incidences):
    with pytest.raises(ValueError, match="no variables left"):
        coin(incidences, minimum=10)


def test_coin_accepts_bool_numpy_and_nullable_columns():
    as_bool = pd.DataFrame({"a": [True, False, True], "b": [True, True, False]})
    as_int64 = as_bool.astype("Int64")
    as_array = as_bool.to_numpy()
    for data in (as_bool, as_int64, as_array):
        assert coin(data).matrix.iloc[0, 1] == 1


def test_coin_sparse_columns():
    dense = pd.DataFrame({"a": [1, 0, 1], "b": [1, 1, 0]})
    sparse = dense.astype(pd.SparseDtype("int", 0))
    assert coin(sparse).matrix.equals(coin(dense).matrix)


def test_coin_weight_length_is_checked(incidences):
    with pytest.raises(ValueError, match="weight"):
        coin(incidences, weight=[1, 2])


def test_total_is_not_duplicated(incidences):
    c = coin(incidences, total=True)
    assert c.labels == ["Total", "feature_1", "feature_2", "feature_3"]
    assert c.matrix.loc["Total", "feature_1"] == 3


def test_total_with_weights_warns(incidences):
    with pytest.warns(UserWarning, match="total"):
        c = coin(incidences, total=True, weight=[1, 1, 1, 1])
    assert "Total" not in c.labels


def test_subsample_drops_cases_without_events():
    data = pd.DataFrame({"a": [1, 0, 0, 1], "b": [0, 0, 0, 1]})
    assert coin(data, subsample=True).n == 2


def test_sort_is_stable_for_ties(incidences):
    c = coin(incidences, sort=True)
    assert c.labels == ["feature_1", "feature_2", "feature_3"]
    assert coin(incidences, sort=True, decreasing=False).labels[0] == "feature_3"


def test_from_matrix_accepts_a_lower_triangle():
    lower = pd.DataFrame([[3, ""], [2, 3]], index=list("ab"), columns=list("ab"))
    c = Coin.from_matrix(lower, n=4)
    assert c.matrix.loc["a", "b"] == 2 and c.n == 4
    assert sim(c, "Jaccard").loc["a", "b"] == pytest.approx(2 / 4)


def test_coin_requires_square_labelled_matrix():
    with pytest.raises(ValueError):
        Coin(pd.DataFrame([[1, 2, 3]]), 3)
    with pytest.raises(ValueError):
        Coin(pd.DataFrame([[1, 0], [0, 1]], index=["a", "b"], columns=["a", "c"]), 2)


# ------------------------------------------------- sim / edge_list inputs ----

def test_square_incidence_table_is_taken_as_incidences():
    rng = np.random.default_rng(1)
    square = pd.DataFrame(rng.integers(0, 2, (4, 4)), columns=list("abcd"))
    got = sim(square, "coincidences")
    sums = square.sum()
    assert np.diag(got).tolist() == sums[sums >= 1].tolist()
    assert sim(coin(square), "coincidences").equals(got)


def test_sim_rejects_other_types():
    with pytest.raises(TypeError):
        sim([[1, 0], [0, 1]])


def test_edge_list_accepts_incidences_and_coin(incidences):
    kwargs = dict(procedures=["Jaccard", "dice", "Haberman"], criteria="Jaccard", min=0)
    from_inc = edge_list(incidences, **kwargs)
    from_coin = edge_list(coin(incidences), **kwargs)
    assert list(from_inc.columns) == ["Source", "Target", "Jaccard", "dice", "Haberman"]
    pd.testing.assert_frame_equal(from_inc, from_coin)


def test_edge_list_does_not_add_coincidences_column(incidences):
    edges = edge_list(incidences, "Jaccard", criteria="Jaccard", min=0, support=2)
    assert "coincidences" not in edges.columns
    assert len(edges) == 1  # only feature_1/feature_2 coincide at least twice


def test_edge_list_bonferroni_counts_variables():
    rng = np.random.default_rng(3)
    data = pd.DataFrame(rng.integers(0, 2, (50, 6)), columns=list("abcdef"))
    plain = edge_list(data, max=0.5)
    strict = edge_list(data, max=0.5, bonferroni=True)
    assert len(strict) <= len(plain)
    assert strict["p(Z)"].max() <= 0.5 / 15


def test_edge_list_sort_is_stable_and_descending(incidences):
    edges = edge_list(incidences, ["Jaccard"], criteria="Jaccard", min=0, sort="Jaccard")
    assert edges["Jaccard"].is_monotonic_decreasing


def test_expected_list_needs_a_coin(incidences):
    with pytest.raises(TypeError):
        expected_list(incidences)
