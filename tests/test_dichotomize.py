"""Behaviour of ``dichotomize`` that is specific to the python port."""
import numpy as np
import pandas as pd
import pytest

from netcoin import coin, dichotomize


@pytest.fixture
def df():
    return pd.DataFrame(
        {"kw": ["a;b", "b", None, "a;;c", ""], "n": [1, 2.0, 2, np.nan, 1]},
        index=list("vwxyz"),
    )


def test_not_a_dataframe():
    with pytest.raises(TypeError):
        dichotomize([1, 2], "x")


def test_unknown_variable(df):
    with pytest.raises(KeyError):
        dichotomize(df, "nope")


def test_bad_text_option(df):
    with pytest.raises(ValueError):
        dichotomize(df, "kw", text="title")


def test_index_is_kept_and_values_are_int(df):
    out = dichotomize(df, "kw", sep=";")
    assert list(out.index) == list(df.index)
    assert all(np.issubdtype(t, np.integer) for t in out.dtypes)


def test_missing_and_blank_rows_are_flagged(df):
    out = dichotomize(df, "kw", sep=";")
    assert out["kw:None"].tolist() == [0, 0, 1, 0, 1]


def test_empty_tokens_are_dropped(df):
    out = dichotomize(df, "kw", sep=";", nas=None)
    assert sorted(out.columns) == ["a", "b", "c"]
    assert out.loc["y"].to_dict() == {"a": 1, "b": 0, "c": 1}


def test_numbers_are_converted_like_as_character(df):
    out = dichotomize(df, "n", nas=None)
    assert sorted(out.columns) == ["1", "2"]  # 2.0 -> "2", not "2.0"


def test_add_keeps_original_columns(df):
    out = dichotomize(df, "kw", sep=";", add=True)
    assert list(out.columns[:2]) == ["kw", "n"]
    assert out["kw"].tolist() == df["kw"].tolist()


def test_nothing_to_dichotomize_returns_the_data():
    data = pd.DataFrame({"x": [None, "", " "]})
    out = dichotomize(data, "x")
    pd.testing.assert_frame_equal(out, data)


def test_stopwords_only_apply_with_text():
    data = pd.DataFrame({"s": ["the cat", "a cat"]})
    plain = dichotomize(data, "s", sep=" ", stopwords=["the"], nas=None)
    assert "the" in plain.columns
    cleaned = dichotomize(data, "s", sep=" ", text="lower", stopwords=["The", "a"], nas=None)
    assert list(cleaned.columns) == ["cat"]


def test_text_removes_punctuation_and_digits():
    data = pd.DataFrame({"s": ["Hi, 2 you!", "hi"]})
    out = dichotomize(data, "s", sep=" ", text="lower", nas=None)
    assert list(out.columns) == ["hi", "you"]


def test_accents_and_case_order_ties_like_r():
    data = pd.DataFrame({"s": ["Zeta;álamo;alamo;Alamo;beta"]})
    out = dichotomize(data, "s", sep=";", nas=None)
    assert list(out.columns) == ["alamo", "Alamo", "álamo", "beta", "Zeta"]


def test_sparse_has_no_na_column_and_sparse_dtype(df):
    out = dichotomize(df, "kw", sep=";", sparse=True)
    assert "kw:None" not in out.columns
    assert all(isinstance(t, pd.SparseDtype) for t in out.dtypes)
    assert out.sparse.to_dense()["a"].tolist() == [1, 0, 0, 1, 0]


def test_sparse_add_warns(df):
    with pytest.warns(UserWarning):
        dichotomize(df, "kw", sep=";", sparse=True, add=True)


def test_separator_sequence_is_recycled():
    data = pd.DataFrame({"a": ["x;y"], "b": ["x;z"]})
    out = dichotomize(data, ["a", "b"], sep=[";"], nas=None)
    assert sorted(out.columns) == ["a:x", "a:y", "b:x", "b:z"]


def test_feeds_coin():
    data = pd.DataFrame({"s": ["a;b", "a", "b;c", "a;c"]})
    inc = dichotomize(data, "s", sep=";", nas=None)
    cooc = coin(inc, print_n=False)
    assert cooc.loc["a", "a"] == 3
