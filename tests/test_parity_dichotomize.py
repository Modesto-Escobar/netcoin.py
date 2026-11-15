"""``dichotomize`` against netCoin 2.1.49."""
import pandas as pd
import pytest

from netcoin import dichotomize
from conftest import known_gap, read_expected, read_input

# case -> (input fixture, python arguments)
CASES = {
    "ess_age": ("ess", dict(variables="Age")),
    "ess_two": ("ess", dict(variables=["Gender", "Age"])),
    "kw": ("txt", dict(variables="kw2", sep=";")),
    "kw_min2": ("txt", dict(variables="kw2", sep=";", min_freq=2)),
    "kw_length2": ("txt", dict(variables="kw2", sep=";", length=2)),
    "kw_nosort": ("txt", dict(variables="kw2", sep=";", sort=False)),
    "kw_nas_null": ("txt", dict(variables="kw2", sep=";", nas=None)),
    "kw_nas_label": ("txt", dict(variables="kw2", sep=";", nas="Missing")),
    "kw_values": ("txt", dict(variables="kw2", sep=";", values=["red", "blue"])),
    "kw_add": ("txt", dict(variables="kw2", sep=";", add=True)),
    "kw_spaces": ("txt", dict(variables="keywords", sep=";")),
    "words_lower": ("txt", dict(variables="sentence", sep=" ", text="lower",
                                stopwords=["the", "a", "and"])),
    "two_seps": ("txt", dict(variables=["kw2", "colour"], sep=[";", ""])),
    "regex_pool": ("txt", dict(variables=["kw2", "sentence"], sep="[;, ]+",
                               regex=True, text="lower")),
}

_WORDS = 'python splits sep="" on whitespace; R keeps the whole value (e.g. "70 and +")'
_ORDER = ("column order: R breaks frequency ties with the locale collation, python with "
          "code points; and R prefixes the missing-value column with the variable name "
          "(\"kw2:None\"), python does not")
GAPS = {
    "ess_age": _WORDS,
    "ess_two": _WORDS,
    "kw": _ORDER,
    "kw_min2": _ORDER,
    "kw_length2": _ORDER,
    "kw_nosort": _ORDER,
    "kw_nas_null": _ORDER,
    "kw_nas_label": _ORDER,
    "kw_values": _ORDER,
    "kw_add": _ORDER,
    "kw_spaces": "R does not strip whitespace around the separator (\" green\"); python does",
    "words_lower": "python dichotomize() has no 'text' argument (lower/upper)",
    "two_seps": _ORDER,
    "regex_pool": "python dichotomize() has no 'text' argument nor pooling of several variables with regex",
}


@pytest.mark.parametrize("case", CASES)
def test_dichotomize(case, request, ess, txt):
    known_gap(request, GAPS, case)
    fixture, kwargs = CASES[case]
    data = {"ess": ess, "txt": txt}[fixture]
    expected = read_expected(f"dicho_{case}")

    actual = dichotomize(data, **kwargs)

    assert list(actual.columns) == list(expected.columns)
    exp_cols = [c for c in expected.columns if c in actual.columns]
    got = actual[exp_cols].reset_index(drop=True)
    for col in exp_cols:
        if pd.api.types.is_numeric_dtype(expected[col]):
            assert got[col].astype(float).tolist() == expected[col].astype(float).tolist(), col
        else:
            assert got[col].astype(str).tolist() == expected[col].astype(str).tolist(), col
