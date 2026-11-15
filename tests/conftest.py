"""Shared helpers for the parity tests against the R package netCoin.

The reference values live in ``tests/reference`` and are produced by
``tests/reference/make_reference.R``; they are committed, so R is not needed
to run the tests.
"""
from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REFERENCE = Path(__file__).parent / "reference"
INPUTS = REFERENCE / "inputs"
EXPECTED = REFERENCE / "expected"

_NA = ["NA", "NaN"]


def read_input(name: str, strings: bool = False) -> pd.DataFrame:
    """Read an input table. ``strings=True`` keeps empty strings as ``""``."""
    path = INPUTS / f"{name}.csv"
    if strings:
        return pd.read_csv(path, encoding="utf-8", keep_default_na=False,
                           na_values=_NA, dtype=str)
    return pd.read_csv(path, encoding="utf-8")


def read_expected(name: str, strings: bool = False) -> pd.DataFrame | None:
    """Read an expected table; ``None`` when R returned an empty result."""
    path = EXPECTED / f"{name}.csv"
    try:
        return pd.read_csv(path, encoding="utf-8", keep_default_na=False,
                           na_values=_NA, dtype=str if strings else None)
    except pd.errors.EmptyDataError:
        return None


def read_matrix(name: str) -> pd.DataFrame:
    """Read a labelled square matrix written by ``write_matrix`` in R."""
    df = read_expected(name)
    df = df.set_index(".row")
    df.index.name = None
    return df


@pytest.fixture(scope="session")
def manifest() -> dict:
    return json.loads((REFERENCE / "manifest.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="session")
def galapagos() -> pd.DataFrame:
    return read_input("galapagos")


@pytest.fixture(scope="session")
def dice() -> pd.DataFrame:
    return read_input("dice")


@pytest.fixture(scope="session")
def txt() -> pd.DataFrame:
    return read_input("txt", strings=True)


@pytest.fixture(scope="session")
def ess() -> pd.DataFrame:
    return read_input("ess", strings=True)


def assert_matrix_close(actual, expected, rtol=1e-9, atol=1e-12, label=""):
    """Same labels (in order) and numerically equal values (NaN == NaN)."""
    actual = pd.DataFrame(actual)
    assert list(actual.index) == list(expected.index), f"{label} row labels differ"
    assert list(actual.columns) == list(expected.columns), f"{label} column labels differ"
    a = actual.to_numpy(dtype=float)
    e = expected.to_numpy(dtype=float)
    ok = np.isclose(a, e, rtol=rtol, atol=atol, equal_nan=True)
    if not ok.all():
        i, j = np.argwhere(~ok)[0]
        raise AssertionError(
            f"{label} {int((~ok).sum())} of {ok.size} cells differ; first at "
            f"[{expected.index[i]!r}, {expected.columns[j]!r}]: "
            f"python={a[i, j]!r} R={e[i, j]!r}"
        )


def lower_to_full(m: pd.DataFrame) -> pd.DataFrame:
    """Rebuild the full symmetric matrix from a lower-triangle one (``""`` blanks)."""
    a = m.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    iu = np.triu_indices_from(a, k=1)
    a[iu] = a.T[iu]
    return pd.DataFrame(a, index=m.index, columns=m.columns)


def known_gap(request, gaps: dict, key: str):
    """Mark the running test as an expected failure when ``key`` is a known gap.

    The marker is *strict*: when the gap is closed the test fails with
    "XPASS(strict)", which reminds us to delete the entry from ``gaps``.
    """
    if key in gaps:
        request.applymarker(pytest.mark.xfail(reason=gaps[key], strict=True))
