"""``glm_coin`` against netCoin 2.1.49 (with its own marginal effects)."""
import json
import warnings

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("statsmodels")
from netcoin import glm_coin
from conftest import EXPECTED, compare_options, read_expected, read_input

GD = read_input("glm_data", strings=True)
GD["cls"] = pd.Categorical(GD["cls"], categories=["low", "mid", "high"])   # R's factor levels, in this order
for _c in ("cweight", "income", "edu", "cnt", "pos", "flag"):
    GD[_c] = GD[_c].astype(float)

# python call of each R call (the R calls are in tests/reference/make_reference.R)
CASES = {
    "lm": lambda: glm_coin("income ~ Gender + Age", GD),
    "lm_inter": lambda: glm_coin("income ~ Gender*Age + edu", GD),
    "lm_numint": lambda: glm_coin("income ~ Gender*edu", GD),
    "logit": lambda: glm_coin("Social participation ~ Gender + income, binomial", GD, pmax=0.5),
    "logit01": lambda: glm_coin("flag ~ Gender + edu + Age, binomial", GD, pmax=0.5),
    "logit_inter": lambda: glm_coin("Political participation ~ Gender*Social participation + edu, binomial", GD, pmax=0.5),
    "poisson": lambda: glm_coin("cnt ~ Gender + income, poisson", GD, pmax=0.5),
    "quasipoisson": lambda: glm_coin("cnt ~ Gender + Age, quasipoisson", GD, pmax=0.5),
    "gamma": lambda: glm_coin("pos ~ Age + edu, Gamma", GD, pmax=0.5),
    "weights": lambda: glm_coin("income ~ Gender + Age", GD, weights="cweight"),
    "weights_logit": lambda: glm_coin("Social participation ~ Gender + edu, binomial", GD, weights="cweight", pmax=0.5),
    "robust": lambda: glm_coin("income ~ Gender + Age", GD, robust=True),
    "robust_logit": lambda: glm_coin("Social participation ~ Gender + income, binomial", GD, robust=True, pmax=0.5),
    "chain": lambda: glm_coin("income ~ Gender + Age\nSocial participation ~ income + Gender, binomial", GD, pmax=0.5),
    "spanish": lambda: glm_coin("income ~ Gender + Age", GD, language="es", pmax=0.5),
    "mn": lambda: glm_coin("cls ~ Gender + edu, multinomial", GD, pmax=0.9),
    "mn_weights": lambda: glm_coin("cls ~ Gender + edu, multinomial", GD, weights="cweight", pmax=0.9),
    "mn_inter": lambda: glm_coin("cls ~ Gender*Social participation + edu, multinomial", GD, pmax=0.9),
    "mn_binary": lambda: glm_coin("Social participation ~ Gender + edu, multinomial", GD, pmax=0.9),
    "mn_chain": lambda: glm_coin("edu ~ Gender + Age\ncls ~ edu + Gender, multinomial", GD, pmax=0.9),
    "options": lambda: glm_coin("income ~ Gender + Age", GD, pmax=0.2, show_arrows=False, main="GLM"),
}
# R stops its iterations when the deviance changes by less than 1e-8, python iterates further
# nnet stops its optimiser at a relative change of 1e-8, which leaves the coefficients 1e-5 from the maximum
MN = ("mn", "mn_weights", "mn_inter", "mn_binary", "mn_chain")
RTOL = {"gamma": 3e-4, **{c: 1e-3 for c in MN}}
ATOL = {c: 2e-3 for c in MN}
# R's eqCoefs finds no weights in the model frame of a multinom, so its sd.x and latent.sd are not weighted;
# python's are (see reference/README.md)
UNWEIGHTED = {"mn_weights"}


def _run(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return CASES[case]()


def _compare(actual: pd.DataFrame, expected: pd.DataFrame, label: str, rtol=1e-5, atol=1e-8):
    assert list(actual.columns) == list(expected.columns), f"{label}: {list(actual.columns)}"
    assert len(actual) == len(expected), f"{label}: {len(actual)} rows, R has {len(expected)}"
    for col in expected.columns:
        a, e = actual[col].reset_index(drop=True), expected[col].reset_index(drop=True)
        if pd.api.types.is_numeric_dtype(e) and not pd.api.types.is_bool_dtype(e):
            np.testing.assert_allclose(a.astype(float), e.astype(float), rtol=rtol, atol=atol, equal_nan=True,
                                       err_msg=f"{label}/{col}")
        else:
            assert a.astype(str).tolist() == e.astype(str).tolist(), f"{label}/{col}"


@pytest.mark.parametrize("case", CASES)
def test_nodes_match_r(case):
    net = _run(case)
    nodes = net.nodes.drop(columns=["hidden"], errors="ignore")
    _compare(nodes, read_expected(f"net_gl_{case}_nodes", strings=False).drop(columns=["hidden"], errors="ignore"),
             f"{case}/nodes")


@pytest.mark.parametrize("case", CASES)
def test_links_match_r(case):
    net = _run(case)
    ref = json.loads((EXPECTED / f"net_gl_{case}.json").read_text(encoding="utf-8"))
    assert net.lcount() == ref["nlinks"]
    if ref["nlinks"]:
        _compare(net.links, read_expected(f"net_gl_{case}_links"), f"{case}/links", rtol=RTOL.get(case, 1e-5),
                 atol=ATOL.get(case, 1e-8))


@pytest.mark.parametrize("case", CASES)
def test_options_match_r(case):
    net = _run(case)
    ref = json.loads((EXPECTED / f"net_gl_{case}.json").read_text(encoding="utf-8"))
    compare_options(net.options, ref["options"], case)


@pytest.mark.parametrize("case", CASES)
def test_equations_match_r(case):
    net = _run(case)
    _compare(net.equations, read_expected(f"gl_{case}_equations"), f"{case}/equations", rtol=RTOL.get(case, 1e-6),
             atol=0.011)


@pytest.mark.parametrize("case", CASES)
def test_coefficients_match_r(case):
    net = _run(case)
    ref = json.loads((EXPECTED / f"gl_{case}_coeffs.json").read_text(encoding="utf-8"))
    assert len(net.coeffs) == len(ref)
    for table, (name, entry) in zip(net.coeffs, ref.items()):
        want = pd.DataFrame(entry["table"])
        want = want.where(want.notna(), np.nan)
        skip = ["sd.x"] if case in UNWEIGHTED else []
        _compare(table.drop(columns=skip), want.drop(columns=skip), f"{case}/{name}", rtol=RTOL.get(case, 1e-5),
                 atol=ATOL.get(case, 1e-8))
        assert table.attrs["family"] == entry["family"] and table.attrs["link"] == entry["link"]
        sd_y = entry["sd.y"]
        if sd_y is None:
            assert np.isnan(table.attrs["sd.y"])
        else:
            assert table.attrs["sd.y"] == pytest.approx(sd_y, rel=1e-6)
        for key, value in ({} if case in UNWEIGHTED else entry["latent.sd"]).items():
            got = table.attrs["latent.sd"][key]
            assert (np.isnan(got) if value is None else got == pytest.approx(value, rel=RTOL.get(case, 1e-6)))
