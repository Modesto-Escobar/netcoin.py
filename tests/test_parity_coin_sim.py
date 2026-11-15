"""``coin`` and ``sim`` against netCoin 2.1.49."""
import pandas as pd
import pytest

from netcoin import coin, sim
from conftest import (assert_matrix_close, known_gap, lower_to_full, read_input, read_matrix)

# ----------------------------------------------------------------- coin ----

COIN_CASES = {
    # case: (input file, kwargs for the python coin())
    "galapagos": ("galapagos", {}),
    "dice": ("dice", {}),
    "galapagos_min3": ("galapagos", {"minimum": 3}),
    "galapagos_sort": ("galapagos", {"sort": True}),
    "galapagos_sort_inc": ("galapagos", {"sort": True, "decreasing": False}),
    "galapagos_weight": ("galapagos", {"weight": "weights"}),
    "galapagos_total": ("galapagos", {"total": True}),
    "galapagos_na_pairwise": ("galapagos_na", {"pairwise": True}),
    "galapagos_na_listwise": ("galapagos_na", {}),
}


COIN_GAPS = {
    "galapagos_min3": "coin() has no minimum/maximum arguments (phase 2)",
    "galapagos_sort": "coin() has no sort/decreasing arguments (phase 2)",
    "galapagos_sort_inc": "coin() has no sort/decreasing arguments (phase 2)",
    "galapagos_weight": "coin() has no weight argument (phase 2)",
    "galapagos_total": "coin() has no total argument (phase 2)",
    "galapagos_na_pairwise": "coin() has no pairwise argument (phase 2)",
    "galapagos_na_listwise": "R drops incomplete cases (na.omit); python turns NA into 0",
}


@pytest.mark.parametrize("case", COIN_CASES)
def test_coin_matrix(case, request):
    known_gap(request, COIN_GAPS, case)
    source, kwargs = COIN_CASES[case]
    kwargs = dict(kwargs)
    if kwargs.get("weight") == "weights":
        kwargs["weight"] = read_input("weights")["weight"].tolist()
    result = coin(read_input(source), print_n=False, **kwargs)
    assert_matrix_close(lower_to_full(result), read_matrix(f"coin_{case}"), label=case)


# ------------------------------------------------------------------ sim ----

# R similarity -> (python procedure key, names of the matrices R returns)
SIM_PROCS = {
    "matching": "m", "Rogers": "t", "Gower": "g", "Sneath": "s", "Anderberg": "and",
    "Jaccard": "j", "dice": "d", "antiDice": "a", "Ochiai": "o", "Kulczynski": "k",
    "Hamann": "ham", "Yule": "y", "Pearson": "p", "tetrachoric": None, "odds": "od",
    "Russell": "r", "expected": "e", "confidence": "con", "Haberman": "h", "Z": "z",
    "Fisher": "hyp", "coincidences": "cc", "relative": "x", "sConditional": None,
    "tConditional": None, "c.conditional": "i", "c.probable": "cp",
}
# R file names (R object names, '.' kept) produced by each procedure
SIM_FILES = {"confidence": ["conf.L", "confidence", "conf.U"]}
DATASETS = {"galapagos": "galapagos", "dice": "dice"}


SIM_GAPS = {
    "tetrachoric": "procedure not implemented in python",
    "sConditional": "procedure not implemented in python (python 'i' is c.conditional)",
    "tConditional": "procedure not implemented in python",
}
# R bug (netCoin 2.1.49): ``sqrt(N)[is.na(x)]`` fills only the first NA of
# Haberman with sqrt(N) and leaves NaN in the rest; python fills all of them.
R_HABERMAN_BUG = "R leaves NaN where python returns sqrt(N): sqrt(N)[is.na(x)] in sim()"
# a coin matrix built by python does not carry n, so it is guessed from the diagonal
NO_N = "python coin matrix has no n attribute: n is guessed as max(diag)"


def _py_sim(data, key, **kwargs):
    if key is None:
        raise NotImplementedError
    return sim(data, key, **kwargs)


def _as_dict(result, name):
    if isinstance(result, dict):
        return result
    return {name: result}


# the python keys of matrices returned by sim()
_PY_RESULT_KEY = {
    "matching": "matching", "Rogers": "rogers", "Gower": "gower", "Sneath": "sneath",
    "Anderberg": "anderberg", "Jaccard": "jaccard", "dice": "dice", "antiDice": "antidice",
    "Ochiai": "ochiai", "Kulczynski": "kulczynski", "Hamann": "hamann", "Yule": "yule",
    "Pearson": "pearson", "odds": "odds", "Russell": "russell", "expected": "expected",
    "Haberman": "haberman", "Z": "z", "Fisher": "fisher", "coincidences": "coincidences",
    "relative": "relative", "c.conditional": "sconditional", "c.probable": "tconditional",
}
_PY_CONF_KEYS = {"conf.L": "conf_l", "confidence": "confidence", "conf.U": "conf_u"}


def _sim_check(ds, proc, py_input, **kwargs):
    try:
        result = _py_sim(py_input, SIM_PROCS[proc], **kwargs)
    except NotImplementedError:
        pytest.fail("procedure not implemented in python")
    names = SIM_FILES.get(proc, [proc])
    for name in names:
        f = "sim_%s_%s" % (ds, name.replace("-", "_"))
        key = _PY_CONF_KEYS[name] if proc == "confidence" else _PY_RESULT_KEY[proc]
        got = _as_dict(result, key)[key]
        assert_matrix_close(got, read_matrix(f), rtol=1e-8, atol=1e-10, label=f)


# these do not depend on n, so a matrix without n is enough
INDEPENDENT_OF_N = {"Jaccard", "dice", "antiDice", "Ochiai", "Kulczynski", "coincidences",
                    "c.probable"}


def _sim_gaps(ds, proc, from_coin):
    if proc in SIM_GAPS:
        return SIM_GAPS[proc]
    if proc == "Haberman" and ds == "galapagos":
        return R_HABERMAN_BUG
    if from_coin and ds == "dice" and proc not in INDEPENDENT_OF_N:
        return NO_N
    return None


@pytest.mark.parametrize("proc", SIM_PROCS)
@pytest.mark.parametrize("ds", DATASETS)
def test_sim_from_incidences(ds, proc, request):
    reason = _sim_gaps(ds, proc, from_coin=False)
    known_gap(request, {"x": reason} if reason else {}, "x")
    _sim_check(ds, proc, read_input(DATASETS[ds]))


@pytest.mark.parametrize("proc", SIM_PROCS)
@pytest.mark.parametrize("ds", DATASETS)
def test_sim_from_coin_matrix(ds, proc, request):
    reason = _sim_gaps(ds, proc, from_coin=True)
    known_gap(request, {"x": reason} if reason else {}, "x")
    cm = coin(read_input(DATASETS[ds]), print_n=False)
    _sim_check(ds, proc, cm)


def test_sim_confidence_level_90():
    result = sim(read_input("galapagos"), "con", level=0.90)
    assert_matrix_close(result["conf_l"], read_matrix("sim_galapagos_conf.L_l90"), 1e-8, 1e-10)
    assert_matrix_close(result["confidence"], read_matrix("sim_galapagos_confidence_l90"), 1e-8, 1e-10)
    assert_matrix_close(result["conf_u"], read_matrix("sim_galapagos_conf.U_l90"), 1e-8, 1e-10)


def test_sim_minimum_on_incidences():
    result = sim(read_input("galapagos"), "j", minimum=3)
    assert_matrix_close(result, read_matrix("sim_galapagos_Jaccard_min3"), 1e-8, 1e-10)


def test_sim_haberman_matches_r_wherever_r_is_defined():
    """Haberman: identical to R except where R's bug leaves NaN (python gives sqrt(N))."""
    import numpy as np
    got = sim(read_input("galapagos"), "h")
    ref = read_matrix("sim_galapagos_Haberman")
    defined = ~ref.isna().to_numpy()
    assert np.allclose(got.to_numpy()[defined], ref.to_numpy()[defined], rtol=1e-9)
    assert np.allclose(got.to_numpy()[~defined], np.sqrt(17))
