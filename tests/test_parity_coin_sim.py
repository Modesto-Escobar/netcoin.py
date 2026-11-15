"""``coin`` and ``sim`` against netCoin 2.1.49."""
import numpy as np
import pandas as pd
import pytest

from netcoin import Coin, coin, sim
from conftest import assert_matrix_close, known_gap, read_expected, read_input, read_matrix

# ----------------------------------------------------------------- coin ----

COIN_CASES = {
    # case: (input file, kwargs for coin())
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


def _coin(case):
    source, kwargs = COIN_CASES[case]
    kwargs = dict(kwargs)
    if kwargs.get("weight") == "weights":
        kwargs["weight"] = read_input("weights")["weight"].tolist()
    return coin(read_input(source), **kwargs)


def _drop_duplicate_total(m: pd.DataFrame) -> pd.DataFrame:
    """R (total = TRUE) lists "Total" twice; keep the first."""
    keep = [i for i in range(len(m)) if i != 1]
    return m.iloc[keep, keep]


@pytest.mark.parametrize("case", COIN_CASES)
def test_coin_matrix(case):
    result = _coin(case)
    expected = read_matrix(f"coin_{case}")
    if case == "galapagos_total":
        expected = _drop_duplicate_total(expected)
        expected.columns = expected.index = list(result.matrix.index)  # 'Total.1' -> 'Total'
    assert isinstance(result, Coin)
    assert_matrix_close(result.matrix, expected, label=case)


@pytest.mark.parametrize("case", COIN_CASES)
def test_coin_n(case):
    attrs = read_expected(f"coin_{case}_attrs")
    assert _coin(case).n == pytest.approx(float(attrs["n"][0]))


@pytest.mark.parametrize("case", ["galapagos_na_pairwise"])
def test_coin_pairwise_attributes(case):
    result = _coin(case)
    assert result.pairwise
    assert_matrix_close(result.m, read_matrix(f"coin_{case}_m"), label="m")
    assert_matrix_close(result.x, read_matrix(f"coin_{case}_x"), label="x")


def test_coin_is_symmetric_and_integer_when_unweighted():
    c = _coin("galapagos")
    assert (c.matrix.to_numpy() == c.matrix.to_numpy().T).all()
    assert all(np.issubdtype(t, np.integer) for t in c.matrix.dtypes)


# ------------------------------------------------------------------ sim ----

# R similarity -> the names of the matrices R returns for it
SIM_PROCS = {
    "matching": ["matching"], "Rogers": ["Rogers"], "Gower": ["Gower"], "Sneath": ["Sneath"],
    "Anderberg": ["Anderberg"], "Jaccard": ["Jaccard"], "dice": ["dice"],
    "antiDice": ["antiDice"], "Ochiai": ["Ochiai"], "Kulczynski": ["Kulczynski"],
    "Hamann": ["Hamann"], "Yule": ["Yule"], "Pearson": ["Pearson"],
    "tetrachoric": ["tetrachoric"], "odds": ["odds"], "Russell": ["Russell"],
    "expected": ["expected"], "confidence": ["conf.L", "confidence", "conf.U"],
    "Haberman": ["Haberman"], "Z": ["Z"], "Fisher": ["Fisher"],
    "coincidences": ["coincidences"], "relative": ["relative"],
    "sConditional": ["sConditional"], "tConditional": ["tConditional"],
    "c.conditional": ["c.conditional"], "c.probable": ["c.probable"],
}
DATASETS = {"galapagos": "galapagos", "dice": "dice"}

# R bug (netCoin 2.1.49): ``sqrt(N)[is.na(x)]`` fills only the first NA of
# Haberman with sqrt(N) and leaves NaN in the rest; python fills all of them.
R_HABERMAN_BUG = "R leaves NaN where python returns sqrt(N): sqrt(N)[is.na(x)] in sim()"


def _as_dict(result, name):
    return result if isinstance(result, dict) else {name: result}


def _check_sim(ds, proc, data, suffix="", rtol=1e-8, atol=1e-10, **kwargs):
    result = _as_dict(sim(data, proc, **kwargs), proc)
    for name in SIM_PROCS[proc]:
        ref = read_matrix(f"sim_{ds}_{name.replace('-', '_')}{suffix}")
        assert_matrix_close(result[name], ref, rtol, atol, label=f"{ds}/{name}{suffix}")


def _sim_gap(ds, proc):
    return {"Haberman": R_HABERMAN_BUG}.get(proc) if ds == "galapagos" else None


@pytest.mark.parametrize("proc", SIM_PROCS)
@pytest.mark.parametrize("ds", DATASETS)
def test_sim_from_incidences(ds, proc, request):
    known_gap(request, {"x": _sim_gap(ds, proc)} if _sim_gap(ds, proc) else {}, "x")
    _check_sim(ds, proc, read_input(DATASETS[ds]))


@pytest.mark.parametrize("proc", SIM_PROCS)
@pytest.mark.parametrize("ds", DATASETS)
def test_sim_from_coin(ds, proc, request):
    known_gap(request, {"x": _sim_gap(ds, proc)} if _sim_gap(ds, proc) else {}, "x")
    _check_sim(ds, proc, coin(read_input(DATASETS[ds])))


def test_sim_confidence_level_90():
    result = sim(read_input("galapagos"), "confidence", level=0.90)
    for name, suffix in [("conf.L", "_l90"), ("confidence", "_l90"), ("conf.U", "_l90")]:
        assert_matrix_close(result[name], read_matrix(f"sim_galapagos_{name}{suffix}"), 1e-8, 1e-10)


def test_sim_minimum_on_incidences():
    result = sim(read_input("galapagos"), "J", minimum=3)
    assert_matrix_close(result, read_matrix("sim_galapagos_Jaccard_min3"), 1e-8, 1e-10)


def test_sim_distance():
    """R returns a 'dist' object (1 - similarity, zero diagonal); python the full matrix."""
    result = sim(coin(read_input("galapagos")), "Jaccard", distance=True)
    assert_matrix_close(result, read_matrix("sim_galapagos_Jaccard_dist"), 1e-9, 1e-12)


def test_sim_haberman_matches_r_wherever_r_is_defined():
    """Haberman: identical to R except where R's bug leaves NaN (python gives sqrt(N))."""
    got = sim(read_input("galapagos"), "Haberman")
    ref = read_matrix("sim_galapagos_Haberman")
    defined = ~ref.isna().to_numpy()
    assert np.allclose(got.to_numpy()[defined], ref.to_numpy()[defined], rtol=1e-9)
    assert np.allclose(got.to_numpy()[~defined], np.sqrt(17))


# ------------------------------------------------- pairwise and weighted ----

PAIRWISE = ["Jaccard", "Haberman", "Z", "expected", "confidence", "Fisher", "c.conditional",
            "Pearson"]


@pytest.mark.parametrize("proc", PAIRWISE)
def test_sim_pairwise(proc):
    c = coin(read_input("galapagos_na"), pairwise=True)
    _check_sim("galapagos_na", proc, c, suffix="_pairwise", pairwise=True)


@pytest.mark.parametrize("proc", PAIRWISE)
def test_sim_pairwise_from_incidences(proc):
    _check_sim("galapagos_na", proc, read_input("galapagos_na"), suffix="_pairwise", pairwise=True)


WEIGHTED = ["Jaccard", "Haberman", "Z", "expected", "confidence", "Fisher"]


@pytest.mark.parametrize("proc", WEIGHTED)
def test_sim_weighted(proc, request):
    known_gap(request, {"Haberman": R_HABERMAN_BUG}, proc)
    weights = read_input("weights")["weight"].tolist()
    _check_sim("galapagos", proc, read_input("galapagos"), suffix="_weight", weight=weights)


# ----------------------------------------------------- procedure names ----

@pytest.mark.parametrize("name,code", [
    ("J", "Jaccard"), ("Ham", "Hamann"), ("CC", "c.conditional"), ("CP", "c.probable"),
    ("Hyp", "Fisher"), ("fis", "Fisher"), ("Con", "confidence"), ("Od", "odds"),
    ("Ro", "Rogers"), ("And", "Anderberg"), ("Tet", "tetrachoric"), ("II", "tConditional"),
    ("I", "sConditional"), ("Coi", "coincidences"), ("c.conditional", "c.conditional"),
    ("sConditional", "sConditional"), ("haberman", "Haberman"), ("Frequencies", "coincidences"),
])
def test_procedure_names_follow_r(name, code):
    c = coin(read_input("galapagos"))
    got, want = _as_dict(sim(c, name), "x"), _as_dict(sim(c, code), "x")
    assert got.keys() == want.keys()
    for key in want:
        assert got[key].equals(want[key])


def test_unknown_procedure():
    with pytest.raises(ValueError):
        sim(coin(read_input("galapagos")), "?")
