"""``log_coin``: the pieces and the behaviour that need no R reference."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from netcoin import log_coin
from netcoin.logcoin import _CASCADE, _comb, _expand, _ipf, _loglm, deparse_formula


@pytest.fixture
def survey():
    rng = np.random.default_rng(4)
    n = 400
    a = rng.choice(["x", "y", "z"], n)
    b = np.where((a == "x") & (rng.random(n) < 0.6), "p", rng.choice(["p", "q"], n))
    c = rng.choice(["u", "v"], n)
    return pd.DataFrame({"A": a, "B": b, "C": c, "w": rng.uniform(0.5, 2, n)})


# ---------------------------------------------------------------- formulas ----

def test_comb_makes_the_products_of_the_variables():
    assert _comb(["a", "b", "c"], None, 2) == [["a", "b"], ["a", "c"], ["b", "c"]]
    assert _comb(["a", "b", "c"], ["g"], 2) == [["a", "b"], ["a", "c"], ["b", "c"], ["a", "g"], ["b", "g"], ["c", "g"]]
    assert _comb(["a", "b", "c"], ["g"], 3) == [["a", "b", "c"], ["a", "b", "g"], ["a", "c", "g"], ["b", "c", "g"]]
    assert _comb(["a", "b"], ["g"], 3) == [["a", "g"], ["b", "g"]] or _comb(["a", "b"], ["g"], 3) == [["a", "b", "g"]]
    assert _comb(["a"], ["g"], 2) == [["a", "g"]]


def test_the_terms_of_a_product_follow_rs_expansion():
    assert _expand([["a", "b", "c"]], ["a", "b", "c"]) == [("a",), ("b",), ("c",), ("a", "b"), ("a", "c"),
                                                              ("b", "c"), ("a", "b", "c")]


def test_deparse_breaks_long_formulas_as_r_does():
    short = deparse_formula("weight", [["a", "b"], ["a", "c"]])
    assert short == "weight ~ a * b + a * c"
    long = deparse_formula("weight", [[f"variable{i}", f"variable{j}"] for i in range(1, 5) for j in range(i + 1, 6)])
    assert long == ("weight ~ variable1 * variable2 + variable1 * variable3 + variable1 *      variable4 + "
                    "variable1 * variable5 + variable2 * variable3 +      variable2 * variable4 + variable2 * "
                    "variable5 + variable3 *      variable4 + variable3 * variable5 + variable4 * variable5")


def test_the_cascade_of_contrasts_grows_with_the_pivots():
    assert [len(_CASCADE[k]) for k in range(1, 6)] == [1, 3, 7, 15, 31]


# ------------------------------------------------------- independent checks ----

def test_ipf_of_the_independence_model_is_the_expected_table():
    table = np.array([[10.0, 20, 5], [30, 5, 25]])
    fit = _ipf(table, [(0,), (1,)])
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / table.sum()
    np.testing.assert_allclose(fit, expected)
    lrt, df = _loglm(table, ["r", "c"], [("r",), ("c",), ("r", "c")][:2])
    g2 = stats.chi2_contingency(table, lambda_="log-likelihood", correction=False)
    assert lrt == pytest.approx(g2.statistic) and df == g2.dof


def test_the_saturated_model_fits_exactly():
    table = np.random.default_rng(0).integers(1, 30, size=(2, 3, 2)).astype(float)
    lrt, df = _loglm(table, ["a", "b", "c"], [("a",), ("b",), ("c",), ("a", "b"), ("a", "c"), ("b", "c"), ("a", "b", "c")])
    assert lrt == pytest.approx(0, abs=1e-9) and df == 0


def test_two_variables_with_order_two_is_the_saturated_model(survey):
    note = log_coin(survey, ["A", "B"]).options["note"]
    assert "L2= 0.00; d.f.= 0" in note and "Covariance structures: 6 (6)" in note and "weight ~ A * B" in note


def test_the_fit_of_three_variables_matches_the_poisson_deviance(survey):
    smf = pytest.importorskip("statsmodels.formula.api")
    sm = pytest.importorskip("statsmodels.api")
    counts = survey.groupby(["A", "B", "C"]).size().rename("n").reset_index()
    exact = smf.glm("n ~ A*B + A*C + B*C", data=counts, family=sm.families.Poisson()).fit()
    note = log_coin(survey, ["A", "B", "C"]).options["note"]
    reported = float(note.split("L2= ")[1].split(";")[0])
    assert reported == pytest.approx(exact.deviance, abs=0.1)       # IPF stops when the margins differ by < 0.1
    assert f"d.f.= {int(exact.df_resid)}" in note


def test_main_effects_are_deviations_that_add_up_to_zero(survey):
    net = log_coin(survey, ["A", "B", "C"], pmax=0.5)
    nodes = net.nodes.set_index("name")["Estimate"]
    for variable in ("A", "B", "C"):
        assert nodes[[n for n in nodes.index if n.startswith(variable)]].sum() == pytest.approx(0, abs=2e-3)


def test_with_two_categories_the_estimates_are_opposite(survey):
    nodes = log_coin(survey, ["A", "B", "C"]).nodes.set_index("name")["Estimate"]
    assert nodes["Bp"] == pytest.approx(-nodes["Bq"], abs=1e-3)


# ------------------------------------------------------------------- inputs ----

def test_spaces_in_names_become_underscores(survey):
    survey = survey.rename(columns={"A": "A a"})
    net = log_coin(survey, ["A a", "B"])
    assert "A_ax" in net.nodes["name"].tolist() and "A_a:x" in net.nodes["label"].tolist()


def test_weights_scale_the_counts(survey):
    unit = log_coin(survey, ["A", "B"], frequency=True)
    weighted = log_coin(survey, ["A", "B", "w"], weight="w", frequency=True)
    assert weighted.nodes["freq"].sum() != unit.nodes["freq"].sum()
    vector = log_coin(survey, ["A", "B"], weight=survey["w"].to_numpy(), frequency=True)
    pd.testing.assert_frame_equal(vector.nodes, weighted.nodes)


def test_weights_need_one_value_per_case(survey):
    with pytest.raises(ValueError, match="Weights"):
        log_coin(survey, ["A", "B"], weight=[1.0, 2.0])


def test_more_than_five_pivots_are_refused():
    rng = np.random.default_rng(1)
    frame = pd.DataFrame({f"v{i}": rng.choice(["a", "b"], 50) for i in range(6)})
    with pytest.raises(ValueError, match="more than 5 variables"):
        log_coin(frame)
    log_coin(frame, no_first_cat="v0")           # five pivots are fine


def test_unknown_variables_are_an_error(survey):
    with pytest.raises(ValueError, match="not in data: nothing"):
        log_coin(survey, ["A", "nothing"])


def test_missing_values_leave_the_case_out_of_the_table(survey):
    survey.loc[:19, "A"] = np.nan
    net = log_coin(survey, ["A", "B"], frequency=True)
    assert net.nodes.loc[net.nodes["var"] == "A", "freq"].sum() == 380 and net.options["scenarios"] == 400


def test_user_size_is_respected_and_defaults_follow_the_options(survey):
    assert log_coin(survey, ["A", "B"]).options["nodeSize"] == "Estimate"
    assert log_coin(survey, ["A", "B"], frequency=True).options["nodeSize"] == "freq"
    assert log_coin(survey, ["A", "B"], frequency=True, percentage=True).options["nodeSize"] == "prop"
    assert log_coin(survey, ["A", "B"], size="cat").options["nodeSize"] == "cat"


def test_one_sided_links_are_positive_and_two_sided_ones_are_bipolar(survey):
    one = log_coin(survey, ["A", "B", "C"], pmax=0.3)
    two = log_coin(survey, ["A", "B", "C"], pmax=0.3, twotails=True)
    assert (one.links["Estimate"] > 0).all() and "Pr(>z)" in one.links
    assert two.options["linkBipolar"] is True and two.options["linkColor"] == "Estimate"
    assert (two.links["Estimate"] < 0).any() and "Pr(>|z|)" in two.links


def test_links_go_from_the_later_variable_to_the_earlier(survey):
    net = log_coin(survey, ["A", "B", "C"], pmax=0.5)
    first = net.nodes.set_index("name")["var"]
    order = {"A": 0, "B": 1, "C": 2}
    assert all(order[first[s]] > order[first[t]] for s, t in zip(net.links["Source"], net.links["Target"]))


def test_higher_order_interactions_go_through_extra_nodes(survey):
    net = log_coin(survey, ["A", "B", "C"], order=3, pmax=0.5)
    extra = net.nodes[net.nodes["var"].str.startswith(".I")]
    assert len(extra) > 0 and extra["name"].str.match(r"u\d+").all()
    assert set(net.links["Source"]) & set(extra["name"])


def test_igraph_is_returned_when_asked(survey):
    igraph = pytest.importorskip("igraph")
    assert isinstance(log_coin(survey, ["A", "B"], igraph=True), igraph.Graph)


def test_spanish_names_the_columns(survey):
    net = log_coin(survey, ["A", "B"], language="es")
    assert {"nombre", "etiqueta"} <= set(net.nodes.columns) and net.options["nodeLabel"] == "etiqueta"
