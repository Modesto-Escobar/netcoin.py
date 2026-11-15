"""``glm_coin``: parsing, the marginal effects checked by hand, and what the object keeps."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

pytest.importorskip("statsmodels")
import statsmodels.api as sm
import statsmodels.formula.api as smf

from netcoin import glm_coin
from netcoin.glmcoin import (avg_slope, clean_variables, contrasts_gw, fit_model, hc1, parse_specs,
                             split_top_level, summary_nodes_links)


@pytest.fixture
def data():
    rng = np.random.default_rng(21)
    n = 400
    g = rng.choice(["f", "m"], n)
    a = rng.choice(["young", "mid", "old"], n, p=[.3, .4, .3])
    x = rng.normal(0, 1, n)
    y = 1 + 0.8 * (g == "m") + 0.5 * x + 0.6 * (a == "old") + rng.normal(0, 1, n)
    flag = (rng.random(n) < 1 / (1 + np.exp(-(0.2 * y - 1)))).astype(int)
    return pd.DataFrame({"g": g, "age": a, "x": x, "y": y, "flag": flag, "w": rng.uniform(.5, 2, n),
                         "cnt": rng.poisson(np.exp(0.3 + 0.2 * x))})


# ------------------------------------------------------------------- parsing ----

def test_the_top_level_comma_separates_the_family():
    assert split_top_level("y ~ x, binomial") == ("y ~ x", "binomial")
    assert split_top_level("y ~ log(x, 2), poisson") == ("y ~ log(x, 2)", "poisson")
    assert split_top_level("y ~ x") == ("y ~ x", None)


def test_equations_are_read_line_by_line():
    specs = parse_specs("# a comment\ny ~ x\n\nz ~ y + x, binomial\nw ~ x, family = Gamma()")
    assert [(f, fam) for _, f, fam in specs] == [("y ~ x", None), ("z ~ y + x", "binomial"), ("w ~ x", "Gamma")]
    with pytest.raises(ValueError, match="Family not allowed: 'negbin'"):
        parse_specs("y ~ x, negbin")
    with pytest.raises(ValueError, match="formula"):
        parse_specs("just words")


def test_names_are_made_safe_like_r_does():
    assert clean_variables(["Social participation", "edad (años)", "ok_1"]) == ["Social.participation", "edad..anos.", "ok_1"]


# --------------------------------------------------- the effects, by hand ----

def test_the_effect_of_a_factor_is_the_average_prediction_minus_the_weighted_mean(data):
    model = fit_model("y ~ g + age", "y ~ g + age", None, data, None, False)
    table = contrasts_gw(model, ["g"])
    ols = smf.ols("y ~ g + age", data).fit()
    averages = []
    for level in ("f", "m"):
        counterfactual = data.assign(g=level)
        averages.append(ols.predict(counterfactual).mean())
    shares = data["g"].value_counts(normalize=True).reindex(["f", "m"]).to_numpy()
    want = np.array(averages) - (np.array(averages) * shares).sum()
    np.testing.assert_allclose(table["estimate"], want, rtol=1e-9)
    assert table["Source"].tolist() == ["g:f", "g:m"] and (table["Target"] == "y").all()
    # for a linear model the contrast of two categories is their coefficient
    assert table["estimate.base"].iloc[1] == pytest.approx(ols.params["g[T.m]"], rel=1e-9)
    assert table["std.error.base"].iloc[1] == pytest.approx(ols.bse["g[T.m]"], rel=1e-6)
    assert np.isnan(table["std.error.base"].iloc[0]) and table["estimate.base"].iloc[0] == 0


def test_the_slope_of_a_numeric_variable_in_a_linear_model_is_its_coefficient(data):
    model = fit_model("y ~ g + x", "y ~ g + x", None, data, None, False)
    slope = avg_slope(model, "x")
    ols = smf.ols("y ~ g + x", data).fit()
    assert slope["estimate"].iloc[0] == pytest.approx(ols.params["x"], rel=1e-7)
    assert slope["std.error"].iloc[0] == pytest.approx(ols.bse["x"], rel=1e-5)
    assert slope["p.value"].iloc[0] == pytest.approx(2 * stats.norm.sf(abs(ols.params["x"] / ols.bse["x"])), rel=1e-4)


def test_in_a_logit_the_effects_are_percentage_points_of_the_probability(data):
    net = glm_coin("flag ~ g + x, binomial", data, pmax=1)
    model = fit_model("flag ~ g + x", "flag ~ g + x", "binomial", data, None, False)
    table = contrasts_gw(model, ["g"])
    manual = [model.sm_result.predict(sm.add_constant(pd.DataFrame({"g": [1.0 if l == "m" else 0.0] * len(data), "x": data["x"]}),
                                                      has_constant="add").to_numpy()).mean() * 100 for l in ("f", "m")]
    # scaled to percent (a quasibinomial is fitted)
    assert table["estimate.base"].iloc[1] == pytest.approx(manual[1] - manual[0], rel=1e-6)
    assert (net.links["Target"] == "flag:1").all()


def test_the_covariance_can_be_robust(data):
    plain = glm_coin("y ~ g + age", data, pmax=1)
    robust = glm_coin("y ~ g + age", data, pmax=1, robust=True)
    assert plain.links["estimate"].tolist() == pytest.approx(robust.links["estimate"].tolist())
    assert plain.links["std.error"].tolist() != robust.links["std.error"].tolist()
    ols = smf.ols("y ~ g + age", data).fit(cov_type="HC1")
    model = fit_model("y ~ g + age", "y ~ g + age", None, data, None, True)
    np.testing.assert_allclose(np.sqrt(np.diag(model.cov)), np.asarray(ols.bse), rtol=1e-6)


def test_hc1_of_a_logit_matches_statsmodels(data):
    model = fit_model("flag ~ g + x", "flag ~ g + x", "binomial", data, None, False)
    X = model.X
    fit = sm.GLM(data["flag"].to_numpy(float), X, family=sm.families.Binomial()).fit(cov_type="HC1")
    n, k = X.shape                          # sandwich (and R) scale HC0 by n / (n - k); statsmodels does not
    np.testing.assert_allclose(np.sqrt(np.diag(hc1("binomial", X, model.y, model.coef, np.ones(len(X))))),
                               np.asarray(fit.bse) * np.sqrt(n / (n - k)), rtol=1e-5)


def test_weights_count_in_the_nodes(data):
    net = glm_coin("y ~ g + x", data, weights="w", pmax=1)
    freq = net.nodes.set_index("name")["Frequency"]
    assert freq["g:m"] == pytest.approx(data.loc[data["g"] == "m", "w"].sum())
    assert freq["x"] == pytest.approx(data["w"].sum())
    mean = net.nodes.set_index("name")["Mean"]["x"]
    assert mean == pytest.approx((data["x"] * data["w"]).sum() / data["w"].sum())


def test_missing_values_leave_the_case_out_of_the_equation(data):
    data.loc[:9, "x"] = np.nan
    net = glm_coin("y ~ g + x", data, pmax=1)
    assert net.equations["n_obs"].iloc[0] == 390 and net.nodes.set_index("name")["Frequency"]["x"] == 390


def test_a_factor_with_one_level_gives_no_links(data):
    data = data.assign(k="same")
    net = glm_coin("y ~ k + x", data, pmax=1)
    assert "k:same" not in net.nodes["name"].tolist() or not (net.links["Source"] == "k:same").any()


def test_names_with_spaces_come_back_in_the_nodes_and_equations(data):
    data = data.rename(columns={"g": "the group", "age": "age group"})
    net = glm_coin("y ~ the group + age group", data, pmax=1)
    assert {"the group:f", "age group:old"} <= set(net.nodes["name"])
    assert set(net.nodes["variable"]) >= {"the group", "age group"}
    assert net.equations["model"].iloc[0] == "y ~ the group + age group"
    assert net.coeffs[0]["Source"].tolist()[1].startswith("the group:") and net.links["Model"].iloc[0] == "y ~ the group + age group"


def test_several_equations_stack_their_links_and_share_the_nodes(data):
    net = glm_coin("y ~ g + x\nflag ~ y + g, binomial", data, pmax=1)
    assert net.links["Model"].nunique() == 2 and net.nodes["name"].is_unique
    assert {"y", "flag:1"} <= set(net.nodes["name"]) and len(net.equations) == 2 and len(net.coeffs) == 2


def test_link_defaults_are_those_of_glm_coin(data):
    net = glm_coin("y ~ g + x", data)
    assert net.options["showArrows"] is True and net.options["linkBipolar"] is True
    assert net.options["linkWidth"] == "s.value" and net.options["linkColor"] == "estimate"
    assert net.options["nodeColor"] == "variable"
    hidden = net.links["hidden"] if "hidden" in net.links else pd.Series(False, index=net.links.index)
    assert hidden.tolist() == (net.links["p.value"] >= 0.05).tolist()


def test_multinomial_matches_statsmodels_and_gives_one_group_per_class(data):
    from netcoin.glmcoin import fit_multinomial, _softmax_probs
    model = fit_model("age ~ x + g", "age ~ x + g", "multinomial", data, None, False)
    assert model.classes == ["mid", "old", "young"]
    reference = sm.MNLogit(model.y, pd.DataFrame(model.X)).fit(disp=0)
    ours = model.coef.reshape(2, 3)
    np.testing.assert_allclose(ours.T, np.asarray(reference.params), rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(np.sqrt(np.diag(model.cov)), np.asarray(reference.bse).T.ravel(), rtol=1e-5)
    probs = _softmax_probs(model.X, model.coef, 3)
    np.testing.assert_allclose(probs.sum(axis=1), 1.0)
    table = contrasts_gw(model, ["g"])
    assert sorted(set(table["Target"])) == ["age:mid", "age:old", "age:young"] and len(table) == 6


def test_a_multinomial_network_has_a_node_and_a_link_per_class(data):
    net = glm_coin("age ~ g + x, multinomial", data, pmax=1)
    assert {"age:young", "age:mid", "age:old"} <= set(net.nodes["name"])
    assert set(net.links["Target"]) == {"age:young", "age:mid", "age:old"}
    assert net.equations["family"].iloc[0] == "multinomial" and net.equations["df.residual"].isna().iloc[0]
    assert set(net.coeffs[0]["Target"]) == {"age:old", "age:young"}                 # the first class is the reference


def test_the_weights_of_a_multinomial_count_even_with_missing_values(data):
    data.loc[:9, "x"] = np.nan
    net = glm_coin("age ~ g + x, multinomial", data, weights="w", pmax=1)
    kept = data.dropna(subset=["x"])
    assert net.nodes.set_index("name")["Frequency"]["g:m"] == pytest.approx(kept.loc[kept["g"] == "m", "w"].sum())


def test_the_response_of_a_multinomial_must_be_a_factor(data):
    with pytest.raises(ValueError, match="must be a factor"):
        glm_coin("cnt ~ x, multinomial", data)


def test_aliased_terms_are_left_out_of_the_fit_instead_of_failing(data):
    data = data.assign(x2=2 * data["x"])
    with pytest.warns(UserWarning, match="aliased"):
        net = glm_coin("y ~ g + x + x2", data, pmax=1)
    row = net.links[net.links["Source"] == "x2"].iloc[0]
    assert np.isnan(row["estimate"]) and "x2" not in net.coeffs[0]["term"].tolist()
    assert net.equations["df.residual"].iloc[0] == 400 - 3
    assert net.links.loc[net.links["Source"] == "x", "estimate"].notna().all()


def test_unknown_data_problems_are_reported(data):
    with pytest.raises(ValueError, match="not found in data"):
        glm_coin("y ~ x", data, weights="nope")
    with pytest.raises(ValueError, match="'weights' does not match"):
        glm_coin("y ~ x", data, weights=[1.0, 2.0])


def test_the_coefficient_tables_and_fit_are_kept(data):
    net = glm_coin("y ~ g + x", data, pmax=1)
    table = net.coeffs[0]
    assert list(table["term"]) == ["(Intercept)", "gm", "x"] and list(table["Source"]) == ["(Intercept)", "g:m", "x"]
    ols = smf.ols("y ~ g + x", data).fit()
    np.testing.assert_allclose(table["estimate"], ols.params.to_numpy(), rtol=1e-9)
    np.testing.assert_allclose(table["p.value"], ols.pvalues.to_numpy(), rtol=1e-6)
    assert table.attrs["family"] == "gaussian" and table.attrs["sd.y"] == pytest.approx(data["y"].std(ddof=0))
    row = net.equations.iloc[0]
    assert row["r2"] == pytest.approx(ols.rsquared, abs=6e-4) and row["aic"] == pytest.approx(ols.aic + 2, abs=0.011)   # R counts the variance as a parameter
    assert row["logLik"] == pytest.approx(ols.llf, abs=0.011) and row["df.residual"] == ols.df_resid


def test_igraph_returns_a_graph_with_the_fit(data):
    igraph = pytest.importorskip("igraph")
    graph = glm_coin("y ~ g + x", data, igraph=True)
    assert isinstance(graph, igraph.Graph) and "equations" in graph.attributes()


def test_the_page_can_be_written(data, tmp_path):
    glm_coin("y ~ g + x", data, directory=str(tmp_path / "glm"))
    assert (tmp_path / "glm" / "index.html").is_file()
