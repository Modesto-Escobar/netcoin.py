"""``sur_coin`` and ``cob_coin``: inputs, warnings and outputs that need no R reference."""
import warnings

import numpy as np
import pandas as pd
import pytest

from netcoin import Coin, cob_coin, sur_coin


@pytest.fixture
def survey():
    rng = np.random.default_rng(3)
    n = 300
    sex = rng.choice(["Man", "Woman"], n)
    age = rng.choice(["young", "middle", "old"], n)
    works = np.where((age == "middle") & (rng.random(n) < 0.8), "yes", np.where(rng.random(n) < 0.3, "yes", "no"))
    income = 20 + 5 * (works == "yes") + rng.normal(0, 3, n)
    return pd.DataFrame({"sex": sex, "age": age, "works": works, "income": income,
                         "w": rng.uniform(0.5, 1.5, n)})


def test_every_category_is_a_node_named_after_its_variable(survey):
    net = sur_coin(survey, ["sex", "age"])
    assert net.nodes["name"].tolist() == ["sex:Man", "sex:Woman", "age:middle", "age:old", "age:young"]
    assert net.nodes["variable"].tolist() == ["sex"] * 2 + ["age"] * 3
    assert net.options["nodeSize"] == "%" and net.options["scenarios"] == 300


def test_missing_values_are_cases_without_a_category(survey):
    survey.loc[:9, "age"] = np.nan
    net = sur_coin(survey, ["sex", "age"], frequency=True, percentages=False)
    assert net.options["scenarios"] == 300
    assert net.nodes["frequency"][net.nodes["name"].str.startswith("age")].sum() == 290


def test_a_constant_dichotomy_is_ignored(survey):
    survey["always"] = "yes"
    net = sur_coin(survey, ["sex", "age"], dichotomies=["always", "works"], value_dicho="yes")
    assert "always" not in net.nodes["name"].tolist() and "works" in net.nodes["name"].tolist()


def test_one_value_of_dicho_for_each_dichotomy(survey):
    survey["other"] = np.where(survey["sex"] == "Man", "a", "b")
    net = sur_coin(survey, ["age"], dichotomies=["works", "other"], value_dicho=["yes", "b"],
                   frequency=True, percentages=False)
    got = dict(zip(net.nodes["name"], net.nodes["frequency"]))
    assert got["works"] == (survey["works"] == "yes").sum() and got["other"] == (survey["sex"] == "Woman").sum()


def test_numbers_can_be_the_value_of_a_dichotomy(survey):
    survey["flag"] = (survey["works"] == "yes").astype(int)
    net = sur_coin(survey, ["age"], dichotomies="flag", value_dicho=1, frequency=True, percentages=False)
    assert dict(zip(net.nodes["name"], net.nodes["frequency"]))["flag"] == (survey["flag"] == 1).sum()


def test_the_nodes_table_restricts_orders_and_adds_columns(survey):
    nodes = pd.DataFrame({"name": ["age:old", "sex:Man", "nope"], "group": ["x", "y", "z"]})
    with pytest.warns(UserWarning, match="nope"):
        net = sur_coin(survey, ["sex", "age"], nodes=nodes, color="group")
    assert sorted(net.nodes["name"]) == ["age:old", "sex:Man"]
    assert net.nodes.set_index("name")["group"].to_dict() == {"age:old": "x", "sex:Man": "y"}


def test_matrix_layout_follows_the_rows_of_the_nodes_not_their_final_order(survey):
    nodes = pd.DataFrame({"name": ["sex:Woman", "age:old", "sex:Man"]})
    layout = np.array([[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]])
    net = sur_coin(survey, ["sex", "age"], nodes=nodes, layout=layout, max_link=1)
    got = net.nodes.set_index("name")[["fx", "fy"]].to_dict("index")
    assert got["sex:Woman"] == {"fx": 0.0, "fy": 1.0} and got["sex:Man"] == {"fx": 4.0, "fy": 5.0}


def test_matrix_layout_needs_a_coordinate_per_node(survey):
    nodes = pd.DataFrame({"name": ["age:old", "sex:Man"]})
    with pytest.warns(UserWarning, match="coordinate per node"):
        sur_coin(survey, ["sex", "age"], nodes=nodes, layout=np.zeros((5, 2)), max_link=1)
    with pytest.warns(UserWarning, match="nodes variable"):
        sur_coin(survey, ["sex", "age"], layout=np.zeros((5, 2)), max_link=1)


def test_the_coordinates_of_mc_and_pc_go_with_the_name_of_each_node(survey):
    nodes = pd.DataFrame({"name": ["age:young", "sex:Woman", "age:old", "sex:Man", "age:middle"]})
    ordered = sur_coin(survey, ["sex", "age"], layout="mc", max_link=1).nodes.set_index("name")
    chosen = sur_coin(survey, ["sex", "age"], nodes=nodes, layout="mc", max_link=1).nodes.set_index("name")
    for name in chosen.index:
        assert chosen.loc[name, ["fx", "fy"]].tolist() == pytest.approx(ordered.loc[name, ["fx", "fy"]].tolist())


def test_one_qualitative_category_is_not_enough(survey):
    with pytest.raises(ValueError, match="net_corr"):
        sur_coin(survey, ["sex"], dichotomies=None, nodes=pd.DataFrame({"name": ["sex:Man"]}))


def test_unknown_variables_are_an_error(survey):
    with pytest.raises(ValueError, match="not in data: nothing"):
        sur_coin(survey, ["sex", "nothing"])


def test_the_weights_need_one_value_per_case(survey):
    with pytest.raises(ValueError, match="Weights"):
        sur_coin(survey, ["sex", "age"], weight=[1.0, 2.0])


def test_weights_are_kept_with_their_cases_when_cases_are_left_out(survey):
    """R forgets to leave out the weights of the incomplete cases; here they go with them."""
    survey.loc[[0, 5, 6], "age"] = np.nan
    by_vector = sur_coin(survey, ["sex", "age"], weight=survey["w"].to_numpy(), max_link=1)
    by_name = sur_coin(survey, ["sex", "age", "w"], weight="w", max_link=1)
    pd.testing.assert_frame_equal(by_vector.links, by_name.links)
    assert by_vector.options["scenarios"] == 297


def test_coin_returns_the_coincidences_of_the_categories(survey):
    result = sur_coin(survey, ["sex", "age"], coin=True)
    assert isinstance(result, Coin) and result.n == 300
    assert result.labels == ["sex:Man", "sex:Woman", "age:middle", "age:old", "age:young"]


def test_no_links_left_gives_a_network_without_links(survey):
    with pytest.warns(UserWarning):
        net = sur_coin(survey, ["sex", "age"], min_link=50)
    assert net.lcount() == 0 and net.links is None


def test_link_attributes_are_added_to_the_procedures(survey):
    net = sur_coin(survey, ["sex", "age"], lwidth="Jaccard", lcolor="Z", max_link=1)
    assert {"Haberman", "Jaccard", "p(Z)"} <= set(net.links.columns)
    assert net.options["linkWidth"] == "Jaccard" and net.options["linkColor"] == "p(Z)"


def test_exogenous_variables_hide_the_links_between_them(survey):
    net = sur_coin(survey, ["sex", "age", "works"], exogenous=["sex", "age"], max_link=1)
    chained = net.links["chaine"] == "Yes"
    source = net.links["Source"].str.split(":").str[0]
    target = net.links["Target"].str.split(":").str[0]
    assert (chained == ~(source.isin(["sex", "age"]) & target.isin(["sex", "age"]))).all()
    assert net.links["hidden"].tolist() == (~chained).tolist()


def test_exogenous_variables_work_with_a_filter_of_your_own(survey):
    rule = lambda links: links["Haberman"] > 0
    net = sur_coin(survey, ["sex", "age", "works"], exogenous="sex", link_filter=rule, max_link=1)
    assert (net.links["hidden"] == ~((net.links["Haberman"] > 0) & (net.links["chaine"] == "Yes"))).all()


def test_metric_variables_are_nodes_and_their_correlations_links(survey):
    net = sur_coin(survey, ["sex", "works"], metric="income", max_link=1)
    row = net.nodes.set_index("name").loc["income"]
    assert row["mean"] == pytest.approx(survey["income"].mean()) and row["min"] == survey["income"].min()
    assert row["%"] == pytest.approx((survey["income"].mean() - survey["income"].min())
                                     / (survey["income"].max() - survey["income"].min()) * 100)
    link = net.links[(net.links["Source"] == "works:yes") & (net.links["Target"] == "income")].iloc[0]
    r = np.corrcoef(survey["income"], survey["works"] == "yes")[0, 1]
    assert link["Haberman"] == pytest.approx(r * np.sqrt(300))
    assert link["p(Z)"] < 1e-6


def test_metric_procedures_are_limited_to_pearson_haberman_and_z(survey):
    net = sur_coin(survey, ["sex", "works"], metric="income", procedures=["Jaccard", "Pearson"],
                   criteria="Pearson", min_link=0, max_link=np.inf)
    assert list(net.links.columns) == ["Source", "Target", "Pearson", "p(Z)"]


def test_metric_correlations_with_weights_and_pairwise(survey):
    survey.loc[:4, "income"] = np.nan
    weighted = sur_coin(survey, ["sex", "works", "w"], metric="income", weight="w", max_link=1)
    pairwise = sur_coin(survey, ["sex", "works"], metric="income", pairwise=True, max_link=1)
    plain = sur_coin(survey, ["sex", "works"], metric="income", max_link=1)
    assert weighted.lcount() > 0
    # with no other missing value, the pairwise and the listwise correlations are the same
    pd.testing.assert_frame_equal(pairwise.links, plain.links)


def test_mc_and_pc_have_no_coordinates_for_the_metric_nodes(survey):
    """As in R: the analysis covers the categories only, so the layout is refused."""
    with pytest.warns(UserWarning, match="pair of coordinates per node"):
        net = sur_coin(survey, ["sex", "works"], metric="income", layout="mc", max_link=1)
    assert "fx" not in net.nodes.columns


def test_a_matrix_layout_with_metric_nodes_follows_the_nodes(survey):
    nodes = pd.DataFrame({"name": ["income", "works:yes", "sex:Man", "sex:Woman"]})
    layout = np.array([[9.0, 9.0], [1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    net = sur_coin(survey, ["sex", "works"], metric="income", nodes=nodes, layout=layout, max_link=1)
    got = net.nodes.set_index("name")[["fx", "fy"]].to_dict("index")
    assert got["income"] == {"fx": 9.0, "fy": 9.0} and got["sex:Man"] == {"fx": 3.0, "fy": 4.0}


def test_igraph_is_returned_when_asked(survey):
    igraph = pytest.importorskip("igraph")
    graph = sur_coin(survey, ["sex", "age"], igraph=True)
    assert isinstance(graph, igraph.Graph) and graph.vcount() == 5


def test_the_page_is_written_to_the_directory(survey, tmp_path):
    folder = tmp_path / "survey"
    sur_coin(survey, ["sex", "age"], directory=str(folder))
    assert (folder / "index.html").is_file()


def test_old_camel_case_names_work_with_a_warning(survey):
    with pytest.warns(DeprecationWarning, match="showArrows"):
        net = sur_coin(survey, ["sex", "age"], showArrows=True)
    assert net.options["showArrows"] is True


def test_the_label_keeps_the_variable_only_for_the_common_ones(survey):
    net = sur_coin(survey, ["sex", "age"], commonlabel=["sex"])
    labels = net.nodes.set_index("name")["label"].to_dict()
    assert labels["sex:Man"] == "sex:Man" and labels["age:old"] == "old"
    assert net.options["nodeLabel"] == "label"


def test_spanish_names_the_columns(survey):
    net = sur_coin(survey, ["sex", "age"], language="es", frequency=True)
    assert {"nombre", "frecuencia"} <= set(net.nodes.columns)


# --------------------------------------------------------------- cob_coin ----

def test_cob_coin_keeps_the_links_between_variables_that_are_significant(survey):
    net = cob_coin(survey, ["sex", "age", "works"], significance=0.2)
    source = net.links["Source"].str.split(":").str[0]
    target = net.links["Target"].str.split(":").str[0]
    assert (source != target).all()
    assert ((net.links["p(Z)"] < 0.1) | (net.links["p(Z)"] > 0.9)).all()
    assert net.options["nodeColor"] == "variable" and net.options["nodeGroup"] == "variable"
    assert net.options["linkBipolar"] is True and net.options["groupText"] is True


def test_cob_coin_puts_the_nodes_on_a_circle(survey):
    net = cob_coin(survey, ["sex", "age", "works"], degree=90)
    x, y = net.nodes["fx"].to_numpy(), net.nodes["fy"].to_numpy()
    assert np.allclose(np.hypot(x, y), 1, atol=1e-3)
    assert (x[0], y[0]) == pytest.approx((0, 1), abs=1e-3)


def test_cob_coin_writes_the_page_with_the_filtered_links(survey, tmp_path):
    folder = tmp_path / "coweb"
    net = cob_coin(survey, ["sex", "age", "works"], directory=str(folder))
    html = (folder / "index.html").read_text(encoding="utf-8")
    assert '"linknames"' in html and net.lcount() > 0


def test_cob_coin_igraph(survey):
    igraph = pytest.importorskip("igraph")
    assert isinstance(cob_coin(survey, ["sex", "age", "works"], igraph=True), igraph.Graph)


def test_cob_coin_needs_the_p_value_of_the_links(survey):
    with pytest.raises(ValueError, match="p-value"):
        cob_coin(survey, ["sex", "age"], criteria="Jaccard", procedures="Jaccard", min_link=0)
