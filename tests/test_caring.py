"""``caring``: reading files, building the equations and the page of networks."""
import json
import os
import warnings

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("statsmodels")
pytest.importorskip("sklearn")
from netcoin import caring, caring_create_graphs, caring_read_file
from netcoin.caring import glm_equations


@pytest.fixture
def survey():
    rng = np.random.default_rng(3)
    n = 300
    return pd.DataFrame({"Sex": rng.choice(["Man", "Woman"], n), "Age group": rng.choice(["young", "mid", "old"], n),
                         "Votes": rng.integers(0, 2, n), "Works": rng.integers(0, 2, n),
                         "Income": rng.normal(50, 10, n).round(1), "W": rng.uniform(.5, 2, n).round(2)})


def _page(folder):
    """The data of every graph of a written page, in order."""
    out = []
    index = 0
    while os.path.isdir(os.path.join(folder, "data", str(index))):
        html = open(os.path.join(folder, "data", str(index), "index.html"), encoding="utf-8").read()
        out.append(json.loads(html.split('id="data">')[1].split("</script>")[0]))
        index += 1
    return out


# ------------------------------------------------------------------ reading ----

def test_a_csv_is_read_with_categories_and_numbers(survey, tmp_path):
    survey.to_csv(tmp_path / "s.csv", index=False)
    data = caring_read_file(str(tmp_path / "s.csv"))
    assert isinstance(data["Sex"].dtype, pd.CategoricalDtype) and pd.api.types.is_numeric_dtype(data["Income"])


def test_text_made_of_numbers_becomes_numeric_and_quotes_are_changed(tmp_path):
    pd.DataFrame({'say "x"': ["1", "2", "3"], "name": ["a", "b", "c"]}).to_csv(tmp_path / "s.csv", index=False)
    data = caring_read_file(str(tmp_path / "s.csv"))
    assert list(data.columns) == ["say 'x'", "name"] and pd.api.types.is_numeric_dtype(data["say 'x'"])


def test_a_latin1_text_file_is_read(tmp_path):
    (tmp_path / "s.csv").write_bytes("edad;sexo\n20;hombre\n30;mujer\n40;mujer\n".encode("latin-1"))
    data = caring_read_file(str(tmp_path / "s.csv"))
    assert data["edad"].tolist() == [20, 30, 40]


def test_excel_and_stata_files_are_read(survey, tmp_path):
    pytest.importorskip("openpyxl")
    survey.to_excel(tmp_path / "s.xlsx", index=False)
    assert list(caring_read_file(str(tmp_path / "s.xlsx")).columns) == list(survey.columns)
    survey.to_stata(tmp_path / "s.dta", write_index=False, variable_labels={"Sex": "Sex of the person"})
    data = caring_read_file(str(tmp_path / "s.dta"))
    assert "Sex of the person" in data.columns and "Income" in data.columns


# -------------------------------------------------------------- the equations ----

def test_the_equations_follow_rs_ordering_of_variables_and_dummies(survey):
    frame, formulas = glm_equations(survey, ["Sex", "Age group", "Votes", "Income"], ["Sex"], [0, 1, 2, 3],
                                    ["gaussian", "binomial", "binomial", "gaussian"], ["Income"], ["Votes"], [1])
    assert formulas == ("Income ~ Votes:1+Agegroup:young+Agegroup:old+Agegroup:mid+Sex , gaussian\n"
                        "Votes:1 ~ Agegroup:young+Agegroup:old+Agegroup:mid+Sex , binomial\n"
                        "Agegroup:young ~ Sex , binomial\n"
                        "Agegroup:old ~ Sex , binomial\n"
                        "Agegroup:mid ~ Sex , binomial")
    assert {"Votes:1", "Agegroup:young", "Agegroup:old", "Agegroup:mid", "Sex", "Income"} == set(frame.columns)
    assert set(frame["Votes:1"].unique()) <= {0.0, 1.0}


def test_a_dichotomy_value_is_recycled_for_every_dichotomy(survey):
    frame, formulas = glm_equations(survey, ["Sex", "Votes", "Works", "Income"], ["Sex"], [0, 1, 1, 2],
                                    [None, "binomial", "binomial", "gaussian"], ["Income"], ["Votes", "Works"], [1])
    assert "Works:1" in formulas and "Works:NA" not in formulas      # R takes valueDicho[2], which is NA


def test_exogenous_variables_explain_and_are_never_explained(survey):
    _, formulas = glm_equations(survey, ["Sex", "Votes", "Income"], ["Sex", "Votes"], [0, 0, 1],
                                [None, None, "gaussian"], ["Income"], None, None)
    assert formulas == "Income ~ Sex+Votes , gaussian" or formulas == "Income ~ Votes+Sex , gaussian"
    assert formulas.count("~") == 1


# ----------------------------------------------------------------- the page ----

def test_each_plot_is_a_graph_of_the_page(survey, tmp_path):
    arguments = {"variables": ["Sex", "Age group", "Votes", "Income"], "weight": "W", "metric": ["Income"],
                 "exogenous": ["Sex"], "chaine": [0, 1, 2, 3], "valueDicho": 1, "dichotomies": ["Votes"],
                 "family": ["gaussian", "binomial", "binomial", "gaussian"], "nclustersMin": 2, "nclustersMax": 3,
                 "plot": ["surCoin", "Coincidences", "surScat", "Scattergram", "logCoin2", "Log-linear",
                          "glmCoin", "Regression"],
                 "dir": str(tmp_path / "page")}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        folder = caring_create_graphs(survey, arguments)
    assert folder == str(tmp_path / "page") and os.path.isfile(os.path.join(folder, "index.html"))
    pages = _page(folder)
    assert len(pages) == 4
    names = json.loads(json.dumps(pages[0]["options"]["multigraph"]))["names"]
    assert names == ["Coincidences", "Scattergram", "Log-linear", "Regression"]
    assert "links" not in pages[1] and "fx" in pages[1]["nodenames"]                 # the scattergram has no links
    assert "</br>" in pages[3]["options"]["note"] and "Income ~" in pages[3]["options"]["note"]


def test_the_scattergram_follows_the_numbers_of_groups_of_the_dialog(survey, tmp_path):
    arguments = {"variables": ["Sex", "Age group", "Votes"], "plot": ["surScat", "S"], "nclustersMin": 2,
                 "nclustersMax": 4, "nclustersDef": 3, "dir": str(tmp_path / "p")}
    caring_create_graphs(survey, arguments)
    nodes = _page(str(tmp_path / "p"))[0]
    columns = [n for n in nodes["nodenames"] if n.startswith("Groups(")]
    assert columns == ["Groups(2)", "Groups(4)", "Groups(3)"]


def test_only_metric_variables_make_a_network_of_correlations(survey, tmp_path):
    survey["Age"] = np.random.default_rng(0).normal(40, 10, 300)
    arguments = {"variables": ["Income", "Age"], "metric": ["Income", "Age"], "plot": ["surCoin", "Correlations"],
                 "weight": "W", "dir": str(tmp_path / "p")}
    caring_create_graphs(survey, arguments)
    page = _page(str(tmp_path / "p"))[0]
    assert page["nodenames"][:3] == ["name", "mean", "std"] and page["options"]["linkBipolar"] is True


def test_a_weight_that_is_not_numeric_is_ignored_and_dichotomies_are_found(survey, tmp_path):
    arguments = {"variables": ["Sex", "Votes", "Works"], "weight": "Sex", "plot": ["surCoin", "Net"],
                 "dir": str(tmp_path / "p")}
    caring_create_graphs(survey, arguments)
    nodes = _page(str(tmp_path / "p"))[0]
    names = nodes["nodes"][nodes["nodenames"].index("name")]
    assert "Votes" in names and "Works" in names          # the 0/1 columns are dichotomies, nodes of their own


def test_the_folder_is_temporary_when_none_is_given(survey):
    arguments = {"variables": ["Sex", "Age group"], "plot": ["surCoin", "Net"]}
    folder = caring_create_graphs(survey, arguments)
    assert os.path.isfile(os.path.join(folder, "index.html"))


def test_allnet_needs_incidences(survey, tmp_path):
    ok = {"variables": ["Votes", "Works"], "plot": ["allNet", "Inc"], "dir": str(tmp_path / "a")}
    caring_create_graphs(survey, ok)
    assert len(_page(str(tmp_path / "a"))) == 1
    with pytest.raises(ValueError, match="allNet needs incidences"):
        caring_create_graphs(survey, {"variables": ["Sex", "Votes"], "plot": ["allNet", "Inc"], "dir": str(tmp_path / "b")})


def test_the_regression_needs_independent_variables_and_a_plot_is_required(survey, tmp_path):
    with pytest.raises(ValueError, match="missing independent variables"):
        caring_create_graphs(survey, {"variables": ["Sex", "Income"], "plot": ["glmCoin", "R"], "chaine": [0, 1],
                                      "family": ["gaussian", "gaussian"]})
    with pytest.raises(ValueError, match="plot must name"):
        caring_create_graphs(survey, {"variables": ["Sex", "Income"]})


def test_caring_reads_the_file_and_writes_the_page(survey, tmp_path):
    survey.to_csv(tmp_path / "s.csv", index=False)
    folder = caring(str(tmp_path / "s.csv"), {"variables": ["Sex", "Age group"], "plot": ["surCoin", "Net"],
                                             "dir": str(tmp_path / "out")})
    assert os.path.isfile(os.path.join(folder, "index.html"))
