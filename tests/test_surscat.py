"""``sur_scat``: validations, patterns, samples and what the object keeps."""
import numpy as np
import pandas as pd
import pytest

from netcoin import sur_scat
from netcoin.surscat import cluster_shares, hartigan_refine, mca_rows, pattern_groups, pca_scores


@pytest.fixture
def survey():
    rng = np.random.default_rng(7)
    n = 240
    group = rng.choice(["a", "b", "c"], n)
    frame = pd.DataFrame({
        "q1": np.where(group == "a", rng.choice(["yes", "no"], n, p=[.9, .1]), rng.choice(["yes", "no"], n)),
        "q2": np.where(group == "b", rng.choice(["hi", "lo"], n, p=[.9, .1]), rng.choice(["hi", "lo"], n)),
        "q3": rng.choice(["x", "y", "z"], n),
        "age": rng.integers(18, 80, n).astype(float),
        "inc": rng.normal(50, 10, n),
        "w": rng.uniform(.5, 2, n),
    })
    return frame


QUALI = ["q1", "q2", "q3"]


# ---------------------------------------------------------------- analyses ----

def test_mca_gives_the_non_trivial_factors_and_their_eigenvalues():
    m = np.array([[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 0, 1], [0, 1, 1, 0], [1, 0, 1, 0]], float)
    coords, eig = mca_rows(m, 2)
    assert coords.shape == (5, 2) and len(eig) >= 2 and (eig >= -1e-12).all()


def test_pca_scores_are_uncorrelated_and_have_the_eigenvalues_as_variance():
    x = np.random.default_rng(0).normal(size=(60, 3)) @ np.array([[1, .5, 0], [0, 1, .2], [0, 0, 1]])
    scores, eig, std = pca_scores(x, None, 3)
    np.testing.assert_allclose(np.var(scores, axis=0, ddof=1), eig, rtol=1e-9)
    assert abs(np.corrcoef(scores.T)[0, 1]) < 1e-9


def test_a_weighted_pca_with_equal_weights_is_the_plain_one():
    x = np.random.default_rng(1).normal(size=(50, 3))
    plain, eig, _ = pca_scores(x, None, 2)
    weighted, eig_w, _ = pca_scores(x, np.full(50, 2.0), 2)
    np.testing.assert_allclose(eig_w, eig, rtol=1e-9)
    # the weighted standardization divides by the population deviation, the plain one by n - 1
    np.testing.assert_allclose(np.abs(weighted), np.abs(plain) * np.sqrt(50 / 49), rtol=1e-6)


def test_the_refinement_never_raises_the_sum_of_squares():
    rng = np.random.default_rng(2)
    x = rng.normal(size=(120, 2))
    start = rng.integers(0, 3, 120)

    def sse(l):
        return sum(((x[l == g] - x[l == g].mean(axis=0)) ** 2).sum() for g in range(3))
    assert sse(hartigan_refine(x, start, 3)) <= sse(start)


def test_pattern_groups_lead_with_the_groups_and_follow_with_combinations():
    out = pattern_groups(["G: 2|G: 1", "G: 1", "G: 3"], ["G: 1", "G: 2", "G: 3"])
    assert list(out.categories) == ["G: 1", "G: 2", "G: 3", "G: 2|G: 1"] and out.ordered


def test_cluster_shares_add_the_cases_of_each_group_in_each_node():
    group = pd.Categorical(["u", "v", "u", "u"], categories=["u", "v"])
    shares = cluster_shares(np.array([1, 1, 2, 2]), 2, group, np.array([1.0, 2.0, 3.0, 4.0]))
    assert shares.to_numpy().tolist() == [[1.0, 2.0], [7.0, 0.0]]


# -------------------------------------------------------------- the scattergram ----

def test_cases_are_nodes_with_coordinates_groups_and_no_links(survey):
    net = sur_scat(survey, QUALI, nclusters=3)
    assert net.ncount() == 240 and net.links is None
    assert {"q1", "q2", "q3", "Groups(3)", "name", "fx", "fy"} <= set(net.nodes.columns)
    assert net.options["nodeColor"] == "Groups(3)" and net.options["showAxes"] is True
    assert net.options["axesLabels"][0].startswith("F1 (") and net.scat["plane_name"] == "MCA"
    assert list(net.nodes["Groups(3)"].cat.categories) == ["Group: 1", "Group: 2", "Group: 3"]


def test_groups_are_numbered_from_left_to_right(survey):
    net = sur_scat(survey, QUALI, nclusters=3)
    means = net.nodes.groupby("Groups(3)", observed=True)["fx"].mean()
    assert means.is_monotonic_increasing
    unsorted = sur_scat(survey, QUALI, nclusters=3, sort_clusters=False)
    assert unsorted.nodes["Groups(3)"].nunique() == 3


def test_the_type_is_pca_when_the_active_variables_are_numbers(survey):
    net = sur_scat(survey, ["age", "inc", "q1"], active=["age", "inc"], nclusters=2)
    assert net.scat["plane_name"] == "PCA" and net.options["axesLabels"][0].startswith("PC1 (")
    with pytest.raises(ValueError, match="type=\"pca\" requires at least two"):
        sur_scat(survey, ["age", "q1"], active=["age"], type="pca")


def test_observed_variables_can_be_the_axes(survey):
    net = sur_scat(survey, ["age", "inc", "q1"], active=["age", "inc"], xaxis="age", yaxis="inc")
    assert net.options["axesLabels"] == ["age", "inc"]
    assert not any(c.startswith("Groups(") for c in net.nodes.columns)            # no groups by default
    standard = sur_scat(survey, ["age", "inc", "q1"], active=["age", "inc"], xaxis=("age", True), yaxis="inc", nclusters=2)
    assert abs(standard.nodes["fx"].mean()) < 1e-9 and standard.nodes["fx"].std() == pytest.approx(1)
    assert standard.nodes["fy"].to_numpy() == pytest.approx(survey["inc"].to_numpy())


def test_a_nominal_variable_cannot_be_an_axis(survey):
    with pytest.raises(ValueError, match="nominal variable"):
        sur_scat(survey, QUALI, xaxis="q1")


def test_factor_axes_are_read_in_several_ways(survey):
    net = sur_scat(survey, QUALI, nfactors=3, xaxis="F1", yaxis="dim3", nclusters=2)
    assert net.options["axesLabels"][1].startswith("F3 (") and "plane_name" not in net.scat
    with pytest.raises(ValueError, match="factorial axis"):
        sur_scat(survey, QUALI, xaxis="zzz")


def test_nclusters_is_checked(survey):
    with pytest.raises(ValueError, match="stand alone"):
        sur_scat(survey, QUALI, nclusters=[0, 2])
    with pytest.raises(ValueError, match="whole numbers"):
        sur_scat(survey, QUALI, nclusters=2.5)
    with pytest.warns(UserWarning, match="nclusters=1 makes no groups"):
        net = sur_scat(survey, QUALI, nclusters=1)
    assert not any(c.startswith("Groups(") for c in net.nodes.columns)
    with pytest.warns(UserWarning, match="1 group left out"):
        sur_scat(survey, QUALI, nclusters=[1, 2])


def test_several_numbers_make_several_columns_and_the_last_colours(survey):
    net = sur_scat(survey, QUALI, nclusters=[2, 4, 3])
    assert [c for c in net.nodes.columns if c.startswith("Groups(")] == ["Groups(2)", "Groups(4)", "Groups(3)"]
    assert net.options["nodeColor"] == "Groups(3)"


def test_auto_chooses_by_the_index_and_keeps_the_table(survey):
    net = sur_scat(survey, QUALI, nclusters="auto")
    index = net.scat["cluster_index"]
    assert list(index["k"]) == list(range(2, 11)) and index["optimal"].sum() == 1
    assert [c for c in net.nodes.columns if c.startswith("Groups(")] == [f"Groups({index.loc[index['optimal'], 'k'].item():02d})"
                                                                          if False else f"Groups({int(index.loc[index['optimal'], 'k'].iloc[0])})"]
    both = sur_scat(survey, QUALI, nclusters=["auto", 2, 3], cluster_index="silhouette")
    assert [c for c in both.nodes.columns if c.startswith("Groups(")] == ["Groups(2)", "Groups(3)"]
    assert both.scat["cluster_index"].attrs["index"] == "silhouette"


def test_the_bic_falls_back_to_calinski_harabasz(survey):
    with pytest.warns(UserWarning, match="not available yet"):
        net = sur_scat(survey, QUALI, nclusters=["auto", 2, 3], cluster_index="bic")
    assert net.scat["cluster_index"].attrs["index"] == "Calinski-Harabasz"


def test_what_is_not_ported_yet_is_said(survey):
    with pytest.warns(UserWarning, match="layouts"):
        sur_scat(survey, QUALI, nclusters=2, layouts=["pca"])


def test_weights_and_missing_values(survey):
    survey.loc[:5, "q1"] = np.nan
    net = sur_scat(survey, QUALI, weight="w", nclusters=2)
    assert net.ncount() == 234 and net.nodes["name"].iloc[0] == "7"
    with pytest.raises(ValueError, match="weight column"):
        sur_scat(survey, QUALI, weight="nope")
    with pytest.raises(ValueError, match="Weights"):
        sur_scat(survey, QUALI, weight=[1.0, 2.0])


def test_too_few_complete_cases_are_reported(survey):
    survey.loc[:150, "q1"] = np.nan
    with pytest.warns(UserWarning, match="only 89 of 240 cases are complete"):
        sur_scat(survey, QUALI, nclusters=2)


def test_patterns_merge_the_cases_with_the_same_answers(survey):
    net = sur_scat(survey, QUALI, patterns=True, nclusters=2)
    assert net.ncount() == survey[QUALI].drop_duplicates().shape[0] and "n" in net.nodes.columns
    assert net.nodes["n"].sum() == 240 and net.options["nodeSize"] == "n" and net.options["nodeText"] == "ntext"
    assert net.scat["case_to_pattern"].max() == net.ncount()
    assert net.nodes["name"].iloc[0] == "001" or net.nodes["name"].iloc[0] == "01"
    assert "<b>q1</b>" in net.nodes["ntext"].iloc[0]


def test_patterns_on_some_variables_keep_the_other_columns_as_summaries(survey):
    net = sur_scat(survey, QUALI + ["age"], v_patterns=["q1", "q2"], nclusters=2)
    assert net.ncount() == survey[["q1", "q2"]].drop_duplicates().shape[0]
    assert (net.nodes["age"].between(18, 80)).all()
    assert net.nodes["q3"].str.contains(r"x|y|z").all()
    assert net.nodes.columns[:3].tolist() == ["q1", "q2", "n"]


def test_a_pattern_can_span_several_groups(survey):
    net = sur_scat(survey, QUALI, v_patterns=["q1"], nclusters=3)
    column = net.nodes["Groups(3)"].astype(str)
    assert column.str.contains(r"\|").any()
    assert set(net.scat["cluster_shares"]["Groups(3)"].columns) <= {"Group: 1", "Group: 2", "Group: 3"}


def test_a_sample_of_the_nodes_is_drawn_when_there_are_too_many(survey):
    net = sur_scat(survey, QUALI, nclusters=2, max_n=50)
    assert net.ncount() == 50 and net.scat["sampled_from"] == 240 and len(net.scat["sampled_nodes"]) == 50
    assert list(net.scat["sampled_nodes"]) == sorted(net.scat["sampled_nodes"])
    again = sur_scat(survey, QUALI, nclusters=2, max_n=50)
    assert list(again.scat["sampled_nodes"]) == list(net.scat["sampled_nodes"])


def test_jitter_moves_the_points_a_little(survey):
    plain = sur_scat(survey, QUALI, nclusters=2)
    moved = sur_scat(survey, QUALI, nclusters=2, jitter=True)
    shift = (moved.nodes["fx"] - plain.nodes["fx"]).abs()
    assert 0 < shift.max() <= 0.02 * (plain.nodes["fx"].max() - plain.nodes["fx"].min()) / 2 + 1e-12


def test_centroids_are_the_means_of_the_groups_on_the_plane(survey):
    net = sur_scat(survey, QUALI + ["age"], active=QUALI, nclusters=3)
    assert sur_scat(survey, QUALI, nclusters=3).cluster_means is None       # no numeric variable to average
    centroids = net.clusters["Groups(3).MCA"]
    by_hand = net.nodes.groupby("Groups(3)", observed=True)[["fx", "fy"]].mean()
    np.testing.assert_allclose(centroids.to_numpy(), by_hand.to_numpy())
    means = net.cluster_means["Groups(3)"]
    assert list(means.columns) == ["age"]
    np.testing.assert_allclose(means["age"], net.nodes.groupby("Groups(3)", observed=True)["age"].mean())


def test_a_seed_makes_the_result_repeatable(survey):
    a, b = sur_scat(survey, QUALI, nclusters=3, seed=5), sur_scat(survey, QUALI, nclusters=3, seed=5)
    pd.testing.assert_frame_equal(a.nodes, b.nodes)


def test_columns_decide_the_order_of_the_node_table(survey):
    net = sur_scat(survey, QUALI + ["age"], active=QUALI, columns=["_nclusters", "age"], nclusters=2)
    assert net.nodes.columns[:2].tolist() == ["Groups(2)", "age"]
    with pytest.warns(UserWarning, match="not a column"):
        sur_scat(survey, QUALI, columns=["nothing"], nclusters=2)


def test_the_language_names_groups_and_cases(survey):
    net = sur_scat(survey, QUALI, nclusters=2, language="es", patterns=True)
    assert "Grupos(2)" in net.nodes.columns and "nombre" in net.nodes.columns
    assert net.nodes["Grupos(2)"].iloc[0].startswith("Grupo: ") and "<b>Casos</b>" in net.nodes["ntext"].iloc[0]


def test_repeated_variables_are_kept_once(survey):
    with pytest.warns(UserWarning, match="kept once: q1"):
        net = sur_scat(survey, ["q1", "q2", "q1", "q3"], nclusters=2)
    assert net.nodes.columns.tolist().count("q1") == 1


def test_the_page_can_be_written(survey, tmp_path):
    sur_scat(survey, QUALI, nclusters=2, directory=str(tmp_path / "scat"))
    assert (tmp_path / "scat" / "index.html").is_file()
