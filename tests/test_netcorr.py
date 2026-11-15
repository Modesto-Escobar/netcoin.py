"""``net_corr``, ``d_net_corr`` and ``layout_fact``: inputs, errors and independent checks."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from netcoin import d_net_corr, layout_fact, net_corr
from netcoin.netcorr import _component_labels, _corr_criteria, _corr_methods, corr_stat


@pytest.fixture
def data():
    rng = np.random.default_rng(11)
    n = 80
    base = rng.normal(size=(n, 2))
    return pd.DataFrame({
        "a": base[:, 0] + rng.normal(0, .5, n),
        "b": base[:, 0] + rng.normal(0, 1, n),
        "c": base[:, 1] + rng.normal(0, .5, n),
        "d": base[:, 1] - base[:, 0] + rng.normal(0, 1, n),
    })


def _link(net, a, b):
    links = net.links
    row = links[((links["Source"] == a) & (links["Target"] == b)) | ((links["Source"] == b) & (links["Target"] == a))]
    return row.iloc[0]


# -------------------------------------------------------------- arguments ----

def test_methods_can_be_abbreviated():
    assert _corr_methods("p") == ["pearson"] and _corr_methods(["s", "ken", "PA", "d"]) == \
        ["spearman", "kendall", "partial", "dcor"]
    assert _corr_methods(["pearson", "p"]) == ["pearson"]
    with pytest.raises(ValueError, match="not available"):
        _corr_methods("pe_")                      # no method starts like that
    with pytest.raises(ValueError, match="not available"):
        _corr_methods("x")


def test_criteria_name_a_column_of_the_links():
    assert _corr_criteria("P", ["pearson"]) == "p"
    assert _corr_criteria("pearson", ["pearson", "kendall"]) == "value"
    assert _corr_criteria("z(pearson)", ["pearson"]) == "z"
    assert _corr_criteria("p(kendall)", ["pearson", "kendall"]) == "p(kendall)"
    with pytest.raises(ValueError, match="criteria"):
        _corr_criteria("z(kendall)", ["pearson"])


def test_the_variables_must_be_numeric_and_enough(data):
    with pytest.raises(ValueError, match="quantitative: t"):
        net_corr(data.assign(t="x"))
    with pytest.raises(ValueError, match="at least two"):
        net_corr(data[["a"]])
    with pytest.raises(ValueError, match="three complete"):
        net_corr(data.head(2))


def test_weights_are_checked(data):
    with pytest.raises(ValueError, match="weight must be"):
        net_corr(data, weight="nope")
    with pytest.raises(ValueError, match="weight must be"):
        net_corr(data, weight=[1, 2, 3])
    with pytest.raises(ValueError, match="negative"):
        net_corr(data, weight=-np.ones(len(data)))
    with pytest.raises(ValueError, match="not available for the methods"):
        net_corr(data, weight=np.ones(len(data)), method="kendall")


def test_a_weight_column_is_not_a_variable(data):
    net = net_corr(data.assign(w=1.0), weight="w", max_link=1)
    assert net.nodes["name"].tolist() == ["a", "b", "c", "d"]


def test_missing_values_leave_out_cases_or_pairs(data):
    data.loc[0, "a"] = np.nan
    data.loc[1, "b"] = np.nan
    listwise = net_corr(data, frequency=True, max_link=1)
    pairwise = net_corr(data, pairwise=True, max_link=1)
    assert listwise.nodes["frequency"].tolist() == [78] * 4 and "n" not in listwise.links
    assert pairwise.nodes["frequency"].tolist() == [79, 79, 80, 80]
    assert _link(pairwise, "a", "b")["n"] == 78 and _link(pairwise, "a", "c")["n"] == 79


def test_cases_with_a_missing_weight_are_left_out(data):
    w = np.ones(len(data))
    w[:5] = np.nan
    net = net_corr(data, weight=w, frequency=True, max_link=1)
    assert net.nodes["frequency"].tolist() == [75] * 4


def test_links_are_selected_by_criteria_and_range(data):
    everything = net_corr(data, max_link=1)
    assert everything.lcount() == 6
    strong = net_corr(data, criteria="value", min_link=0.5)
    assert (strong.links["value"].abs() >= 0.5).all() and strong.lcount() < 6
    one_sided = net_corr(data, link_bipolar=False)
    assert (one_sided.links["p"] <= 0.5).all()
    bonferroni = net_corr(data, bonferroni=True, max_link=0.05)
    assert (bonferroni.links["p"] <= 0.05 / 6).all()


def test_links_can_be_sorted(data):
    net = net_corr(data, max_link=1, sort_links="value", decreasing_links=False)
    assert net.links["value"].is_monotonic_increasing


def test_nodes_can_be_sorted_by_cases_then_mean(data):
    data.loc[:9, "c"] = np.nan
    net = net_corr(data, pairwise=True, sort=True, frequency=True, max_link=1)
    assert net.nodes["frequency"].tolist() == sorted(net.nodes["frequency"], reverse=True)
    assert net.nodes["name"].iloc[-1] == "c"


def test_the_defaults_make_bipolar_links_coloured_by_value(data):
    net = net_corr(data)
    assert net.options["linkBipolar"] is True
    assert (net.options["linkWidth"], net.options["linkColor"], net.options["linkWeight"]) == ("value",) * 3


def test_nothing_selected_gives_a_network_without_links(data):
    net = net_corr(data, criteria="value", min_link=2)
    assert net.links is None and net.ncount() == 4


def test_pc_layout_follows_the_names_of_the_nodes(data):
    nodes = pd.DataFrame({"name": ["d", "a"], "g": ["x", "y"]})
    plain = net_corr(data, layout="pc").nodes.set_index("name")
    chosen = net_corr(data, nodes=nodes, layout="pc").nodes.set_index("name")
    assert chosen.index.tolist() == ["a", "d", "b", "c"]          # those in nodes go first
    for name in chosen.index:
        assert chosen.loc[name, ["fx", "fy"]].tolist() == pytest.approx(plain.loc[name, ["fx", "fy"]].tolist())


def test_igraph_is_returned_when_asked(data):
    igraph = pytest.importorskip("igraph")
    assert isinstance(net_corr(data, igraph=True), igraph.Graph)


def test_the_page_is_written_to_the_directory(data, tmp_path):
    net_corr(data, directory=str(tmp_path / "corr"))
    assert (tmp_path / "corr" / "index.html").is_file()


# --------------------------------------------- independent checks of the numbers ----

def test_pearson_matches_scipy(data):
    net = net_corr(data, max_link=1)
    row = _link(net, "a", "b")
    r, p_two = stats.pearsonr(data["a"], data["b"])
    assert row["value"] == pytest.approx(r) and row["p"] == pytest.approx(p_two)
    assert _link(net_corr(data, link_bipolar=False, max_link=1), "a", "b")["p"] == pytest.approx(p_two / 2)


def test_spearman_and_kendall_match_scipy_with_many_cases(data):
    net = net_corr(data, method=["spearman", "kendall"], max_link=1)
    row = _link(net, "a", "d")
    assert row["value"] == pytest.approx(stats.spearmanr(data["a"], data["d"]).statistic)
    assert row["p"] == pytest.approx(stats.spearmanr(data["a"], data["d"]).pvalue, rel=1e-6)
    assert row["kendall"] == pytest.approx(stats.kendalltau(data["a"], data["d"]).statistic)
    assert row["p(kendall)"] == pytest.approx(stats.kendalltau(data["a"], data["d"]).pvalue, rel=0.02)


def test_kendall_exact_test_for_small_samples():
    rng = np.random.default_rng(1)
    small = pd.DataFrame(rng.normal(size=(8, 3)), columns=list("xyz"))
    net = net_corr(small, method="kendall", max_link=1)
    for a, b in (("x", "y"), ("x", "z"), ("y", "z")):
        want = stats.kendalltau(small[a], small[b], method="exact")
        assert _link(net, a, b)["p"] == pytest.approx(want.pvalue, rel=1e-9)


def test_spearman_exact_test_for_small_samples():
    rng = np.random.default_rng(2)
    small = pd.DataFrame(rng.normal(size=(7, 3)), columns=list("xyz"))
    net = net_corr(small, method="spearman", max_link=1)
    for a, b in (("x", "y"), ("x", "z"), ("y", "z")):
        y = small[b].to_numpy()
        want = stats.permutation_test((small[a].to_numpy(),), lambda v: stats.spearmanr(v, y).statistic,
                                      permutation_type="pairings", n_resamples=np.inf)
        assert _link(net, a, b)["p"] == pytest.approx(want.pvalue, rel=1e-9)


def test_partial_correlation_is_the_inverse_of_the_correlation_matrix(data):
    net = net_corr(data, method="partial", max_link=1)
    precision = np.linalg.inv(np.corrcoef(data.to_numpy(), rowvar=False))
    expected = -precision[0, 1] / np.sqrt(precision[0, 0] * precision[1, 1])
    row = _link(net, "a", "b")
    assert row["value"] == pytest.approx(expected)
    t = expected * np.sqrt((len(data) - 4) / (1 - expected ** 2))
    assert row["p"] == pytest.approx(2 * stats.t.sf(abs(t), len(data) - 4))


def test_weights_of_one_change_nothing_and_integer_weights_repeat_cases(data):
    plain = net_corr(data, max_link=1)
    ones = net_corr(data, weight=np.ones(len(data)), max_link=1)
    assert ones.links["value"].tolist() == pytest.approx(plain.links["value"].tolist())
    counts = np.arange(1, len(data) + 1) % 3 + 1
    repeated = pd.DataFrame(np.repeat(data.to_numpy(), counts, axis=0), columns=data.columns)
    weighted = net_corr(data, weight=counts.astype(float), max_link=1)
    again = net_corr(repeated, max_link=1)
    assert weighted.links["value"].tolist() == pytest.approx(again.links["value"].tolist())
    assert weighted.nodes["mean"].tolist() == pytest.approx(again.nodes["mean"].tolist(), abs=0.006)


def test_distance_correlation_of_a_dependent_pair_is_high():
    x = np.linspace(-3, 3, 60)
    frame = pd.DataFrame({"x": x, "y": x ** 2, "z": np.sin(np.arange(60) * 7.0)})
    net = net_corr(frame, method="dcor", max_link=1)
    row = _link(net, "x", "y")
    assert row["value"] > 0.15 and row["p"] < 0.05       # Pearson sees nothing in it:
    assert abs(stats.pearsonr(frame["x"], frame["y"]).statistic) < 1e-9


def test_corr_stat_gives_one_value_per_pair_column_by_column():
    X = np.random.default_rng(0).normal(size=(20, 4))
    result = corr_stat(X, None, "pearson", False, True)
    assert len(result["p"]) == 6 and result["R"].shape == (4, 4)
    assert (result["n"] == 20).all()


# ------------------------------------------------------------------ d_net_corr ----

@pytest.fixture
def blocks():
    rng = np.random.default_rng(5)
    lat = rng.normal(size=(60, 3))
    cols = {}
    for k, sds in enumerate(([.5, .8, 1.2], [.5, .9, 1.0], [.6, 1.1])):
        for j, sd in enumerate(sds):
            cols[f"{'abc'[k]}{j + 1}"] = lat[:, k] + rng.normal(0, sd, 60)
    return pd.DataFrame(cols)


def test_component_labels_are_numbered_by_size():
    assert _component_labels([1, 1, 2, 3, 3, 3]) == ["C-2", "C-2", "C-3", "C-1", "C-1", "C-1"]
    many = list(range(1, 11)) + [1]
    labels = _component_labels(many)
    assert labels[0] == "C- 1" and labels[1] == "C- 2" and labels[-1] == "C- 1"


def test_d_net_corr_makes_frames_that_follow_the_components(blocks):
    net = d_net_corr(blocks, sequence=np.round(np.arange(20, 101, 5) / 100, 2))
    frames = net.options["frames"]
    assert len(frames) >= 2 and frames[0] == "C0.2"
    assert net.options["speed"] == 50 and net.options["frame"] == 0
    assert {"component", "Degree", "ratio", "closeness", "betweenness", "eigen"} <= set(net.nodes.columns) | set(net.links.columns)
    assert set(net.links["_frame_"]) == set(range(len(frames)))
    first = net.links[net.links["_frame_"] == 0]
    assert (first["ratio"] >= 1).all()


def test_d_net_corr_keeps_every_frame_without_components(blocks):
    seq = [0.3, 0.5, 0.7]
    every = d_net_corr(blocks, components=False, sequence=seq)
    assert every.options["frames"] == ["C0.3", "C0.5", "C0.7"]


def test_d_net_corr_needs_at_least_two_frames(blocks):
    with pytest.raises(ValueError, match="only one graph"):
        d_net_corr(blocks, sequence=[0.3])


def test_d_net_corr_text_shows_the_strong_values_only(blocks):
    seq = [0.3, 0.6]
    net = d_net_corr(blocks, components=False, sequence=seq, text_filter=(0.5, 0.1), means=False)
    links = net.links
    for frame, threshold in enumerate(seq):
        part = links[links["_frame_"] == frame]
        limit = min(0.5, threshold + 0.1)
        assert (part["text"] != "").tolist() == (part["value"].abs() >= limit).tolist()
        assert part.loc[part["text"] != "", "text"].tolist() == [f"{v:.2f}" for v in part.loc[part["text"] != "", "value"]]


def test_d_net_corr_warns_if_layout_and_factorial_are_given(blocks):
    layout = np.random.default_rng(0).normal(size=(blocks.shape[1], 2))
    with pytest.warns(UserWarning, match="incompatible"):
        d_net_corr(blocks, layout=layout, factorial="pc", components=False, sequence=[0.3, 0.5])


def test_d_net_corr_writes_the_page(blocks, tmp_path):
    d_net_corr(blocks, components=False, sequence=[0.3, 0.5], directory=str(tmp_path / "evolving"))
    assert (tmp_path / "evolving" / "index.html").is_file()


# --------------------------------------------------------------- layout_fact ----

def test_principal_components_are_unit_eigenvectors(data):
    plane = layout_fact(data, "pc")
    assert list(plane.columns) == ["Comp.1", "Comp.2"]
    assert (plane ** 2).sum().tolist() == pytest.approx([1, 1])


def test_layout_fact_leaves_out_incomplete_cases(data):
    data.loc[0, "a"] = np.nan
    assert layout_fact(data, "pc").equals(layout_fact(data.dropna(), "pc"))


def test_rotated_factors_keep_the_communalities(data):
    plain, varimax = layout_fact(data, "nf"), layout_fact(data, "vf")
    assert (plain ** 2).sum(axis=1).tolist() == pytest.approx((varimax ** 2).sum(axis=1).tolist(), rel=1e-6)
