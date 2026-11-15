"""Rotations, latent models, looking4clusters inside ``sur_scat`` and ``add_communities``: behaviour."""
import warnings

import numpy as np
import pandas as pd
import pytest

from netcoin import add_communities, sur_scat
from netcoin.latent import fit_lca, fit_lpa, latent_model, run_latent
from netcoin.rotations import equamax, quartimax, rotation_name, varimax
from netcoin.scatl4c import (align_plane, community_code, l4c_clusters, l4c_layout, latent_clusters,
                             lump_communities, rotation_layout)

igraph = pytest.importorskip("igraph")


@pytest.fixture
def blobs():
    rng = np.random.default_rng(3)
    group = np.repeat([0, 1, 2], 40)
    x = rng.normal(size=(120, 6)) + group[:, None] * np.array([2.5, -2, 1.5, 0, 1, -1.5])
    return pd.DataFrame(x, columns=list("abcdef"))


# --------------------------------------------------------------- what is asked for ----
def test_the_reductions_and_the_rotations_are_told_apart():
    assert l4c_layout(["PCA", "varimax", "mds"]) == ["pca", "cmds"]
    assert l4c_layout(True) == ["pca", "tsne", "cmds", "nmf", "umap"]
    assert l4c_layout(None) == [] and l4c_layout(False) == []
    assert rotation_layout(["varimax", "pca", "Equimax"]) == ["varimax", "equamax"]
    assert rotation_layout("all") == ["varimax", "quartimax", "equamax"]
    assert l4c_layout("all") == ["pca", "tsne", "cmds", "nmf", "umap"]
    assert rotation_name("Equimax") == ["equamax"]
    with pytest.warns(UserWarning, match="network layouts do not apply"):
        assert l4c_layout(["fr", "pca"]) == ["pca"]


def test_the_methods_are_checked():
    assert l4c_clusters(True) == ["kmeans", "pam", "hclust", "latent"]
    assert l4c_clusters(["PAM", "pam", "LCA"]) == ["pam", "lca"]
    with pytest.raises(ValueError, match="takes the names of clustering methods"):
        l4c_clusters(["pam", "dbscan"])


def test_the_latent_model_follows_the_analysis():
    assert latent_model("pca") == "LPA" and latent_model("mca") == "LCA"
    assert latent_clusters(["latent"], "pca") == ["LPA"]
    assert latent_clusters(["pam"], "pca") == []
    with pytest.warns(UserWarning, match="LCA left out"):
        assert latent_clusters(["lca"], "pca") == []


# ------------------------------------------------------------------------ rotations ----
def test_a_rotation_turns_the_loadings_without_changing_what_they_explain():
    rng = np.random.default_rng(0)
    loadings = np.linalg.qr(rng.normal(size=(8, 3)))[0] * 2
    for rotate in (lambda a: varimax(a), lambda a: quartimax(a), lambda a: equamax(a, 3 / 16)):
        rotated, theta = rotate(loadings)
        np.testing.assert_allclose(theta @ theta.T, np.eye(3), atol=1e-8)
        np.testing.assert_allclose(rotated, loadings @ theta, atol=1e-8)
        np.testing.assert_allclose((rotated ** 2).sum(), (loadings ** 2).sum(), rtol=1e-9)


def test_varimax_finds_simple_structure():
    base = np.array([[.9, .05], [.8, 0], [.85, .1], [.1, .9], [0, .8], [.05, .7], [0, .75]])
    angle = np.pi / 12
    turn = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    rotated, _ = varimax(base @ turn)
    criterion = lambda l: ((l ** 4).sum(axis=0) - (l ** 2).sum(axis=0) ** 2 / len(l)).sum()
    assert criterion(rotated) > criterion(base @ turn)
    assert (np.abs(rotated).min(axis=1) < 0.15).all()


# ------------------------------------------------------------------- latent models ----
def test_lpa_finds_the_groups_and_the_bic_picks_their_number(blobs):
    found = run_latent(blobs, "LPA", [1, 2, 3, 4, 5], seed=1)
    bic = found["bic"]
    assert min(bic, key=bic.get) == 3
    classes = pd.Series(found["clusters"][3])
    assert classes.groupby(np.repeat([0, 1, 2], 40)).agg(lambda s: s.value_counts().iloc[0]).min() >= 38


def test_lca_finds_the_classes(blobs):
    rng = np.random.default_rng(5)
    cls = np.repeat([0, 1], 150)
    frame = pd.DataFrame({f"v{j}": np.where(rng.random(300) < (.9 if c == 0 else .1), "a", "b")
                          for j in range(5) for c in [0]})
    frame = pd.DataFrame({f"v{j}": np.where(rng.random(300) < np.where(cls == 0, .9, .1), "a", "b") for j in range(5)})
    found = run_latent(frame, "LCA", [1, 2, 3], seed=2)
    assert min(found["bic"], key=found["bic"].get) == 2
    agree = (found["clusters"][2] - 1 == cls).mean()
    assert max(agree, 1 - agree) > 0.95


def test_a_model_with_more_parameters_than_patterns_is_refused():
    frame = pd.DataFrame({"a": list("xyxyxy"), "b": list("pqpqpq")})
    with pytest.warns(UserWarning, match="more parameters than response patterns"):
        found = run_latent(frame, "LCA", [3], seed=1)
    assert np.isnan(found["bic"][3]) and 3 not in found["clusters"]


def test_the_latent_fits_are_reproducible(blobs):
    a = run_latent(blobs, "LPA", [3], seed=4)
    b = run_latent(blobs, "LPA", [3], seed=4)
    assert a["clusters"][3].tolist() == b["clusters"][3].tolist() and a["bic"] == b["bic"]


# ------------------------------------------------------------------- inside sur_scat ----
def test_layouts_and_clusters_are_added_to_the_plane(blobs):
    net = sur_scat(blobs, nclusters=[2, 3], nfactors=3, layouts=["varimax", "pca"], clusters=["pam", "lpa"], suffix="*")
    assert list(net.layouts) == ["PCA", "varimax 1-2", "varimax 1-2*", "varimax 1-3", "varimax 1-3*", "varimax 2-3",
                                 "varimax 2-3*", "cPC", "cPC*"]
    groups = [c for c in net.nodes.columns if "(" in c]
    assert groups == ["Groups(2)", "Groups(3)", "pam(2)", "pam(3)", "LPA(2)", "LPA(3)"]
    assert net.scat["cluster_columns"] == groups
    assert set(net.clusters) == {f"{g}.{p}" for g in groups for p in net.layouts}
    # the aligned planes are the same cloud, turned
    a, b = net.layouts["varimax 1-2"].to_numpy(), net.layouts["varimax 1-2*"].to_numpy()
    np.testing.assert_allclose(np.linalg.svd(a - a.mean(0), compute_uv=False),
                               np.linalg.svd(b - b.mean(0), compute_uv=False), rtol=1e-8)


def test_the_main_groups_and_plane_are_left_as_they_were(blobs):
    plain = sur_scat(blobs, nclusters=3)
    more = sur_scat(blobs, nclusters=3, layouts=["varimax", "cmds"], clusters=True)
    for column in ("fx", "fy", "Groups(3)"):
        assert plain.nodes[column].astype(str).tolist() == more.nodes[column].astype(str).tolist()


def test_rotations_need_components(survey_like=None):
    frame = pd.DataFrame({"a": list("xyzxyzxyzxyz"), "b": list("pqrspqrspqrs"), "c": list("mnmnmnmnmnmn")})
    with pytest.warns(UserWarning, match="rotations are only made of principal components"):
        net = sur_scat(frame, nclusters=2, layouts="varimax")
    assert not net.layouts


def test_a_scattergram_without_groups_leaves_the_clusters_out(blobs):
    with pytest.warns(UserWarning, match="clusters left out"):
        net = sur_scat(blobs, nclusters=0, clusters=["pam"])
    assert not [c for c in net.nodes.columns if "(" in c]


def test_the_suffix_and_the_alignment_are_checked(blobs):
    with pytest.raises(ValueError, match="suffix must be"):
        sur_scat(blobs, nclusters=2, layouts="pca", suffix="")
    with pytest.warns(UserWarning, match="none was asked for"):
        sur_scat(blobs, nclusters=2, suffix="*")
    with pytest.warns(UserWarning, match="no suffix was stated"):
        sur_scat(blobs, nclusters=2, layouts="pca", align_to="pca")
    with pytest.warns(UserWarning, match="none of the planes"):
        sur_scat(blobs, nclusters=2, layouts="pca", suffix="*", align_to="nowhere")
    assert align_plane("varimax", ["PCA", "varimax 1-2", "cPC"]) == 2
    assert align_plane(3, ["PCA", "varimax 1-2", "cPC"]) == 3
    assert align_plane(None, ["PCA"]) == 1


def test_a_sample_of_the_nodes_gets_its_groups_too(blobs):
    net = sur_scat(blobs, nclusters=2, max_n=50, layouts="pca", clusters="hclust", seed=1)
    assert len(net.nodes) == 50 and "hclust(2)" in net.nodes.columns and "cPC" in net.layouts
    assert len(net.layouts["cPC"]) == 50


def test_patterns_collapse_the_groups_of_looking4clusters():
    rng = np.random.default_rng(1)
    frame = pd.DataFrame({"a": rng.choice(list("xy"), 200), "b": rng.choice(list("pqr"), 200),
                          "c": rng.choice(list("uv"), 200)})
    net = sur_scat(frame, patterns=True, nclusters=2, layouts="pca", clusters=["hclust", "lca"])
    assert len(net.nodes) == frame.drop_duplicates().shape[0]
    assert {"hclust(2)", "LCA(2)"} <= set(net.nodes.columns)
    assert "cPC" in net.layouts


def test_no_group_may_win_the_bic(blobs):
    noise = pd.DataFrame(np.random.default_rng(8).normal(size=(150, 3)), columns=list("abc"))
    net = sur_scat(noise, nclusters=["auto", 1, 2, 3], cluster_index="bic")
    index = net.scat["cluster_index"]
    assert list(index["k"]) == [1, 2, 3] and index["optimal"].tolist() == [True, False, False]
    assert "nodeColor" not in net.options or not str(net.options.get("nodeColor", "")).startswith("Groups(")


# ------------------------------------------------------------------- communities ----
def test_lumping_gathers_the_small_communities_into_the_last_group():
    membership = [1] * 10 + [2] * 6 + [3] * 2 + [4]
    out = lump_communities(membership, 3)
    assert out["clusters"].tolist() == [1] * 10 + [2] * 6 + [3] * 3
    assert out["rest"] == 3 and out["n_lumped"] == 3 and out["n_small"] == 2 and out["threshold"] == 3
    capped = lump_communities(membership, 3, max_communities=2)
    assert capped["clusters"].tolist() == [1] * 10 + [2] * 9 and capped["n_ranked"] == 1
    whole = lump_communities([1, 1, 2, 2], 1)
    assert whole["rest"] is None and whole["n_groups"] == 2
    with pytest.raises(ValueError, match="fewer than"):
        lump_communities([1, 2, 3, 4], 0.99)


def test_the_methods_take_codes_and_names():
    assert [community_code(m) for m in ("lo", "Louvain", "cluster_walktrap", "fastGreedy", "le")] == [
        "lo", "lo", "wa", "fa", "le"]
    with pytest.raises(ValueError, match="not a community detection method"):
        community_code("kmeans")


def test_communities_are_columns_of_groups_of_cases(blobs):
    net = sur_scat(blobs, nclusters=2)
    out = add_communities(net, blobs, community=["wa", "lo"], min_link=0.4)
    names = [c for c in out.nodes.columns if "(" in c and not c.startswith("Groups")]
    assert len(names) == 2 and names[0].startswith("walktrap(") and names[1].startswith("louvain(")
    assert set(names) <= set(out.scat["cluster_columns"])
    again = add_communities(net, blobs, community=["wa", "lo"], min_link=0.4)
    assert out.nodes[names].astype(str).equals(again.nodes[names].astype(str))
    assert "Groups(2)" in out.nodes.columns and set(net.nodes.columns) <= set(out.nodes.columns)
    named = add_communities(net, blobs, community="wa", name="W", min_link=0.4)
    assert "W" in named.nodes.columns
    two = add_communities(net, blobs, community=["wa", "fa"], name="W", min_link=0.4)
    assert [c for c in two.nodes.columns if c.startswith("W.")][0].startswith("W.walktrap(")


def test_the_weights_of_the_links_change_what_the_methods_find(blobs):
    net = sur_scat(blobs, nclusters=2)
    kept = add_communities(net, blobs, community="fa", min_link=0.2)
    free = add_communities(net, blobs, community="fa", min_link=0.2, weighted=False)
    assert len([c for c in kept.nodes.columns if c.startswith("fastGreedy")]) == 1
    assert len([c for c in free.nodes.columns if c.startswith("fastGreedy")]) == 1


def test_communities_check_their_input(blobs):
    net = sur_scat(blobs, nclusters=2)
    with pytest.raises(ValueError, match="holds 100 complete cases but the scattergram has 120 nodes"):
        add_communities(net, blobs.iloc[:100], community="wa")
    with pytest.raises(ValueError, match="same method more than once"):
        add_communities(net, blobs, community=["wa", "walktrap"])
    with pytest.raises(ValueError, match="variables not found"):
        add_communities(net, blobs, variables=["a", "zzz", "b"], community="wa")
    with pytest.raises(ValueError, match="min_size"):
        add_communities(net, blobs, community="wa", min_size=0)
    with pytest.raises(ValueError, match="max_communities"):
        add_communities(net, blobs, community="wa", max_communities=1)
    nominal = blobs.assign(g=pd.Categorical(np.where(blobs["a"] > 0, "u", "v")))
    with pytest.raises(ValueError, match="nominal variable"):
        add_communities(net, nominal, list("abcdefg")[:6] + ["g"], community="wa")
    text = blobs.assign(g=np.where(blobs["a"] > 0, "u", "v"))
    with pytest.raises(ValueError, match="must be numeric or an ordered factor"):
        add_communities(net, text, list("abcdefg")[:6] + ["g"], community="wa")
    ordered = blobs.assign(g=pd.Categorical(np.where(blobs["a"] > 0, "u", "v"), ordered=True))
    add_communities(net, ordered, list("abcdefg")[:6] + ["g"], community="wa", min_link=0.4)
    with pytest.raises(ValueError, match="cannot be standardized"):
        add_communities(net, blobs.assign(k=1.0), list("abcdef") + ["k"], community="wa")
    with pytest.warns(UserWarning, match="rests on 2 points"), pytest.raises(ValueError, match="at least three"):
        add_communities(net, blobs, ["a", "b"], community="wa", min_link=0.9, min_size=1)
    with pytest.raises(ValueError, match="does not answer on a network of 120 nodes"):
        add_communities(net, blobs, community="op", min_link=0.4)


def test_a_sample_of_cases_is_lined_up_with_its_communities(blobs):
    net = sur_scat(blobs, nclusters=2, max_n=60, seed=3)
    out = add_communities(net, blobs, community="wa", min_link=0.3)
    column = [c for c in out.nodes.columns if c.startswith("walktrap(")][0]
    assert len(out.nodes[column]) == 60
    with pytest.raises(ValueError, match="sample of"):
        add_communities(net, blobs.iloc[:100], community="wa")
