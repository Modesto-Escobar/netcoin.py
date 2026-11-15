"""``looking4clusters``: the object, the methods that start at random and the checks."""
import numpy as np
import pandas as pd
import pytest

from netcoin import looking4clusters
from netcoin.clusters import (Looking4Clusters, as_factor, calinski_harabasz, cutree, n_components, pam,
                              run_kmeans, run_nmf, run_pca)
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform


@pytest.fixture
def blobs():
    rng = np.random.default_rng(0)
    centres = np.array([[0, 0, 0], [6, 6, 0], [0, 6, 6]])
    return pd.DataFrame(np.vstack([c + rng.normal(0, 0.6, (25, 3)) for c in centres]), columns=["a", "b", "c"])


def test_as_factor_orders_numbers_by_value_and_text_alphabetically():
    assert list(as_factor([10, 2, 2, 1]).categories) == ["1", "2", "10"]
    assert list(as_factor(["pam_2", "pam_10", "pam_1"]).categories) == ["pam_1", "pam_10", "pam_2"]


def test_the_object_describes_itself(blobs):
    text = repr(looking4clusters(blobs, selectedk=3, seed=1))
    assert text.startswith("An object of class looking4clusters\n3 variables across 75 samples")
    assert "clusters added: kmeans, pam, hclust" in text and "dimensional reductions added" in text


def test_a_second_column_for_the_same_method_makes_a_table(blobs):
    obj = Looking4Clusters(blobs.to_numpy(), [f"s{i}" for i in range(75)], list(blobs.columns))
    obj.add_cluster(np.arange(75) % 2, "my groups")
    assert list(obj.clusters) == ["my_groups"] and isinstance(obj.clusters["my_groups"], pd.Series)
    obj.add_cluster(np.arange(75) % 3, "my groups", optim_cluster=True)
    frame = obj.clusters["my_groups"]
    assert list(frame.columns) == ["levels_2", "levels_3"] and obj.optim["my_groups"] == 3


def test_groups_of_your_own_are_kept_apart(blobs):
    obj = looking4clusters(blobs, groups=["x"] * 40 + ["y"] * 35, running_all=False)
    assert obj.options["myGroups"] == "cluster_1" and list(obj.clusters) == ["cluster_1"]
    assert obj.reductions == {} and list(obj.clusters["cluster_1"].cat.categories) == ["x", "y"]


def test_wrong_lengths_are_refused(blobs):
    obj = looking4clusters(blobs, running_all=False)
    with pytest.raises(ValueError, match="one per sample"):
        obj.add_cluster([1, 2, 3])
    with pytest.raises(ValueError, match="one row per sample"):
        obj.add_reduction(np.zeros((3, 2)))
    with pytest.raises(ValueError, match="matrix like"):
        looking4clusters([1, 2, 3])


def test_samples_are_named_after_the_index(blobs):
    assert looking4clusters(blobs, running_all=False).samples[:2] == ["sample_1", "sample_2"]
    named = blobs.set_index(pd.Index([f"case{i}" for i in range(75)]))
    assert looking4clusters(named, running_all=False).samples[0] == "case0"


# ------------------------------------------------------------ the methods ----

def test_the_best_number_of_groups_is_found(blobs):
    obj = looking4clusters(blobs, seed=0)
    assert obj.optim["pam"] == 3 and obj.optim["hclust"] == 3 and obj.optim["kmeans"] == 3
    for method in ("kmeans", "pam", "hclust"):
        frame = obj.clusters[method]
        assert list(frame.columns) == [f"levels_{k}" for k in range(2, 11)]
        assert len(frame["levels_3"].cat.categories) == 3


def test_selectedk_makes_a_range_around_it(blobs):
    obj = looking4clusters(blobs, selectedk=4, seed=0)
    assert list(obj.clusters["pam"].columns) == [f"levels_{k}" for k in range(2, 10)]
    assert obj.optim["pam"] == 4 and obj.optim["hclust"] == 4


def test_the_range_of_groups_is_limited_by_the_cases():
    rng = np.random.default_rng(1)
    small = pd.DataFrame(rng.normal(size=(6, 3)))
    obj = looking4clusters(small, seed=0)
    assert list(obj.clusters["pam"].columns) == ["levels_2", "levels_3", "levels_4", "levels_5"]


def test_kmeans_with_a_seed_can_be_repeated(blobs):
    a = run_kmeans(Looking4Clusters(blobs.to_numpy(), list(map(str, range(75))), ["a", "b", "c"]), seed=5)
    b = run_kmeans(Looking4Clusters(blobs.to_numpy(), list(map(str, range(75))), ["a", "b", "c"]), seed=5)
    pd.testing.assert_frame_equal(a.clusters["kmeans"], b.clusters["kmeans"])


def test_pam_finds_the_medoids_of_clear_groups(blobs):
    dmat = squareform(pdist(blobs.to_numpy()))
    medoids, labels = pam(dmat, 3)
    assert sorted(m // 25 for m in medoids) == [0, 1, 2]
    assert [len(set(labels[i * 25:(i + 1) * 25])) for i in range(3)] == [1, 1, 1]


def test_cutree_follows_the_order_of_the_merges():
    x = np.array([[0.0], [1.0], [5.0], [6.0], [20.0]])
    tree = linkage(pdist(x), "average")
    assert cutree(tree, 3).tolist() == [1, 1, 2, 2, 3] and cutree(tree, 1).tolist() == [1] * 5
    assert cutree(tree, 5).tolist() == [1, 2, 3, 4, 5]


def test_calinski_harabasz_prefers_the_true_groups(blobs):
    condensed = pdist(blobs.to_numpy())
    true = np.repeat([1, 2, 3], 25)
    assert calinski_harabasz(condensed, true) > calinski_harabasz(condensed, np.arange(75) % 3 + 1)


def test_pca_explains_the_variance(blobs):
    obj = run_pca(Looking4Clusters(blobs.to_numpy(), list(map(str, range(75))), ["a", "b", "c"]))
    assert sum(obj.pctvar["pca"]) <= 100 and obj.pctvar["pca"][0] > obj.pctvar["pca"][1]
    assert list(obj.reductions["pca"].columns) == ["PC1", "PC2"]


def test_the_reductions_are_added(blobs):
    pytest.importorskip("sklearn")
    obj = looking4clusters(blobs.abs(), seed=0)
    assert {"pca", "tsne", "cmds", "nmf"} <= set(obj.reductions)
    assert all(frame.shape == (75, 2) for frame in obj.reductions.values())


def test_nmf_refuses_negative_data_with_a_warning(blobs):
    with pytest.warns(UserWarning, match="negative entries"):
        run_nmf(Looking4Clusters(blobs.to_numpy() - 3, list(map(str, range(75))), ["a", "b", "c"]))


def test_umap_is_optional():
    umap = pytest.importorskip("umap")
    rng = np.random.default_rng(2)
    obj = looking4clusters(pd.DataFrame(np.abs(rng.normal(size=(40, 5)))), running_all=True, seed=0)
    assert "umap" in obj.reductions


# --------------------------------------------------------------- components ----

def test_components_need_enough_cases_and_variables():
    with pytest.raises(ValueError, match="five samples"):
        looking4clusters(pd.DataFrame(np.ones((4, 6))), components=True)
    with pytest.raises(ValueError, match="five variables"):
        looking4clusters(pd.DataFrame(np.ones((6, 4))), components=True)


def test_components_replace_the_variables():
    rng = np.random.default_rng(3)
    counts = pd.DataFrame(rng.poisson(3, (30, 8)))
    obj = looking4clusters(counts, components=True, running_all=False)
    assert obj.variables == [f"PC{i + 1}" for i in range(obj.data.shape[1])] and obj.data.shape[0] == 30


def test_n_components_keeps_four_when_the_drop_is_sharp():
    assert n_components([0.5, 0.2, 0.1, 0.08, 0.05, 0.03, 0.02, 0.01, 0.005, 0.005]) == 4
