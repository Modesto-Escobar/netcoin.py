"""``looking4clusters`` (the deterministic parts) against R."""
import json
import warnings

import numpy as np
import pandas as pd
import pytest

from netcoin.clusters import (Looking4Clusters, _create, run_mds, run_pam_hclust, run_pca)
from conftest import EXPECTED, read_expected, read_input

BLOBS = read_input("l4c_blobs")
COUNTS = read_input("l4c_counts")

# R's pam and hclust runs: id -> (data, distance, agglomeration, selectedk)
RUNS = {
    "euclid_complete": (BLOBS, "euclidean", "complete", None),
    "manhattan_ward": (BLOBS, "manhattan", "ward.D", None),
    "euclid_ward2": (BLOBS, "euclidean", "ward.D2", None),
    "maximum_single": (BLOBS, "maximum", "single", None),
    "euclid_mcquitty": (BLOBS, "euclidean", "mcquitty", None),
    "euclid_average_k4": (BLOBS, "euclidean", "average", 4),
    "euclid_centroid": (BLOBS, "euclidean", "centroid", None),
    "euclid_median": (BLOBS, "euclidean", "median", None),
    "canberra_average": (COUNTS, "canberra", "average", None),
    "binary_complete": (COUNTS, "binary", "complete", None),
    "counts_ward2": (COUNTS, "euclidean", "ward.D2", 3),
    "canberra_average_base": (COUNTS, "canberra", "average", None),
}
# parallelDist, which R uses when it is installed, multiplies a Canberra distance with terms 0/0 by
# (n + 1) / n instead of scaling it to all the variables as dist() does: the "_base" runs are
# made with dist(), as python does
PARALLEL_DIST = {"canberra_average": "pam and hclust of R use parallelDist, whose Canberra distance is "
                                     "scaled differently from dist()",
                 "binary_complete": "binary distances of 12 variables take 9 values: ties decide the groups"}


def _run(case):
    data, distance, agglomeration, k = RUNS[case]
    return run_pam_hclust(_create(data, False), distance, agglomeration, k)


KNOWN_OK = {("pam", "canberra_average")}


@pytest.mark.parametrize("case", RUNS)
@pytest.mark.parametrize("method", ["pam", "hclust"])
def test_clusters_match_r(case, method, request):
    if case in PARALLEL_DIST and (method, case) not in KNOWN_OK:
        pytest.skip("depends on parallelDist (see PARALLEL_DIST)")
    obj = _run(case)
    got = obj.clusters[method].astype(object).astype(str).reset_index(drop=True)
    want = read_expected(f"l4c_{case}_{method}", strings=True)
    assert list(got.columns) == list(want.columns)
    pd.testing.assert_frame_equal(got, want, check_dtype=False)
    optim = json.loads((EXPECTED / f"l4c_{case}_optim.json").read_text(encoding="utf-8"))
    assert obj.optim[method] == optim[method]


@pytest.mark.parametrize("distance", ["canberra", "binary", "euclidean", "maximum", "manhattan", "minkowski"])
def test_dissimilarities_match_r_dist(distance):
    from netcoin.clusters import dissimilarity
    want = read_expected(f"l4c_dist_{distance}")["d"].to_numpy()
    np.testing.assert_allclose(dissimilarity(COUNTS.to_numpy(), distance), want, rtol=1e-9, atol=1e-12)


def test_mds_of_canberra_distances_matches_base_r():
    obj = run_mds(_create(COUNTS, False))
    want = read_expected("l4c_counts_cmds_base")
    for j in range(2):
        a, e = obj.reductions["cmds"].iloc[:, j].to_numpy(), want.iloc[:, j].to_numpy()
        sign = 1 if np.abs(a - e).sum() <= np.abs(a + e).sum() else -1
        np.testing.assert_allclose(sign * a, e, rtol=1e-6, atol=1e-8)
    pct = json.loads((EXPECTED / "l4c_counts_cmds_base_pctvar.json").read_text(encoding="utf-8"))
    np.testing.assert_allclose(obj.pctvar["cmds"], pct["cmds"], rtol=1e-7)


@pytest.mark.parametrize("name,data", [("blobs", BLOBS)])
def test_pca_and_mds_match_r(name, data):
    obj = run_mds(run_pca(_create(data, False)))
    pct = json.loads((EXPECTED / f"l4c_{name}_pctvar.json").read_text(encoding="utf-8"))
    for key, frame in (("pca", obj.reductions["pca"]), ("cmds", obj.reductions["cmds"])):
        want = read_expected(f"l4c_{name}_{key}")
        assert list(frame.columns) == list(want.columns)
        for j in range(2):                       # the sign of an axis is arbitrary
            a, e = frame.iloc[:, j].to_numpy(), want.iloc[:, j].to_numpy()
            sign = 1 if np.abs(a - e).sum() <= np.abs(a + e).sum() else -1
            np.testing.assert_allclose(sign * a, e, rtol=1e-6, atol=1e-8, err_msg=f"{name}/{key}/{j}")
        np.testing.assert_allclose(obj.pctvar[key], pct[key], rtol=1e-7)


def test_components_match_r():
    obj = _create(COUNTS, True)
    want = read_expected("l4c_components")
    assert obj.data.shape == want.shape
    for j in range(want.shape[1]):
        a, e = obj.data[:, j], want.iloc[:, j].to_numpy()
        sign = 1 if np.abs(a - e).sum() <= np.abs(a + e).sum() else -1
        np.testing.assert_allclose(sign * a, e, rtol=1e-6, atol=1e-8)
