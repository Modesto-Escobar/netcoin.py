# From the R package to Python

`netcoin` ports netCoin 2.1.49 for R. Results are checked against R: `tests/test_parity_*.py` compare
the tables, the options and the numbers of almost 400 calls with what R returns (the values are in
`tests/reference/`; R is only needed to regenerate them). This page tells how to read an R script as a
Python one and where the two differ.

## Names

Functions and arguments are written in `snake_case`; the rest of the name is the R one.

| R | Python |
|---|---|
| `netCoin`, `multigraphCreate` | `netCoin`, `multigraphCreate` (kept, as classes of the package) |
| `barCoin`, `timeCoin`, `pieCoin`, `incTime`, `mobileEdges` | `bar_coin`, `time_coin`, `pie_coin`, `inc_time`, `mobile_edges` |
| `dichotomize`, `coin`, `coocur`, `sim`, `edgeList`, `expectedList`, `propCoin`, `asNodes` | `dichotomize`, `coin`, `coocur`, `sim`, `edge_list`, `expected_list`, `prop_coin`, `as_nodes` |
| `allNet`, `surCoin`, `cobCoin`, `netCorr`, `dNetCorr`, `logCoin`, `glmCoin` | `all_net`, `sur_coin`, `cob_coin`, `net_corr`, `d_net_corr`, `log_coin`, `glm_coin` |
| `surScat`, `addAxes`, `addClusters`, `replaceClusters`, `layoutsAlign`, `addCommunities` | `sur_scat`, `add_axes`, `add_clusters`, `replace_clusters`, `layouts_align`, `add_communities` |
| `clusterMatch`, `clusterAgree`, `looking4clusters` | `cluster_match`, `cluster_agree`, `looking4clusters` |
| `caring` (and its helpers) | `caring`, `caring_create_graphs`, `caring_read_file` |
| `layoutMCA`, `layoutPCA`, `layoutFact`, `layoutCircle`, `layoutGrid` | `layout_mca`, `layout_pca`, `layout_fact`, `layout_circle`, `layout_grid` |
| `toIgraph`, `fromIgraph`, `calCentr`, `savePajek`, `saveGhml` | `to_igraph`, `from_igraph`, `cal_centr`, `save_pajek`, `save_graphml` |

Arguments: `minL`, `maxL`, `linkBipolar`, `valueDicho`, `nclusters`, `clusterIndex`, `maxN`... are
`min_link`, `max_link`, `link_bipolar`, `value_dicho`, `nclusters`, `cluster_index`, `max_n`... In
`netCoin(...)` the camelCase names of R and of earlier versions (`showArrows`, `nodeFilter`...) still
work, with a warning. `caring` and `caring_create_graphs` take the keys of the R function, and of its
graphical interface, as they are.

Values: a `data.frame` is a `DataFrame`; a factor is a categorical column (an ordered factor is an *ordered*
one); `NULL` is `None`; a named list is a `dict`. Positions that a user gives are **1-based**, as in R
(`Source` and `Target` as numbers, `which`, `to`, the case-to-pattern map kept in `net.scat`).

Objects: every network function returns a `netCoin` whose tables are `net.nodes` and `net.links`, the options
the page reads are `net.options`, and the planes of a scattergram `net.layouts`. The attributes R
keeps with `attr(x, "...")` live in `net.scat` (`case_to_pattern`, `cluster_columns`, `cluster_index`...); the
centroids of the groups are `net.clusters` and `net.cluster_means`. `net.add(...)` makes a copy with
some arguments changed.

## What is not ported

The galleries (`gallery`, `netGallery`, `netExhibit`) and `pathCoin` (which needs `lavaan`) are not
ported. Everything else of the list above is.

## Libraries behind it

| R | Python |
|---|---|
| igraph | python-igraph (extra `igraph`): communities, centralities, layouts, `to_igraph` |
| `glm`, `nnet::multinom`, `loglin` | statsmodels, patsy (extra `models`); the log-linear fit (IPF) and the multinomial are written in the package |
| `marginaleffects` | written in the package (delta method, against the weighted mean of every variable) |
| `cluster::pam`, `hclust`, `kmeans`, `Rtsne`, `NMF`, `umap` | PAM and the distances are written in the package; hclust is scipy's; k-means, t-SNE and NMF are scikit-learn's; UMAP is umap-learn (extra `cluster`) |
| mclust, poLCA | written in the package (an EM for each; see below) |
| reading `.sav`, `.dta` and `.xlsx` (`caring`) | pyreadstat, openpyxl (extra `files`) |

## Where Python differs on purpose

A suspected bug of R is **not** copied, except when copying it is what makes the objects equal. The
complete list, with the R behaviour of every case, is in `tests/reference/README.md`; these are the ones
you can notice:

- The statistic of the correlations of metric variables in `sur_coin` is the usual
  `t = r·sqrt((n-2)/(1-r²))`; R divides by `1-r`, so `p(Z)` of those links differs.
- Weights given as a vector with incomplete cases are matched with their cases (R shifts them), in
  `sur_coin` and `log_coin`.
- A node table given with a layout `"mc"` or `"pc"` (`sur_coin`, `all_net`, `net_corr`) is aligned by name;
  R gives each node the coordinates of another one when it re-orders them.
- `sim(..., "Haberman")`, `edge_list(..., "Fisher")`, `dichotomize` with empty terms, `coin(total=TRUE)`,
  `expected_list(confidence=True)` and `all_net(pairwise=True)` with missing values work as intended, where R
  fails or repeats a value.
- `to_igraph` carries the width, weight, colour and text of the links; R's `toIgraph` loses them. For
  that reason `add_communities` weighs the links by their correlation unless you say `weighted=False` (which
  gives the communities R finds).
- `caring` recycles `valueDicho` for every dichotomy and reads text columns of numbers as numbers, and
  `glm_coin` leaves the aliased terms of collinear models out of the fit instead of stopping.

## Where the results cannot be the same number by number

- **Random numbers.** R's and numpy's generators differ, so the samples that `max_n` draws, the starts of
  k-means, t-SNE, UMAP and NMF, and the random layouts (`fruchterman_reingold`, `drl`, `random`...)
  differ. Where there are clear groups the partition is the same; where there are not, the solutions are
  different local optima (the Python k-means keeps the best of several starts, refined by optimal
  transfers, as R's Hartigan-Wong does).
- **Ties.** hclust and PAM on data with many equal distances, and the optimal assignment that matches
  categories (`cluster_match`), may break ties in another order; the result reports them in `attrs["ties"]`.
- **igraph's version.** R's igraph and python-igraph have different C cores: `kamada_kawai` and `sugiyama`
  give other coordinates, and the community methods that start at random (`louvain`, `label_prop`,
  `spinglass`) cannot be compared one to one. `gem` does not exist in python-igraph.
- **Latent models.** `clusters="lpa"` and `"lca"` are fitted by an EM written in the package (equal diagonal
  variances for profiles, independent categories inside each class), not by mclust and poLCA. The BIC
  agrees to 0.2 % and the groups wherever the model is well identified.
- **The sign of an axis** of a factorial analysis is arbitrary, and so is the numbering of the groups that
  depends on it (they are numbered from left to right on the first axis).
