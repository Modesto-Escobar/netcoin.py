# Reference values from the R package netCoin

`make_reference.R` runs netCoin (R) on the datasets shipped with the package
plus a few synthetic tables and stores the results as CSV. The Python tests in
`tests/test_parity_*.py` compare against them, so **R is not needed to run the
tests**, only to regenerate the references:

```bash
Rscript tests/reference/make_reference.R        # from the repository root
```

* `inputs/`      data given to both implementations
* `expected/`    what netCoin returns (`coin_*`, `sim_*`, `edge_*`, `dicho_*`)
* `manifest.json` R / netCoin / rD3plot versions, collation locale and the R call of each case

Frequency ties in `dichotomize` follow the **locale collation** of the R
session (`collate` in the manifest), so regenerate under the same locale.

## Known gaps

Expected failures are declared in each test file (`GAPS`, `known_gap`) and are
*strict*: when a gap is closed the test reports `XPASS(strict)` and the entry
must be deleted. Run `pytest -rx` to list them with their reasons.

## Behaviour of R 2.1.49 that looks like a bug (not copied)

| Where | What R does |
|---|---|
| `sim(..., "H")` (Haberman) | `sqrt(N)[is.na(x)]` fills only the first `NA` with `sqrt(N)`; the others stay `NaN`. Python fills all of them. |
| `expectedList(confidence = TRUE)` | Fails (`cbind` with 78 and 64 rows): it calls `edgeList(<matrix>, "shape", min = -Inf)` and `edgeList` rewrites `min = -Inf` to `1`. |
| `dichotomize` empty tokens | `"a;;b"`, or `sep = "C"` with a `""` value, keeps an empty term that ends up as an unnamed column (`V3`). Python drops empty tokens. |
| `edgeList(..., criteria = "Fisher")` | The p-value cap (`max = .5`) is applied for `"Z"` and `"Hyp"` but not for `"Fisher"`. |
| `sim(coin, "L")` alone | Returns a single matrix instead of `conf.L`, `confidence`, `conf.U`. |
| `coin(total = TRUE)` | The `Total` variable is listed twice (it is added and also passes the frequency filter). Python lists it once. |
| `coin(weight = w)` | Overwrites a user-given `maximum` with `sum(w)`. Python honours `maximum` and defaults to `sum(w)`. |
| `allNet(..., pairwise = TRUE)` with missing values | Fails: `colSums(incidences) > 0` is `NA` and is used as a column index. Python ignores the `NA` in the sums. |
| `allNet(..., layout = "mc")` | The layout is computed on all the variables, so it does not line up with the nodes if `minimum`, `maximum` or `sort` change them. Python uses the variables that are nodes. |
| `netCoin(node_filter / degreeFilter)` | The helper `hideLinks()` changes a local copy of the links and has no effect: the links of hidden nodes are never marked `hidden` in the object (the page hides them). Python does the same so that both objects are equal. |
| `cluster_edge_betweenness` | The result is the unweighted one even if a `weight` attribute exists. Python follows it. |
| `surCoin(metric = )`, column `p(Z)` | The t statistic of a correlation is computed as `r / sqrt((1 - r) / (n - 2))`; the usual one is `r * sqrt((n - 2) / (1 - r^2))`. Python uses the usual one, so the `p(Z)` of the links of a metric variable differs from R (the tests leave those cells out and check them against a direct calculation). `n` is `nrow(incidences)` in both. |
| `surCoin(weight = <vector>)` with incomplete cases | `na.omit(cbind(data, weight))` drops the rows with a missing value but the vector of weights keeps its length, so the weights no longer match the cases. Python drops the weights with their cases (like `weight = "<column>"`, which R does handle). |
| `surCoin(nodes = , layout = "mc")` (or `"pc"`) | `merge()` sorts the nodes alphabetically but the coordinates stay in the order of the incidences, so each node gets the coordinates of another one. Python aligns them by name (`nodes_mc` is a strict xfail). `allNet` has the same flaw. |
| `surCoin(nodes = , metric = )` | The metric variables listed in `nodes` are reported as "not present in the data frame". Python does not mention them. |
| `surCoin(metric = )` with an undefined correlation | A constant variable gives `NaN`, which turns into rows of `NA` in the links. Python leaves those links out. |
| `cobCoin(dir = )` | The page is written by the inner `surCoin` before the links are filtered. Python writes the final network. |
| `netCorr(nodes = , layout = "pc")` | The coordinates of `layoutPCA` stay in the order of the correlation matrix while the nodes are re-ordered by `merge()`, so each node gets the coordinates of another one. Python aligns them by name (`layout_pc_nodes` is a strict xfail). |
| `toIgraph()` (rD3plot) | For the links it looks for the options `lwidth`, `lweight`, `lcolor` and `ltext`, which are called `linkWidth`... in the object, so the width, weight, colour and text never reach the igraph graph. Python maps them. In `d_netCorr` this makes R's `wdegree` equal to `degree`; the tests leave that column out. |
| `d_netCorr` eigenvector centrality | For frames with several components the eigenvector is not unique (ARPACK returns different values); the tests leave that column out. |
| group numbers of communities | The numbers `G.1`, `G.2`... come from the algorithm's internals; the groups are the same but may be numbered differently (the tests compare the partition). |

### Replicated on purpose

`surCoin` turns a missing value of a categorical variable into a case with none of its categories
(a row of zeros), also with `pairwise = TRUE`, so the cases are not dropped and `pairwise` only
changes the correlations of the metric variables. With a single variable and a dichotomy
the categories get no `variable:` prefix (`Female`, `Male`) and `variable` repeats the name. With
metric variables and `layout = "mc"`/`"pc"` the layout is refused with a warning (there are no
coordinates for the metric nodes).

## Differences that come from the igraph version

R uses igraph 2.3.4 and python-igraph 1.0.0 has a newer C core, so `kamada_kawai` and `sugiyama`
return other coordinates (strict xfails). Layouts that start at random (`fruchterman_reingold`,
`drl`, `lgl`, `random`, `force_atlas2`...) and communities such as `louvain`, `label_prop` or
`spinglass` cannot be compared number by number. `gem` is not in python-igraph.

