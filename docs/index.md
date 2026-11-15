# netcoin documentation

`netcoin` is the Python port of the R package [netCoin](https://CRAN.R-project.org/package=netCoin):
it analyses coincidences (cases that share characteristics, words that appear together, people who
lived in the same years) and draws the result as an interactive network, scattergram, bar chart,
timeline or matrix of pies that you explore in the browser.

| Page | What it is for |
|---|---|
| [User guide](guide.md) | From a table of cases to the networks of a survey, correlations, log-linear and regression models, and a page that gathers them |
| [The scattergram](scattergram.md) | `sur_scat`: planes, groups, rotations, latent models, communities and the comparison of groupings |
| [From the R package to Python](r-to-python.md) | Names, objects, what is not ported and where the results differ from R |
| [Reference](reference.md) | Every function: signature, parameters and what it returns (generated from the docstrings) |

The `README` has the installation and a quick start.

## Which function for what

| I have | I want | Use |
|---|---|---|
| nodes and links | an interactive network | `netCoin` |
| a table of 0/1 columns, or categories in a cell | the incidences | `dichotomize` |
| incidences | counts, similarities, expected values | `coin`, `sim`, `expected_list`, `edge_list` |
| incidences | a complete network | `all_net` |
| a survey | the network of its categories | `sur_coin`, `cob_coin` |
| a survey | its categories as a log-linear / regression network | `log_coin`, `glm_coin` |
| numeric columns | the network of their correlations | `net_corr`, `d_net_corr` |
| a survey | its cases on a plane, in groups | `sur_scat` and the functions that extend it |
| numeric columns | groups of cases and two-dimensional views | `looking4clusters` |
| a data file and a list of analyses | one page with all of them | `caring` |
| several graphs | one page | `multigraphCreate` |
| periods, events, categories, pairs | timeline, bars, pies | `time_coin`, `bar_coin`, `pie_coin` |
| a network | igraph, Pajek, GraphML | `to_igraph`, `from_igraph`, `save_pajek`, `save_graphml` |

## Installing the extras

The core needs `pandas`, `numpy` and `scipy`. The rest is installed by need:

| Extra | Brings | Needed for |
|---|---|---|
| `igraph` | python-igraph | communities, centralities and igraph layouts in `netCoin`; `to_igraph`, `from_igraph`, `save_graphml`; `add_communities` |
| `models` | statsmodels, patsy | `glm_coin` and `caring` with regressions |
| `cluster` | scikit-learn, umap-learn | k-means in `sur_scat`; k-means, t-SNE, NMF and UMAP in `looking4clusters` |
| `files` | openpyxl, pyreadstat | `caring_read_file` of `.xlsx` and `.sav` |

```bash
pip install -e ".[igraph,models,cluster,files]"
```
