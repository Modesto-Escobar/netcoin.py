# netcoin: Interactive Analytic Networks for Coincidences in Python

<!-- badges: start -->
[![License: GPL (>= 2)](https://img.shields.io/badge/license-GPL%20(%3E%3D%202)-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
[![Python](https://img.shields.io/badge/python-%3E%3D3.9-blue.svg)](https://www.python.org/)
<!-- badges: end -->

**netcoin** is the Python port of the R package
[**netCoin**](https://CRAN.R-project.org/package=netCoin). It combines the
data-analysis capabilities of Python (pandas, numpy, scipy) with the
interactive visualization libraries of JavaScript (D3.js) to create networks
of coincidences, co-occurrences and similarities that can be explored in a
browser or embedded as self-contained HTML files.

> The repository is named **`netcoin.py`** to distinguish it from the original
> R package; the importable package itself is just `netcoin`.

## Mission

This project's aim is to integrate traditional statistical techniques with
machine-learning and social-network-analysis tools in order to produce visual
and interactive displays of large datasets. The objectives are:

1. Efficiently combine different statistical techniques by integrating them
   under the study of the coincidence of people, objects, events or
   characteristics across multiple scenarios.
2. Provide open-source software that, under the premise of network
   coincidence analysis, generates interactive graphics for both exploratory
   and confirmatory analysis of vast quantities of information.
3. Apply the above to large databases in fields such as:
   - survey data combined with administrative records;
   - networks created by social-media users and reproduced through their
     messages;
   - semantic maps drawn from scientific abstracts over long time spans;
   - large biographical databases of leading figures in philosophy, science,
     the social sciences and the arts, together with their major works.

## Installation

Install the development version directly from GitHub:

```bash
pip install git+https://github.com/Modesto-Escobar/netcoin.py.git
```

Or clone the repository and install in editable mode for local development:

```bash
git clone https://github.com/Modesto-Escobar/netcoin.py.git
cd netcoin.py
pip install -e ".[dev]"
python -m pytest
```

Runtime dependencies (installed automatically): `pandas`, `numpy`, `scipy`.
Python 3.9 or newer is required.

## Quick start

### Building a network from nodes and links

```python
from netcoin import netCoin

nodes = {
    "name":     ["A", "B", "C", "D"],
    "category": ["x", "x", "y", "y"],
    "size":     [10, 20, 15, 25],
}

# Numeric Source/Target values are 1-based positions in ``nodes``
# (as in the R package); string values are matched against ``name``.
links = {
    "Source": [1, 2, 1, 3],
    "Target": [2, 3, 4, 4],
    "weight": [0.5, 0.8, 0.3, 0.9],
}

net = netCoin(
    nodes=nodes,
    links=links,
    name="name",         # node identifier
    color="category",    # map a node attribute to colour
    size="size",         # map a node attribute to size
    main="My Network",
    language="en",
)
net.plot(directory="output/network")   # writes an interactive HTML file
```

### From raw data to a coincidence network

```python
import pandas as pd
from netcoin import dichotomize, coin, sim, edge_list, netCoin

data = pd.DataFrame({"keywords": ["red;blue", "blue;green", "red", "red;green"]})

incidences = dichotomize(data, "keywords", sep=";", nas=None)  # 0/1 columns
cooc = coin(incidences)          # a Coin object: coincidences + n
print(cooc)                      # n= 4, lower triangle

similarities = sim(cooc, "Jaccard")              # square DataFrame
edges = edge_list(cooc, ["Jaccard", "Haberman"], criteria="Jaccard", min=0)
```

`coin()` returns a `Coin` object (the coincidence matrix in `.matrix` and the
number of cases in `.n`), which `sim()` and `edge_list()` take as input. Both
also accept the table of incidences directly. Procedure names follow the R
package: `"Jaccard"`, `"Haberman"`, `"Z"`, `"confidence"`... or their short
forms (`"J"`, `"H"`, `"Ham"`, `"CC"`, `"Hyp"`).

### A complete network in one call

```python
from netcoin import all_net

net = all_net(incidences, procedures=["Jaccard", "Haberman"], lwidth="Jaccard",
              layout="pc", main="Keywords")      # counts, nodes, links and graph
net = net.add(color="%")                          # a copy with some arguments changed
net.plot()                                        # opens it in the browser
```

### From a survey

`sur_coin` takes a table of cases (one column per question) and makes every category a node
(`"question:answer"`); the links tell which categories go together more or less than expected.
Quantitative variables (`metric`) become nodes too and are correlated with the categories.
`cob_coin` keeps the significant links between different questions and puts the nodes on a
circle (a CoWeb graph).

```python
from netcoin import sur_coin, cob_coin

net = sur_coin(survey, ["sex", "age", "works"], metric="income",
               lwidth="Haberman", lcolor="Haberman", link_bipolar=True, color="variable",
               commonlabel="sex")
web = cob_coin(survey, ["sex", "age", "works"], degree=90, significance=0.05)
web.plot()
```

### From quantitative variables

`net_corr` makes every numeric column a node (with its mean, standard deviation and range) and
every pair a link with the correlation (`value`), its contrast (`z`, `p`) and, with `pairwise`,
the number of cases. `method` is `"pearson"`, `"spearman"`, `"kendall"`, `"partial"` or
`"dcor"` (distance correlation), alone or several at once; weights and missing values are
supported, and small samples use the exact tests. `d_net_corr` shows how the network breaks
into components as the threshold of the correlation rises.

```python
from netcoin import net_corr, d_net_corr

net = net_corr(cars, method=["pearson", "spearman"], layout="pc", main="Correlations")
evolving = d_net_corr(cars, factorial="pc", sequence=[0.3, 0.5, 0.7, 0.9],
                      directory="corr_evolving")
```

### Log-linear associations

`log_coin` fits a log-linear model to the table of several categorical variables and draws the
categories (nodes sized by their main effect) joined by the interactions that are significant.
The note of the graph reports the fit of the model (L2, degrees of freedom and p).

```python
from netcoin import log_coin

net = log_coin(survey, ["gender", "age", "works"], pmax=0.05, frequency=True, percentage=True)
```

### Drawing options

`netCoin` has the arguments of the R function, in `snake_case`; the tables are DataFrames
(`net.nodes`, `net.links`) and `net.options` is what the page reads.

```python
net = netCoin(nodes, links, name="name", label="label", size="size",
              color="group", shape="kind",           # a column, one value, or one per node
              lwidth="weight", lcolor="weight",
              layout="fruchterman_reingold",         # or a matrix, or several planes in a dict
              community="walktrap",                  # detected with igraph (pip install netcoin[igraph])
              node_filter="size > 2", degree_filter=2,
              centrality=["closeness", "betweenness"],
              controls=[1, 2, 3, 4, 5], show_arrows=True, main="My network")
net.plot()
```

The names of earlier versions (`showArrows`, `defaultColor`, `nodeFilter`, `fromMatrix`...) still work
but warn: use `show_arrows`, `default_color`, `node_filter`, `from_matrix`... `orderA`/`orderD` are
now `sort`/`decreasing`.

With python-igraph you can also go back and forth and save for other programs:

```python
from netcoin import to_igraph, from_igraph, save_pajek, save_graphml, layout_circle, layout_grid

graph = to_igraph(net)                      # an igraph.Graph with the node and link attributes
net2 = from_igraph(graph)                   # and back
save_pajek(net, "network.net", arcs=["weight"], partitions=["group"])
save_graphml(net, "network.graphml")
coordinates = layout_circle(net.nodes, name="name")      # also layout_grid("a,b.c.d")
```

### Several graphs in one page

```python
from netcoin import multigraphCreate

multigraphCreate([net1, net2, net3], names=["2019", "2020", "2021"],
                 mode="default")                    # a selector between graphs
multigraphCreate([net1, net2], mode="parallel")     # side by side
multigraphCreate([net1, net2, net3], mode="frame",  # one network that evolves
                 directory="evolution", speed=30, loop=True)
```

An existing folder is only replaced if it holds a graph written by this package.

### Bars, timelines and pies

```python
from netcoin import bar_coin, time_coin, pie_coin, coin

# coincidences of every category with the others (survey data, weighted)
bar = bar_coin(survey, ["Gender", "Age"], dichotomies="_none", weight="cweight",
               expected=True, confidence=True, significance=True)
bar.plot()

# periods (and events on them) on a time axis
time_coin(people, start="birth", end="death", group="school",
          events=works, event_names="works", event_period="name", event_time="date").plot()

# a matrix of pies, one per pair of variables
pie_coin(coin(incidences), main="Coincidences").plot()
```

Each returns an object (`BarCoin`, `TimeCoin`, `PieCoin`) with its tables and options, which
can be written with `directory=` or `.plot()`, and put in a `multigraphCreate`.

## Main features

| Area                              | Key API                                    |
|-----------------------------------|--------------------------------------------|
| Interactive coincidence networks  | `netCoin(...)`, `netCoin.plot(...)`        |
| Build a network from a matrix     | `netCoin.fromMatrix(...)`                  |
| Bar charts of coincidences        | `bar_coin(...)`                            |
| Timelines of periods and events   | `time_coin(...)`, `inc_time(...)`, `mobile_edges(...)` |
| Pie charts of coincidences        | `pie_coin(...)`                            |
| Layouts, igraph, Pajek, GraphML   | `layout_circle`, `layout_grid`, `to_igraph`, `from_igraph`, `save_pajek`, `save_graphml` |
| Multi-graph visualisations        | `multigraphCreate(...)` (default, parallel, frame) |
| Coincidences (`Coin` object)      | `coin(...)`                                |
| Similarity measures (Jaccard, …)  | `sim(...)`                                 |
| Observed vs expected coincidences | `expected_list(...)`                       |
| Edge lists                        | `edge_list(...)`                           |
| Dichotomising categorical data    | `dichotomize(...)`                         |
| Counts and percentages            | `coocur(...)`, `prop_coin(...)`, `as_nodes(...)` |
| Complete network from incidences  | `all_net(...)`, `netCoin.add(...)`         |
| Networks from a survey, CoWeb     | `sur_coin(...)`, `cob_coin(...)`           |
| Log-linear associations           | `log_coin(...)`                            |
| Networks of correlations          | `net_corr(...)`, `d_net_corr(...)`, `layout_fact(...)` |
| Layouts from coincidences         | `layout_mca(...)`, `layout_pca(...)`       |

Visualisations are powered by D3.js. All required web assets (`d3.min.js`,
templates, stylesheets and language files for English, Spanish and Catalan)
ship with the package under `netcoin/www/` and are copied to the output
directory automatically. They are the engine of the R package
rD3plot (version recorded in
`netcoin/www/SOURCE.json`); to update them run
`python tools/sync_www.py --source <path to rD3plot>`.

## Documentation and examples

Runnable examples live in [`python_examples/`](python_examples/) and cover
several real datasets (Galápagos finches, Italian Renaissance families,
artworks, etc.). Run, for instance:

```bash
python python_examples/python_net_Coin_example.py
```

The generated `index.html` can be opened directly in a browser.

For the conceptual background, the methodology and the full feature set, see
the R package documentation at
<https://modesto-escobar.github.io/netCoin/>.

## Citation

If you use **netcoin** in academic work, please cite the original R package
until a dedicated Python citation is available:

> Escobar, M., Barrios, D., Prieto, C., Martínez-Uribe, L.,
> Cabrera-Álvarez, P., & Calvo-López, C. *netCoin: Interactive Analytic
> Networks*. <https://CRAN.R-project.org/package=netCoin>

### Further reading

- Escobar, M. (2015). Studying Coincidences with Network Analysis and Other
  Multivariate Tools. *The Stata Journal*, 15(4), 1118–1156.
- Escobar, M., & Martínez-Uribe, L. (2020). Network Coincidence Analysis:
  The Netcoin R Package. *Journal of Statistical Software*, 93(11), 1–31.
- Escobar, M., & Calvo, C. (2024). *Redes Analíticas. Grafos de Coincidencia
  y Regresión*. Madrid: CIS.

## Authors

- **Modesto Escobar** (Universidad de Salamanca) — creator & maintainer
  [[ORCID]](https://orcid.org/0000-0003-2072-6071)
- David Barrios (Universidad de Salamanca)
- Carlos Prieto (Universidad de Salamanca)
  [[ORCID]](https://orcid.org/0000-0003-2064-4842)
- José Antonio Robles Querencia (Universidad de Salamanca)
- Emilio Portela Navarro (Universidad de Salamanca)

## Acknowledgments

This work has been supported by grants **CSO2013-49278-EXP**,
**PGC2018-093755-B100**, **PDC2022-133355-100** and **PID2023-147358NB-100**,
funded by MICIU/AEI/10.13039/501100011033 and by the European Union
NextGenerationEU/PRTR programme.

## Bug reports and contributions

Please file issues and pull requests at:

<https://github.com/Modesto-Escobar/netcoin.py/issues>

## License

GPL-2 | GPL-3
