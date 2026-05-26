# netcoin: Interactive Analytic Networks for Coincidences in Python

<!-- badges: start -->
[![License: GPL (>= 2)](https://img.shields.io/badge/license-GPL%20(%3E%3D%202)-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
[![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue.svg)](https://www.python.org/)
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
pip install -e .
```

Runtime dependencies (installed automatically): `pandas`, `numpy`, `scipy`.
Python 3.8 or newer is required.

## Quick start

### Building a network from nodes and links

```python
from netcoin import netCoin

nodes = {
    "name":     ["A", "B", "C", "D"],
    "category": ["x", "x", "y", "y"],
    "size":     [10, 20, 15, 25],
}

links = {
    "Source": [0, 1, 0, 2],
    "Target": [1, 2, 3, 3],
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

### Coincidence analysis from a binary DataFrame

```python
import pandas as pd
from netcoin import coin

df = pd.DataFrame({
    "feature_1": [1, 0, 1, 1],
    "feature_2": [1, 1, 0, 1],
    "feature_3": [0, 1, 1, 0],
})

cooc = coin(df)   # co-occurrence matrix (lower-triangular by default)
print(cooc)
```

### Similarity-based edge lists

```python
from netcoin import sim, edge_list

similarities = sim(df, procedures="j")     # Jaccard by default
edges        = edge_list(df, procedures="j")
```

## Main features

| Area                              | Key API                                    |
|-----------------------------------|--------------------------------------------|
| Interactive coincidence networks  | `netCoin(...)`, `netCoin.plot(...)`        |
| Build a network from a matrix     | `netCoin.fromMatrix(...)`                  |
| Multi-graph visualisations        | `multigraphCreate(...)`                    |
| Co-occurrence matrices            | `coin(...)`                                |
| Similarity measures (Jaccard, …)  | `sim(...)`                                 |
| Edge lists from binary data       | `edge_list(...)`                           |
| Dichotomising categorical data    | `dichotomize(...)`                         |

Visualisations are powered by D3.js. All required web assets (`d3.min.js`,
templates, stylesheets and language files for English, Spanish and Catalan)
ship with the package under `netcoin/www/` and are copied to the output
directory automatically.

## Documentation and examples

Runnable examples live in [`python_examples/`](python_examples/) and cover
several real datasets (Galápagos finches, Italian Renaissance families,
artworks, etc.). Run, for instance:

```bash
python python_examples/python_netCoin_example.py
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

## Authors

- **Modesto Escobar** (Universidad de Salamanca) — creator & maintainer
  [[ORCID]](https://orcid.org/0000-0003-2072-6071)
- David Barrios (Universidad de Salamanca)
- Carlos Prieto (Universidad de Salamanca)
  [[ORCID]](https://orcid.org/0000-0003-2064-4842)
- Luis Martínez-Uribe (Universidad de Salamanca)
  [[ORCID]](https://orcid.org/0000-0002-7795-3972)
- Pablo Cabrera-Álvarez (Universidad de Salamanca)
  [[ORCID]](https://orcid.org/0000-0001-8105-5908)
- Cristina Calvo-López (Universidad de Salamanca)
  [[ORCID]](https://orcid.org/0000-0001-5039-1263)

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
