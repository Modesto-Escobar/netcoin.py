# User guide

This guide follows an analysis from the table of a survey to the page you open in the browser. Every
code block runs, one after the other, in the same session (`tests/test_docs.py` runs them).

Everything that draws returns an object (`netCoin`, `BarCoin`...) that holds its tables
(`net.nodes`, `net.links`) and its options (`net.options`). `net.plot()` writes the page in a temporary
folder and opens it in the browser. To write it where you want, without opening it, give `directory="folder"` when
you build the object or ask for a copy of it with `net.add(directory="folder")`; the folder holds an
`index.html` that works by itself.

## The data

A survey is a table with one row per case and one column per question. Text columns are
the categorical questions, numeric columns the quantitative ones. The examples use a made-up one:

```python
import numpy as np
import pandas as pd

rng = np.random.default_rng(2024)
n = 500
age = rng.choice(["18-34", "35-54", "55+"], n, p=[.3, .4, .3])
works = np.where(rng.random(n) < np.select([age == "18-34", age == "35-54"], [.7, .85], .3), "Yes", "No")
gender = rng.choice(["Female", "Male"], n)
survey = pd.DataFrame({
    "gender": gender,
    "age": age,
    "works": works,
    "union": np.where(rng.random(n) < np.where(works == "Yes", .35, .08), "Member", "Not member"),
    "income": np.round(rng.normal(1000, 250, n) + 700 * (works == "Yes") + 120 * (gender == "Male"), 1),
    "hours": np.round(rng.normal(25, 8, n) + 15 * (works == "Yes"), 1),
    "cweight": np.round(rng.uniform(.5, 2, n), 3),
})
survey.head()
```

## Coincidences

The package started as the study of *coincidences*: how many cases have two characteristics at once.
`dichotomize` turns the categories of a column into 0/1 columns (one per category, or per term when
several are in one cell), `coin` counts what coincides, `sim` and `edge_list` measure how
much, and `all_net` does everything and draws it.

```python
from netcoin import all_net, coin, dichotomize, edge_list, sim

incidences = dichotomize(survey, ["gender", "age", "works", "union"], add=False)
c = coin(incidences)                                  # the Coin object: c.matrix, c.n
jaccard = sim(c, "Jaccard")                           # a square table
edges = edge_list(c, ["Jaccard", "Haberman"], criteria="Haberman", min=0)

net = all_net(incidences, procedures=["Jaccard", "Haberman"], lwidth="Jaccard",
              criteria="Jaccard", min_link=0.1, layout="mc", main="Coincidences",
              directory="out/coincidences")
```

`"Jaccard"`, `"Haberman"` (standardized residual), `"Z"`, `"Hyp"` (hypergeometric), `"Kulczynski"`,
`"Yule"`... are the procedures; `edge_list` keeps the pairs that pass `criteria`. `layout="mc"`
puts the nodes by a correspondence analysis of the incidences, `"pc"` by a principal components
analysis of their correlations.

## A network of a survey

`sur_coin` does the same with the questions of a table: every category is a node (`"question:answer"`),
a quantitative variable in `metric` is a node correlated with the categories, and `weight` weighs
the cases. `cob_coin` is the CoWeb version: only the links among different questions that are
significant, and the nodes on a circle.

```python
from netcoin import cob_coin, sur_coin

net = sur_coin(survey, ["gender", "age", "works", "union"], metric="income", weight="cweight",
               lwidth="Haberman", lcolor="Haberman", link_bipolar=True, color="variable",
               criteria="Haberman", min_link=2, main="Survey")
web = cob_coin(survey, ["gender", "age", "works", "union"], degree=90, significance=0.05,
               directory="out/web")
```

## Correlations

`net_corr` makes every numeric column a node and every pair a link; `method` can be `"pearson"`,
`"spearman"`, `"kendall"`, `"partial"` or `"dcor"`, and several at once. `d_net_corr` shows how the network
splits as the correlation required rises (by default it keeps only the steps where the number of
components grows; `components=False` keeps them all).

```python
from netcoin import d_net_corr, net_corr

quantitative = survey[["income", "hours"]].assign(
    age_years=np.round(30 + survey["income"] / 60 + rng.normal(0, 6, n)),
    score=survey["hours"] / 10 + rng.normal(0, 1.2, n))
corr = net_corr(quantitative, method=["pearson", "spearman"], weight=survey["cweight"], criteria="value",
                min_link=0.0, layout="pc", main="Correlations")
steps = d_net_corr(quantitative, sequence=[0.1, 0.3, 0.5], components=False,
                   directory="out/corr_evolving")
```

## Log-linear models

`log_coin` fits a log-linear model to the table of several categorical variables and draws the
categories (nodes sized by their main effect) joined by the interactions that are significant.
The note of the graph reports the fit of the model.

```python
from netcoin import log_coin

loglin = log_coin(survey, ["gender", "age", "works", "union"], pmax=0.05, frequency=True,
                  percentage=True, order=3)
```

## Regression

`glm_coin` fits one or several equations, one per line (`"response ~ terms, family"`), and draws
the categories and the numeric variables with the effects between them. Its links give
the marginal effect of every category against the weighted mean of its variable (and against its
first category); numeric variables carry their average slope. Standard errors come from the delta
method (the robust HC1 covariance with `robust=True`). It needs the `models` extra.

```python
from netcoin import glm_coin

reg = glm_coin("""income ~ gender + age + works
union ~ income + works, binomial""", survey, weights="cweight", pmax=0.5)
reg.equations       # fit of every equation (R2, AIC, deviance...)
reg.coeffs          # the coefficient tables
```

## Looking for groups

`looking4clusters` runs k-means, PAM and hierarchical clustering for 2 to 10 groups (the best by the
Calinski-Harabasz index) and projects the cases on two dimensions with PCA, classical MDS, t-SNE, NMF
and UMAP. The methods whose library is missing (the `cluster` extra) are skipped with a warning.

```python
from netcoin import looking4clusters

found = looking4clusters(survey[["income", "hours"]], seed=1)
found.optim              # the number of groups each method found best
found.reductions["pca"]  # two coordinates per case
```

The groups and the planes can go straight into a scattergram of the cases; see the
[scattergram guide](scattergram.md).

## Several analyses in one page

`caring` reads a data file and writes one page with the networks that a dictionary of arguments
asks for. The keys are those of the R function and of the graphical interface of netCoin (or their
`snake_case` forms); `plot` alternates the kind of network and its title.

```python
from netcoin import caring

survey.to_csv("out/survey.csv", index=False)
folder = caring("out/survey.csv", {
    "variables": ["gender", "age", "works", "union"],
    "weight": "cweight",
    "plot": ["surCoin", "Associations", "logCoin2", "Log-linear", "surScat", "Cases"],
    "dir": "out/page",
})
```

## Saving and going to other programs

```python
from netcoin import save_pajek, to_igraph

save_pajek(net, "out/survey.net")
graph = to_igraph(net)           # needs python-igraph (extra "igraph")
```

Several graphs can share a page with `multigraphCreate` (a selector, side by side, or a network
that evolves frame by frame).
