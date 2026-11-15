# The scattergram of the cases

`sur_scat` draws the cases of a survey (or the patterns of answers, when many cases answer alike) on
a plane and colours them by groups. The nodes keep the variables, so the page can describe and filter
the groups. Around it there are functions to add planes and groups, to align the planes, to find
communities of cases and to compare groupings.

```python
import numpy as np
import pandas as pd

rng = np.random.default_rng(7)
n = 300
group = rng.choice([0, 1, 2], n)
quant = pd.DataFrame(rng.normal(size=(n, 6)) + group[:, None] * np.array([2, -2, 1, 0, 1, -1]),
                     columns=[f"x{i}" for i in range(1, 7)])
quant["w"] = rng.uniform(.5, 2, n)
quali = pd.DataFrame({
    "q1": np.where(group == 0, rng.choice(["yes", "no"], n, p=[.85, .15]), rng.choice(["yes", "no"], n)),
    "q2": np.where(group == 1, rng.choice(["hi", "lo"], n, p=[.85, .15]), rng.choice(["hi", "lo"], n)),
    "q3": rng.choice(["x", "y", "z"], n),
    "q4": np.where(group == 2, rng.choice(["a", "b"], n, p=[.8, .2]), rng.choice(["a", "b"], n)),
})
```

## The plane and its groups

With qualitative variables the plane is that of a multiple correspondence analysis (`type="mca"`),
with quantitative ones that of a principal components analysis (`type="pca"`); by default the type is
`"pca"` if every active variable is a number. The groups are k-means on the first `nfactors` factors
(`cluster_on="variables"` clusters on the variables instead), numbered from left to right.

```python
from netcoin import sur_scat

scat = sur_scat(quali, nclusters=3, main="Cases")
scat.nodes.head()          # the variables, the groups and the coordinates (fx, fy)
scat.clusters              # the centroid of every group on every plane
pca = sur_scat(quant, ["x1", "x2", "x3", "x4", "x5", "x6"], weight="w", nclusters=[2, 3])
```

`nclusters` is a number, a list of them (one column of groups each; the last colours the graph),
`0` for no groups, or `"auto"`: the number is chosen by `cluster_index` among those stated next to it (2 to
10 if none): `"ch"` (Calinski-Harabasz), `"silhouette"` or `"bic"`. The BIC is that of the latent model that
suits the analysis (profiles for `pca`, classes for `mca`); stating a 1 among the numbers also weighs having
no groups at all.

```python
auto = sur_scat(quant, ["x1", "x2", "x3", "x4", "x5", "x6"], nclusters=["auto", 1, 2, 3, 4], cluster_index="bic")
auto.scat["cluster_index"]    # the BIC of every number of groups; the optimal one is marked
```

Patterns, samples and axes:

```python
many = sur_scat(quali, patterns=True, nclusters=2)       # one node per pattern of answers
sample = sur_scat(quant, ["x1", "x2", "x3"], max_n=100, seed=1)   # at most 100 nodes, drawn at random
observed = sur_scat(quant, ["x1", "x2", "x3"], xaxis="x1", yaxis="x2", nclusters=0)
```

## Rotations, other planes, other groups

`layouts` adds planes next to the main one, and `clusters` adds groups next to the k-means ones (one
column for every number in `nclusters`):

| `layouts` | plane |
|---|---|
| `"varimax"`, `"quartimax"`, `"equamax"` | rotation of the first `nfactors` components (only `type="pca"`); with three or more, one plane per pair |
| `"pca"`, `"cmds"`, `"tsne"`, `"nmf"`, `"umap"` | the two-dimensional views of `looking4clusters` (on the variables as they come) |
| `True` or `"all"` | all of them |

| `clusters` | groups |
|---|---|
| `"kmeans"`, `"pam"`, `"hclust"` | those of `looking4clusters` |
| `"lpa"` / `"lca"` | latent profiles (quantitative variables) / latent classes (qualitative ones) |
| `True` | all of them, with the latent model that suits the analysis |

```python
rich = sur_scat(quant, ["x1", "x2", "x3", "x4", "x5", "x6"], nclusters=[2, 3], nfactors=3,
                layouts=["varimax", "pca"], clusters=["pam", "hclust", "lpa"], suffix="*")
list(rich.layouts)         # PCA, varimax 1-2, varimax 1-2*, ... cPC, cPC*
```

With `suffix`, every plane that was added is also shown aligned (Procrustes) to the main plane, under its
name followed by the suffix; `align_to` chooses another plane (its name, or its position).
`layouts` and `clusters` ask for the same things as the R arguments, but a scattergram has no links, so
the layouts of networks (`"fr"`, `"ka"`...) mean nothing here and are ignored with a warning.

## Adding things afterwards

All these functions take a scattergram and return a modified copy.

```python
from netcoin import add_axes, add_clusters, layouts_align, replace_clusters

mine = rng.integers(1, 4, len(rich.nodes))                       # one group per node
step = add_clusters(rich, mine, name="mine")
step = add_axes(step, rng.normal(size=(len(rich.nodes), 2)), name="noise")
step = layouts_align(step, to=1, suffix="~")                     # every plane laid over the first
step = replace_clusters(step, rng.integers(1, 3, len(rich.nodes)), name="only")   # clear the others
list(step.nodes.columns)[-4:]
```

When the nodes are patterns, a vector of one value per *case* is accepted and collapsed into the
patterns (a pattern that spans several groups shows them all, separated by `|`).

`looking4clusters` results go in as they are: `add_clusters(scat, found)` adds every grouping it
found and `add_axes(scat, found)` every reduction.

## Communities of cases

`add_communities` correlates the cases with one another over some variables, takes the communities of
that network (needs the `igraph` extra) and adds them as groups. The links are weighed by their correlation;
communities with fewer than `min_size` cases (a proportion if below 1) are gathered into a residual group.

```python
from netcoin import add_communities

linked = add_communities(rich, quant, ["x1", "x2", "x3", "x4", "x5", "x6"], community=["walktrap", "louvain"],
                         min_link=0.5, min_size=0.05)
[c for c in linked.nodes.columns if c.startswith(("walktrap", "louvain"))]
```

`community` takes codes or names: `"ed"` edgeBetweenness, `"fa"` fastGreedy, `"la"` labelProp, `"le"`
leadingEigen, `"lo"` louvain, `"op"` optimal, `"sp"` spinglass, `"wa"` walktrap. A network of cases is
large, so `edgeBetweenness` and `spinglass` refuse more than 2000 nodes and `optimal` more than 100
(`max_nodes` overrides). `weighted=False` finds the communities of the unweighted network, which is what
R does (see [Differences from R](r-to-python.md)).

## Comparing groupings

`cluster_match` compares a model classification with the others, category by category: the categories
are matched one to one so that the diagonal is as large as possible, and every classification gets its
accuracy, kappa, sensitivity, precision, F1 and specificity. `cluster_agree` measures the agreement
between partitions (adjusted Rand index, normalized mutual information, variation of information),
whatever the number of categories.

```python
from netcoin import cluster_agree, cluster_match

linked = add_clusters(linked, group, name="truth")          # the groups that made the data
cluster_match(linked, "truth", columns=["Groups(3)", "pam(3)", "hclust(3)", "LPA(3)"])
cluster_agree(linked, "truth")
```

The result is a table with one row per classification; its `attrs` hold the confusion tables
(`tables`), the matching (`mapping`) and, for ties, `ties`.
