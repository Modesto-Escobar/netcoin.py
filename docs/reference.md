# Reference

Generated from the docstrings by `python tools/make_reference.py`; do not edit it by hand.

## The network and its page

### `netCoin`

```python
class netCoin(nodes=None, links=None, tree=None, community: 'str | None' = None, layout=None, name: 'str | None' = None, label: 'str | None' = None, group: 'str | None' = None, group_text: 'bool' = False, label_size: 'str | None' = None, size=None, color=None, shape=None, border: 'str | None' = None, legend: 'str | None' = None, sort: 'str | None' = None, decreasing: 'bool' = False, ntext: 'str | None' = None, info: 'str | None' = None, image: 'Sequence[str] | str | None' = None, image_names: 'Sequence[str] | str | None' = None, centrality=None, node_bipolar: 'bool' = False, node_scale_limits: 'Sequence[float] | None' = None, node_filter=None, degree_filter=None, lwidth: 'str | None' = None, lweight: 'str | None' = None, lcolor=None, ltext: 'str | None' = None, intensity: 'str | None' = None, link_bipolar: 'bool' = False, link_scale_limits: 'Sequence[float] | None' = None, link_filter=None, repulsion: 'float' = 25, distance: 'float' = 10, zoom: 'float' = 1, fixed: 'bool | None' = None, limits: 'Sequence[float] | None' = None, main: 'str | None' = None, note: 'str | None' = None, show_coordinates: 'bool' = False, show_arrows: 'bool' = False, show_legend: 'bool' = True, frequencies: 'bool' = False, statistics: 'bool' = False, show_axes: 'bool' = False, axes_labels: 'Sequence[str] | None' = None, nodes_label: 'str | None' = None, links_label: 'str | None' = None, scenarios: 'float | None' = None, help: 'str | None' = None, help_on: 'bool' = False, mode: 'Sequence[str] | str' = ('network', 'heatmap'), rounded_items: 'bool' = False, controls: 'Sequence[int] | int | None' = (1, 2, 3, 4, 5, 6, 7, 8), cex: 'float' = 1, background: 'str | None' = None, default_color: 'str' = '#1f77b4', language: 'str' = 'en', directory: 'str | None' = None, *, source: 'str' = 'Source', target: 'str' = 'Target', node_attr_names: 'Sequence[str] | None' = None, link_attr_names: 'Sequence[str] | None' = None, **legacy)
```

An interactive network: nodes, links and the options of how to draw them.

Port of `netCoin` (R). The tables are DataFrames: `nodes` and `links` (`None`
when there are none), `tree` (`Source`, `Target`) and `options` is the
dict that the web page reads. Draw it with `plot` or `directory=`.

**Parameters**

- `nodes` (DataFrame or dict, optional): One row per node.
- `links` (DataFrame or dict, optional): One row per link, with `Source` and `Target` (names of the nodes; if they are numbers, 1-based positions in `nodes`). If there are no nodes they are taken from the links.
- `tree` (DataFrame or list of pairs, optional): Parent -> child relations among the nodes.
- `community` (str, optional): Detect communities with igraph (`"walktrap"`, `"louvain"`, `"fast_greedy"`, `"edge_betweenness"`, `"label_prop"`, `"leading_eigen"`, `"optimal"` or `"spinglass"`; two letters are enough) and put them in the node column `community`.
- `layout` (str, matrix, DataFrame, list or dict, optional): Coordinates: a name (`"fruchterman_reingold"`, `"kamada_kawai"`, `"circle"`, `"grid"`, `"star"`, `"reingold_tilford"`, `"sugiyama"`, `"mds"`, `"drl"`, `"lgl"`, `"davidson_harel"`, `"random"`, `"force_atlas2"`; two letters are enough), a matrix with two columns, or several of them in a dict (planes the page lets you switch between).
- `name, label, group, group_text, label_size, size, color, shape, border, legend,`
- `sort, decreasing, ntext, info`: Node columns used for the name, label, groups, label size, size (or a number for all), colour, shape (`"Circle"`, `"Square"`, `"Diamond"`, `"Triangle"`, `"Cross"`, `"Star"`, `"Wye"`), border, legend, order, text of the tooltip and the info panel. `color` and `shape` can also be one value, one per node, or a dict `{column: {value: colour}}`.
- `image, image_names` (list of str, optional): Node columns with image files (and the names to show).
- `centrality` (list of str, optional): Add the centralities `"wdegree"`, `"closeness"`, `"betweenness"`, `"eigen"` as node columns (needs igraph).
- `node_bipolar, link_bipolar` (bool): Treat negative and positive values as two poles of the colour scale.
- `node_scale_limits, link_scale_limits` (list of 2 or 3 numbers): Limits of the numeric colour scale of `color` / `lcolor`.
- `node_filter, link_filter` (str, callable, or sequence of bool): Which nodes / links are shown. A text is a pandas expression over the columns (`"size > 10 and group == 'a'"`); a callable gets the table and returns the booleans.
- `degree_filter` (number or (min, max)): Show only the nodes whose degree is in the range.
- `lwidth, lweight, lcolor, ltext, intensity`: Link columns for the width, the weight (for layouts), the colour (also one value or one per link), the text and the intensity of the colour.
- `repulsion, distance, zoom` (number): Behaviour of the force layout (0-100, 0-100, 0.1-10).
- `fixed, limits, show_coordinates, show_arrows, show_axes, axes_labels`: Fix the positions, give the area (x1, y1, x2, y2), show the coordinates in the tables, draw arrows and axes (with their labels).
- `main, note, help, help_on, scenarios, cex, background, default_color, language`: Titles, help text (`help_on` shows it at the start), the number of scenarios, font size, background colour or image file, colour of the nodes, and `"en"`/`"es"`/`"ca"`.
- `show_legend, frequencies, statistics, nodes_label, links_label`: Show the legend, the frequency column, the statistics panel and the labels of the two tables.
- `mode` (list of str): `"network"`, `"heatmap"` or both.
- `controls` (list of int): Which controls are shown: 1 sidebar, 2 selection buttons, 3 tables, 4 sliders and buttons, 5 export buttons, 6, 7, 8 legend (0 hides all, negatives remove).
- `rounded_items` (bool): Rounded icons in the tables.
- `directory` (str, optional): Write the page there.

**Notes**

The camelCase names of earlier versions (`nodeFilter`, `showArrows`, `defaultColor`...) still work but are deprecated; `orderA`/`orderD` became `sort` and `decreasing`.

### `multigraphCreate`

```python
def multigraphCreate(graphs, names=None, mode='default', directory='MultiGraph', show=True, *, mfrow=(1, 2), frame=0, speed=50, loop=False, lineplots=None) -> 'str'
```

Write several graphs into one folder and (optionally) open it.

Port of `multigraphCreate` from the R package netCoin.

**Parameters**

- `graphs` (list or dict): The graphs: `netCoin` objects or folders that already hold a graph (with an `index.html`). A dict gives the names.
- `names` (list of str, optional): Name of each graph (`graph1`, `graph2`... by default).
- `mode` ({"default", "parallel", "frame"}): `default`: a selector to move from one graph to another; `parallel`: all of them at the same time, in a grid of `mfrow`; `frame`: the networks are merged into one that evolves, one frame per graph.
- `directory` (str): Folder to write. An existing folder is only replaced if it is a graph made by this package.
- `show` (bool): Open it in the browser.
- `mfrow` ((rows, columns)): Grid of the parallel mode.
- `frame, speed, loop, lineplots`: Frame mode: the first frame shown, the speed (0-100), whether it repeats, and the variables to draw as line plots (a node or link column, or `"degree"`).

**Returns**

- `str`: The folder written.

### `bar_coin`

```python
def bar_coin(data: 'pd.DataFrame', variables: 'Sequence[str] | None' = None, commonlabel: 'str | Sequence[str] | None' = None, dichotomies: 'str | Sequence[str]' = '_all', value_dicho=1, weight: 'str | Sequence[float] | None' = None, subsample: 'bool' = False, sort: 'str | None' = None, decreasing: 'bool' = True, nodes: 'pd.DataFrame | None' = None, name: 'str | None' = None, select: 'str | None' = None, scalebar: 'bool' = False, note: 'str | None' = None, label: 'str | None' = None, text: 'str | None' = None, color: 'str | None' = None, default_color: 'str' = '#1f77b4', expected: 'bool' = False, confidence: 'bool' = False, level: 'float' = 0.95, significance: 'bool' = False, minimum: 'float' = 1, maximum: 'float | None' = None, percentages: 'bool' = False, criteria: 'str' = 'Z', bonferroni: 'bool' = False, support: 'float' = 1, min_link: 'float' = -inf, max_link: 'float' = 1, language: 'str' = 'en', cex: 'float' = 1.0, directory: 'str | None' = None) -> 'BarCoin'
```

Bar chart of the coincidences of every variable (category) with the others.

Port of `barCoin` from the R package netCoin.

**Parameters**

- `data` (DataFrame): The survey or table of cases.
- `variables` (list of str, optional): Columns to use (all by default).
- `commonlabel` (str or list, optional): Variables whose name is kept in the label of their categories (`"variable:category"`); for the rest only the category is shown.
- `dichotomies` (`"_all"`, `"_none"` or list of str): Which variables are dichotomies (`value_dicho` against the rest, e.g. 0/1 columns) and which are categorical (one bar per category).
- `value_dicho`: The value that counts as 1 in the dichotomies.
- `weight` (str or sequence, optional): Weight of each case or the column that holds it.
- `subsample` (bool): Use only the cases with at least one event.
- `sort, decreasing`: Column of the nodes to order the events by.
- `nodes` (DataFrame, optional): Your own table of nodes; it restricts and orders the events and adds its columns (`label`, `text`, `color` can refer to them).
- `name, select, scalebar, note, label, text, color, default_color, cex, language`: As in the page options.
- `expected, confidence, significance` (bool): Draw the expected coincidences, their confidence limits (at `level`) and the p-value of the link (`criteria`: `"Z"` or `"hyp"`).
- `minimum, maximum` (float): Frequency range of the events.
- `percentages` (bool): Show percentages of the cases instead of counts.
- `bonferroni, support, min_link, max_link`: Filters of the links, as in `edge_list` (`support` is the least number of coincidences and `max_link` bounds the p-value).
- `directory` (str, optional): Write the page there.

**Returns**

- `BarCoin`

### `BarCoin`

```python
class BarCoin(events: 'pd.DataFrame', links: 'pd.DataFrame', name: 'str | None' = None, select: 'str | None' = None, source: 'str | None' = None, target: 'str | None' = None, label: 'str | None' = None, text: 'str | None' = None, color: 'str | None' = None, incidences: 'str | None' = None, coincidences: 'str | None' = None, expected: 'str | None' = None, confidence: 'str | Sequence[str] | None' = None, level: 'float' = 0.95, significance: 'str | None' = None, sort: 'str | None' = None, decreasing: 'bool' = False, scalebar: 'bool' = False, default_color: 'str' = '#1f77b4', note: 'str | None' = None, cex: 'float' = 1, language: 'str' = 'en', directory: 'str | None' = None)
```

A bar chart: for each event, the bars of its coincidences with the others.

Port of rD3plot's `barplot_rd3`. Build it with `bar_coin`.

**Attributes**

- `nodes` (DataFrame): One row per event, with its name and number of incidences.
- `links` (DataFrame): `Source`, `Target` and the statistics of every pair (`coincidences`, `expected`, confidence limits, p-value...).
- `options` (dict): What the page reads (column names, colours, language...).

### `time_coin`

```python
def time_coin(periods: 'pd.DataFrame', name: 'str' = 'name', start: 'str' = 'start', end: 'str' = 'end', group: 'str | None' = None, text: 'str | None' = None, main: 'str | None' = None, note: 'str | None' = None, info: 'str | None' = None, events: 'pd.DataFrame | None' = None, event_names: 'str' = 'name', event_period: 'str' = 'period', event_time: 'str' = 'date', event_color: 'str | None' = None, event_shape: 'str | None' = None, cex: 'float' = 1, language: 'str' = 'en', directory: 'str | None' = None) -> 'TimeCoin'
```

Timeline of periods (and events) from a table.

Port of `timeCoin` from the R package netCoin.

**Parameters**

- `periods` (DataFrame): One row per period: `name`, `start` and `end` (numbers or date-times).
- `group, text, info`: Columns used to colour the bars, to describe them and for the info panel.
- `events` (DataFrame, optional): Events: the `event_names` column names them, `event_period` says which period they belong to (a name of `periods`) and `event_time` when (a number or a date-time).
- `event_color, event_shape`: Columns of `events` that set the colour and the shape of the marks.
- `main, note, cex, language, directory`: Title, footnote, font size, `"en"`/`"es"`/`"ca"` and where to write the page.

**Returns**

- `TimeCoin`

### `TimeCoin`

```python
class TimeCoin(periods: 'pd.DataFrame', name: 'str' = 'name', start: 'str' = 'start', end: 'str' = 'end', group: 'str | None' = None, text: 'str | None' = None, main: 'str | None' = None, note: 'str | None' = None, info: 'str | None' = None, events: 'pd.DataFrame | None' = None, event_names: 'str' = 'name', event_period: 'str' = 'period', event_time: 'str' = 'date', event_color: 'str | None' = None, event_shape: 'str | None' = None, cex: 'float' = 1, language: 'str' = 'en', directory: 'str | None' = None)
```

A timeline: bars for periods (people, schools...), optionally with events on them.

Port of rD3plot's `timeline_rd3`. Build it with `time_coin`.

**Attributes**

- `periods` (DataFrame): One row per period, with its name, start and end.
- `events` (DataFrame or None): Events placed on the periods.
- `options` (dict): What the page reads.

### `pie_coin`

```python
def pie_coin(x: 'Coin', colors: 'Sequence[str]' = ('#000000', '#8dc7e6', '#ffffff', '#005587'), nodes: 'pd.DataFrame | None' = None, links: 'pd.DataFrame | None' = None, name: 'str | None' = None, lcolor=None, expected: 'bool' = True, abline: 'int | None' = None, main: 'str | None' = None, note: 'str | None' = None, show_legend: 'bool' = True, help: 'str | None' = None, help_on: 'bool' = False, cex: 'float' = 1, language: 'str' = 'en', directory: 'str | None' = None) -> 'PieCoin'
```

A matrix of pie charts: every pair of variables, split by who has what.

Port of `pieCoin` from the R package netCoin. Each pie shows the cases that
have both characteristics, only the row one, only the column one and
neither; with `expected` the expected coincidences are marked.

**Parameters**

- `x` (Coin): See `coin`.
- `colors` (four colours): For both, row only, neither and column only (in the order of the pie).
- `nodes, links` (DataFrame, optional): Information about the variables and the pairs (links from `Source` rows to `Target` columns).
- `lcolor`: Colour of the links: a column of `links`, a colour, or one value per link.
- `expected` (bool): Mark the expected coincidences.
- `abline` (int, optional): Draw a line after that many rows and columns.
- `main, note, show_legend, help, help_on, cex, language, directory`: Page options.

**Returns**

- `PieCoin`

### `PieCoin`

```python
class PieCoin(v: 'np.ndarray', w: 'np.ndarray | None' = None, labels: 'Sequence[str] | None' = None, colors: 'Sequence[str] | None' = None, nodes: 'pd.DataFrame | None' = None, links: 'pd.DataFrame | None' = None, name: 'str | None' = None, source: 'str | None' = None, target: 'str | None' = None, lcolor=None, abline_x: 'int | None' = None, abline_y: 'int | None' = None, hide_upper: 'bool' = False, main: 'str | None' = None, note: 'str | None' = None, show_legend: 'bool' = True, help: 'str | None' = None, help_on: 'bool' = False, cex: 'float' = 1, language: 'str' = 'en', rownames: 'Sequence[str] | None' = None, colnames: 'Sequence[str] | None' = None, directory: 'str | None' = None)
```

A matrix of pie charts, one per pair of variables.

Port of rD3plot's `pie_rd3`. Build it with `pie_coin`.

**Attributes**

- `data` (ndarray): `[row, column, part]` counts of each pie.
- `dataw` (ndarray or None): The expected counts, the same way.
- `rownames, colnames` (list of str)
- `nodes, links` (DataFrame or None): Optional information of the variables and of the pairs.
- `options` (dict)

### `inc_time`

```python
def inc_time(data: 'pd.DataFrame', name: 'str' = 'name', beginning: 'str' = 'birth', end: 'str' = 'death') -> 'pd.DataFrame'
```

Which of the people (or periods) are alive in every year.

Port of `incTime`: a table with one row per year, from the first beginning to the
last end, and one 0/1 column per person. An open end (missing) goes on forever.
A missing beginning gives missing values.

**Parameters**

- `data` (DataFrame): One row per person or period.
- `name, beginning, end` (str): Columns with the name, the first year and the last year.

**Returns**

- `DataFrame`: One row per year, one 0/1 column per person.

### `mobile_edges`

```python
def mobile_edges(data: 'pd.DataFrame', name: 'int | str' = 0, number: 'int | str' = 1, difference: 'float' = 0) -> 'pd.DataFrame'
```

Links between the items whose numbers are close (`|a - b| <= difference`).

Port of `mobileEdges`: `sim.` is `1 + difference - |a - b|` (at least 1) and
`dist.` its complement `1 + difference - sim.`. Use it, for instance, to link
people born in the same period.

**Parameters**

- `data` (DataFrame): One row per item.
- `name, number` (int or str): Columns (position or name) with the item and its number.
- `difference` (float): Largest difference of numbers that still links two items.

**Returns**

- `DataFrame`: An edge list (`Source`, `Target`, `sim.`, `dist.`).

## Coincidences

### `dichotomize`

```python
def dichotomize(data: 'pd.DataFrame', variables: 'str | Sequence[str]', sep: 'str | Sequence[str] | dict' = '', regex: 'bool' = False, text: 'str' = 'none', min_freq: 'float' = 1, length: 'int' = 0, sort: 'bool' = True, values: 'Sequence | None' = None, stopwords: 'Iterable[str] | None' = None, nas: 'str | None' = 'None', sparse: 'bool' = False, add: 'bool' = False) -> 'pd.DataFrame'
```

Convert categorical or textual columns into 0/1 indicator columns.

Port of `dichotomize` from the R package netCoin. The columns are named
after the categories (or the words) found in the selected variables.

**Parameters**

- `data` (DataFrame): Table with the variables to convert. A cell can hold a single value or several events delimited by `sep`.
- `variables` (str or list of str): Columns to convert.
- `sep` (str, sequence of str or dict): Separator that divides the events of a cell. `""` (default) takes the whole cell as one category. `"C"` concatenates the variables of each row before dichotomising (requires `nas`). A sequence gives one separator per variable (a shorter one recycles its first element).
- `regex` (bool): Treat `sep` as a regular expression. With several variables, their terms are pooled per row into a single group named after the first variable.
- `text` ({"none", "lower", "upper"}): Case conversion. `"lower"`/`"upper"` also remove punctuation and numbers.
- `min_freq` (float): Keep only terms with at least this frequency. A value in (0, 1) is a proportion of the rows. (`min` in R.)
- `length` (int): Maximum number of columns (most frequent terms first); 0 is no limit.
- `sort` (bool): Order the columns by decreasing frequency. Otherwise they follow the alphabetical order of the terms. Ties are broken alphabetically.
- `values` (sequence, optional): Restrict the conversion to these terms.
- `stopwords` (iterable of str, optional): Words to remove; only applied when `text` is `"lower"` or `"upper"`.
- `nas` (str or None): Name of the column that flags the rows without any term; it is called `"<variable>:<nas>"`. `None` omits it (always omitted if `sparse`).
- `sparse` (bool): Return pandas sparse columns instead of dense ones.
- `add` (bool): Append the new columns to `data` (ignored if `sparse`).

**Returns**

- `DataFrame`: 0/1 integer columns, with the index of `data`. With several independent variables the columns are prefixed with `"<variable>:"`.

**Notes**

The order of terms with the same frequency follows an approximation of the locale collation used by R (case and accents are secondary). Empty tokens (`"a;;b"`) are discarded, whereas R keeps an unnamed column.

### `coin`

```python
def coin(incidences: 'pd.DataFrame | np.ndarray', minimum: 'float' = 1, maximum: 'float | None' = None, sort: 'bool' = False, decreasing: 'bool' = True, total: 'bool' = False, subsample: 'bool' = False, weight: 'Sequence[float] | None' = None, pairwise: 'bool' = False) -> 'Coin'
```

Count the coincidences among the variables of an incidence table.

Port of `coin` from the R package netCoin.

**Parameters**

- `incidences` (DataFrame or ndarray): One row per case, one 0/1 column per variable (see `dichotomize`). Missing values are allowed.
- `minimum, maximum` (float): Keep the variables whose frequency lies in `[minimum, maximum]`. `maximum` defaults to the number of cases (the total weight with `weight`).
- `sort, decreasing` (bool): Order the variables by frequency (stable for ties).
- `total` (bool): Add a first variable `"Total"` that is present in every case. Not applicable with weights.
- `subsample` (bool): Use only the cases with at least one event.
- `weight` (sequence of float, optional): Weight of each case.
- `pairwise` (bool): Treat missing values pairwise (each pair uses the cases where both variables are known) instead of dropping the incomplete cases.

**Returns**

- `Coin`

### `Coin`

```python
class Coin(matrix: 'pd.DataFrame', n: 'float', m: 'pd.DataFrame | None' = None, x: 'pd.DataFrame | None' = None) -> None
```

Coincidences among a set of events or characteristics.

This is what `coin` returns; `sim`, `edge_list` and the
rest of the functions take it as input. It plays the role of the `coin`
class of the R package.

**Attributes**

- `matrix` (DataFrame): Square, symmetric table with the number of cases in which each pair of variables coincides. The diagonal holds the frequency of each variable.
- `n` (float): Number of cases (the sum of the weights, if they were used).
- `m` (DataFrame, optional): Pairwise: number of cases with no missing value in each pair.
- `x` (DataFrame, optional): Pairwise: cases where the row variable is missing and the column variable is present.

### `coocur`

```python
def coocur(ocurrences: 'pd.DataFrame | np.ndarray', minimum: 'float' = 1, maximum: 'float' = inf, sort: 'bool' = False, decreasing: 'bool' = True) -> 'Cooc'
```

Co-occurrence of counts: for each pair, the sum over cases of the smaller count.

Port of `coocur` from the R package netCoin. Unlike `coin`, the
variables need not be 0/1.

**Parameters**

- `ocurrences` (DataFrame or array): One row per case and one column per variable, with counts.
- `minimum, maximum` (float): Keep the variables whose total lies in this range.
- `sort, decreasing` (bool): Sort the variables by their totals.

**Returns**

- `Cooc`

### `Cooc`

```python
class Cooc(matrix: 'pd.DataFrame', n: 'float', m: 'float') -> None
```

Result of `coocur`: occurrences shared by pairs of variables.

`n` is the total number of occurrences and `m` the sum, over the cases,
of the largest occurrence of each case.

**Attributes**

- `matrix` (DataFrame): Square matrix of shared occurrences (the diagonal holds the totals).
- `n, m` (float): Total occurrences and sum of the largest occurrence of each case.

### `sim`

```python
def sim(data: 'Coin | pd.DataFrame | np.ndarray', procedures: 'str | Sequence[str]' = 'Jaccard', level: 'float' = 0.95, distance: 'bool' = False, minimum: 'float' = 1, maximum: 'float' = inf, sort: 'bool' = False, decreasing: 'bool' = False, weight: 'Sequence[float] | None' = None, pairwise: 'bool' = False) -> 'pd.DataFrame | dict[str, pd.DataFrame]'
```

Similarities, associations and statistics for every pair of variables.

Port of `sim` from the R package netCoin.

**Parameters**

- `data` (Coin or DataFrame): A `Coin` object or a table of incidences (see `coin`).
- `procedures` (str or list of str): What to compute: `"matching"`, `"Rogers"`, `"Gower"`, `"Sneath"`, `"Anderberg"`, `"Jaccard"`, `"dice"`, `"antiDice"`, `"Ochiai"`, `"Kulczynski"`, `"Hamann"`, `"Yule"`, `"Pearson"`, `"tetrachoric"`, `"odds"`, `"Russell"`, `"expected"`, `"confidence"` (also `conf.L` and `conf.U`), `"Haberman"`, `"Z"`, `"Fisher"`, `"coincidences"`, `"relative"`, `"sConditional"`, `"tConditional"`, `"c.conditional"` and `"c.probable"`; or `"all"`. R's short forms work too (`"J"`, `"Ham"`, `"CC"`, `"Hyp"`...).
- `level` (float): Confidence level for `"confidence"`.
- `distance` (bool): Return distances (`1 - similarity`) for the similarity coefficients.
- `minimum, maximum` (float): Keep the variables whose frequency lies in that range.
- `sort, decreasing` (bool): Order the variables by frequency.
- `weight, pairwise`: As in `coin` (used when `data` is a table of incidences; `pairwise=True` with a `Coin` needs one built pairwise).

**Returns**

- `DataFrame or dict of DataFrame`: One square matrix per procedure: a DataFrame if only one was asked for, otherwise a dict keyed by the name of the matrix.

### `edge_list`

```python
def edge_list(data: 'Coin | pd.DataFrame | np.ndarray', procedures: 'str | Sequence[str]' = 'Haberman', criteria: 'str' = 'Z', level: 'float' = 0.95, bonferroni: 'bool' = False, min: 'float' = -inf, max: 'float' = inf, support: 'float' = -inf, directed: 'bool' = False, diagonal: 'bool' = False, sort: 'str | None' = None, decreasing: 'bool' = True, pairwise: 'bool' = False) -> 'pd.DataFrame'
```

List of links (with their weights) between pairs of variables.

Port of `edgeList` from the R package netCoin.

**Parameters**

- `data` (Coin, DataFrame or ndarray): A `Coin` object or a table of incidences. With `procedures="shape"` it is a numeric matrix whose cells become the weights; with `"shape"` (or `"tree"`) on a table of text, an adjacency list: the first column is linked to the others.
- `procedures` (str or list of str): Statistics to include (see `sim`), or `"shape"`/`"tree"`.
- `criteria` (str): Statistic that `min` and `max` refer to. For p-values (`"Z"` and `"Fisher"`) `max` defaults to 0.5.
- `level` (float): Confidence level for `"confidence"`.
- `bonferroni` (bool): Divide a p-value `max` by the number of pairs.
- `min, max` (float): Keep the pairs whose `criteria` lies in `[min, max]`.
- `support` (float): Keep the pairs with at least this many coincidences.
- `directed, diagonal` (bool): Include both directions / the diagonal.
- `sort` (str, optional): Statistic to order the links by (`decreasing` sets the direction).
- `pairwise` (bool): As in `coin`, for a table of incidences.

**Returns**

- `DataFrame`: Columns `Source`, `Target` and one per statistic. Empty (with a warning) if no pair meets the conditions. As in R, the list goes column by column of the lower triangle, so `Source` is the variable that comes first.

### `expected_list`

```python
def expected_list(data: 'Coin', names: 'Sequence[str] | None' = None, min: 'float' = 1, confidence: 'bool' = False) -> 'pd.DataFrame'
```

Observed and expected coincidences for every pair.

Port of `expectedList` from the R package netCoin. With
`confidence=True` it adds the upper/lower limit of the expectation
(`confidence`; R uses a fixed 1.64, i.e. a one-sided 95 % level).

**Parameters**

- `data` (Coin): The coincidences.
- `names` (list of str, optional): Restrict the pairs to these variables.
- `min` (float): Smallest number of coincidences of a pair to be listed.
- `confidence` (bool): Add the limits of the expectation.

**Returns**

- `DataFrame`: One row per pair (`Source`, `Target`, `coincidences`, `expected`, and the limits with `confidence`).

### `prop_coin`

```python
def prop_coin(x: 'Coin', margin: 'int' = 0, decimals: 'int' = 1) -> 'pd.DataFrame'
```

Coincidences as percentages.

Port of `propCoin` from the R package netCoin.

**Parameters**

- `x` (Coin)
- `margin` ({0, 1, 2}): 0: percentage of the cases (of the cases with data in both variables, if pairwise); 1: of the frequency of the row variable; 2: of the column's.
- `decimals` (int): Rounding.

### `as_nodes`

```python
def as_nodes(data: 'Coin | pd.DataFrame', frequency: 'bool' = True, percentages: 'bool' = False, language: 'str' = 'en') -> 'pd.DataFrame'
```

Table of nodes (name, frequency and/or percentage) from a `Coin`.

Port of `asNodes` from the R package netCoin. If `data` is an edge list
(with `Source` and `Target` columns) the nodes are the sorted names found
in it.

**Parameters**

- `data` (Coin or DataFrame): The coincidences, or an edge list.
- `frequency, percentages` (bool): Add the frequency and/or the percentage of each variable.
- `language` ({"en", "es", "ca"}): Language of the column names.

**Returns**

- `DataFrame`: One row per node, ready for `netCoin`.

### `all_net`

```python
def all_net(incidences: 'pd.DataFrame', weight: 'str | Sequence[float] | None' = None, subsample: 'bool' = False, pairwise: 'bool' = False, minimum: 'float' = 1, maximum: 'float | None' = None, sort: 'bool' = False, decreasing: 'bool' = True, frequency: 'bool' = False, percentages: 'bool' = True, procedures: 'str | Sequence[str]' = 'Haberman', criteria: 'str' = 'Z', bonferroni: 'bool' = False, level: 'float' = 0.95, support: 'float' = -inf, min_link: 'float' = -inf, max_link: 'float' = inf, directed: 'bool' = False, diagonal: 'bool' = False, sort_links: 'str | None' = None, decreasing_links: 'bool' = True, igraph: 'bool' = False, **netcoin_args)
```

Build a `netCoin` network from a table of incidences.

Port of `allNet` from the R package netCoin: it counts the coincidences
(`coin`), turns the variables into nodes (`as_nodes`), the
pairs into links (`edge_list`) and draws them with `netCoin`.

**Parameters**

- `incidences` (DataFrame): One row per case, one 0/1 column per variable (see `dichotomize`).
- `weight` (str or sequence, optional): Weight of each case, or the name of the column that holds it.
- `subsample, pairwise, minimum, maximum, sort, decreasing`: As in `coin`.
- `frequency, percentages` (bool): Put the frequency and/or the percentage of each variable in the nodes. With `percentages` the nodes are sized by `"%"`.
- `procedures, criteria, bonferroni, level, support`: As in `edge_list`; `procedures` also gets the link attributes that `lwidth`, `lweight`, `lcolor` and `ltext` refer to.
- `min_link, max_link` (float): Range of `criteria` kept as links (`minL` and `maxL` in R).
- `directed, diagonal` (bool): As in `edge_list`.
- `sort_links, decreasing_links`: Order of the links (`sortL` and `decreasingL` in R).
- `igraph` (bool): Return an `igraph.Graph` instead (needs the optional `igraph` package).
- `**netcoin_args`: Anything else is passed to `netCoin` (`name`, `nodes`, `language`, `color`, `layout`...). `layout` can also be `"mc"` (correspondence analysis) or `"pc"` (principal components of the correlations); `directory` writes the web page.

## Surveys

### `sur_coin`

```python
def sur_coin(data: 'pd.DataFrame', variables: 'Sequence[str] | None' = None, commonlabel: 'str | Sequence[str] | None' = None, dichotomies: 'str | Sequence[str] | None' = None, value_dicho=1, metric: 'str | Sequence[str] | None' = None, exogenous: 'str | Sequence[str] | None' = None, weight: 'str | Sequence[float] | None' = None, subsample: 'bool' = False, pairwise: 'bool' = False, minimum: 'float' = 1, maximum: 'float | None' = None, sort: 'bool' = False, decreasing: 'bool' = True, frequency: 'bool' = False, percentages: 'bool' = True, procedures: 'str | Sequence[str]' = 'Haberman', criteria: 'str' = 'Z', bonferroni: 'bool' = False, support: 'float' = -inf, min_link: 'float' = -inf, max_link: 'float' = inf, directed: 'bool' = False, diagonal: 'bool' = False, sort_links: 'str | None' = None, decreasing_links: 'bool' = True, igraph: 'bool' = False, coin: 'bool' = False, directory: 'str | None' = None, **netcoin_args)
```

A network of the categories of a survey (port of `surCoin` of the R package).

Every category of the `variables` becomes a node (`"variable:category"`), and the
links tell which categories go together more (or less) than expected: the pairs are
counted with `coin` and described with `edge_list`.

**Parameters**

- `data` (DataFrame): One row per case; the variables can be text, categories or numbers. A missing value is a case that has none of the categories of that variable.
- `variables` (list of str, optional): The qualitative variables (all the columns by default).
- `commonlabel` (str or list, optional): Variables whose name stays in the label of their categories (`"variable:category"`); in the others only the category is shown. Creates the column `label`.
- `dichotomies` (str or list, optional): Variables to take as yes/no: 1 where they equal `value_dicho`. They make one node each, named after the variable.
- `value_dicho`: The value that counts as yes in the `dichotomies` (one for all, or one per variable).
- `metric` (str or list, optional): Quantitative variables. Each makes a node (its mean, rescaled to 0-100, is its `"%"`) and is correlated with the rest, with the Pearson coefficient, the Haberman residual and the p-value of the correlation (`procedures` and `criteria` are limited to those three).
- `exogenous` (str or list, optional): Variables (or `metric`/`dichotomies`) that cannot be explained one by another: the links between two of them are hidden (`chaine == "No"`).
- `weight` (str or sequence, optional): The weight of each case, or the name of the column that holds it. With it, the cases with a missing value are left out (unless `pairwise`).
- `subsample, pairwise, minimum, maximum, sort, decreasing`: As in `coin`; `minimum` also limits the categories. `maximum` is the number of cases by default.
- `frequency, percentages` (bool): Put the frequency and/or the percentage of each category in the nodes. With `percentages` the nodes are sized by `"%"`.
- `procedures, criteria, bonferroni, support`: As in `edge_list`; `procedures` also takes the link attributes that `lwidth`, `lweight`, `lcolor` and `ltext` refer to.
- `min_link, max_link` (float): Range of `criteria` kept as links (`minL` and `maxL` in R). For a p-value `max_link` is 0.5 unless you say otherwise.
- `directed, diagonal` (bool): As in `edge_list`.
- `sort_links, decreasing_links`: Order of the links (`sortL` and `decreasingL` in R).
- `igraph` (bool): Return an `igraph.Graph` instead of a `netCoin`.
- `coin` (bool): Return the `Coin` object of the categories and stop there.
- `directory` (str, optional): Write the web page there.
- `**netcoin_args`: The rest goes to `netCoin` (`nodes`, `name`, `language`, `level`, `color`, `layout`, `lwidth`, `main`...). `nodes` restricts and orders the categories and adds its columns; `layout` can also be `"mc"` (correspondence analysis) or `"pc"` (principal components of the correlations), or a matrix with one row per row of `nodes`.

**Returns**

- `netCoin, Coin or igraph.Graph`

### `cob_coin`

```python
def cob_coin(data: 'pd.DataFrame', variables: 'Sequence[str] | None' = None, degree: 'float' = 0, significance: 'float' = 0.05, **kwargs)
```

A CoWeb network (Upton 2000): the categories on a circle, linked if they go together.

Port of `cobCoin`. It builds the network with `sur_coin`, keeps the links between
categories of different variables whose p-value is below `significance / 2` (they go
together more than expected) or above `1 - significance / 2` (less), and puts the nodes
on a circle.

**Parameters**

- `data` (DataFrame): One row per case.
- `variables` (list of str, optional): The variables to use (all the columns by default).
- `degree` (float): Angle in degrees of the first node of the circle.
- `significance` (float): Two-sided significance level of the links kept.
- `**kwargs`: As in `sur_coin`. By default nodes are coloured and grouped by `variable`, and links are as wide and coloured as their `Haberman` residual (`link_bipolar`).

### `net_corr`

```python
def net_corr(variables: 'pd.DataFrame', weight: 'str | Sequence[float] | None' = None, pairwise: 'bool' = False, sort: 'bool' = False, decreasing: 'bool' = True, frequency: 'bool | None' = None, means: 'bool' = True, method: 'str | Sequence[str]' = 'pearson', criteria: 'str' = 'p', bonferroni: 'bool' = False, min_link: 'float' = 0, max_link: 'float' = inf, sort_links: 'str | None' = None, decreasing_links: 'bool' = True, igraph: 'bool' = False, **netcoin_args)
```

A network of the correlations among quantitative variables (port of `netCorr`).

Every column is a node (with its mean, standard deviation, minimum and maximum) and every
pair a link with the coefficient (`value`), the statistic of its contrast (`z`) and its
p-value (`p`). The links are bipolar: they are coloured by their sign and the contrast is
two-sided, unless you say `link_bipolar=False`.

**Parameters**

- `variables` (DataFrame): One numeric column per variable (at least two) and one row per case.
- `weight` (str or sequence, optional): The weight of each case, or the name of the column that holds it (not for `"kendall"` and `"dcor"`). The cases with a missing weight are left out.
- `pairwise` (bool): Each pair uses the cases that have both values; otherwise only the complete cases are used. The number of cases of each link is then kept in the column `n`.
- `sort, decreasing` (bool): Order the nodes by number of valid cases and, within them, by mean.
- `frequency` (bool, optional): Put the number of valid cases in the nodes (by default only with `pairwise`).
- `means` (bool): Put the mean, standard deviation, minimum and maximum in the nodes.
- `method` (str or list of str): `"pearson"`, `"spearman"`, `"kendall"`, `"partial"` or `"dcor"` (distance correlation, bias corrected); they can be abbreviated. With several, the first gives `value`, `z` and `p` and the others their own columns (`spearman`, `z(spearman)`, `p(spearman)`).
- `criteria` (str): What `min_link` and `max_link` refer to: `"value"`, `"z"`, `"p"` (default), another method, or `"z(method)"`/`"p(method)"`. For a p-value `max_link` is 0.5 (1 with bipolar links) unless you say otherwise; for a coefficient of a bipolar link its absolute value counts.
- `bonferroni` (bool): Divide a p-value `max_link` by the number of pairs.
- `min_link, max_link` (float): Range of `criteria` kept as links.
- `sort_links, decreasing_links`: Column and direction to order the links by.
- `igraph` (bool): Return an `igraph.Graph`.
- `**netcoin_args`: The rest goes to `netCoin` (`nodes`, `name`, `language`, `layout`, `lwidth`...). `layout` can also be `"pc"`: principal components of the correlations. By default `lwidth`, `lweight` and `lcolor` are `"value"`.

**Returns**

- `netCoin or igraph.Graph`

### `d_net_corr`

```python
def d_net_corr(variables: 'pd.DataFrame', nodes: 'pd.DataFrame | None' = None, weight: 'str | Sequence[float] | None' = None, pairwise: 'bool' = False, frequency: 'bool | None' = None, means: 'bool' = True, method: 'str | Sequence[str]' = 'pearson', criteria: 'str' = 'value', bonferroni: 'bool' = False, min_link: 'float' = 0, max_link: 'float' = inf, sort_links: 'str | None' = None, decreasing_links: 'bool' = True, factorial: 'str | Sequence[str]' = 'null', components: 'bool' = True, backcomponents: 'bool' = False, sequence: 'Sequence[float] | None' = None, text_filter: 'Sequence[float]' = (1, 0.99), speed: 'float' = 50, directory: 'str | None' = None, **netcoin_args)
```

A network of correlations that evolves as the threshold of the correlation rises.

Port of `d_netCorr`. It builds `net_corr` and then one frame for every value of
`sequence`: the network of the links with `abs(value)` above it, with the centrality of
each node, the component it belongs to (`C-1` is the largest) and the `ratio` of each
link to the threshold. Only the frames where the number of components grows are kept
(every frame if `components=False`; also those where it falls with `backcomponents`).

**Parameters**

- `variables, nodes, weight, pairwise, frequency, means, method, criteria, bonferroni, min_link, max_link, sort_links, decreasing_links`: As in `net_corr` (`criteria` is `"value"` here).
- `factorial` (str or list of str): Layout from the loadings: `"pc"`, `"nf"`, `"vf"`, `"of"` (see `layout_fact`) or `"null"` for none. Not compatible with `layout`.
- `components, backcomponents` (bool): Which frames are kept (see above).
- `sequence` (sequence of float): The thresholds (0.20 to 1 by 0.01 by default).
- `text_filter` ((float, float)): The links whose `abs(value)` reaches `min(text_filter[0], threshold + text_filter[1])` show their value as text.
- `speed` (float): Speed of the animation (0-100).
- `directory` (str, optional): Write the page there.
- `**netcoin_args`: The rest goes to `net_corr`.

**Returns**

- `netCoin`: The frames merged into one evolving network (see `multigraphCreate`).

**Notes**

The eigenvector centrality of a graph that is not connected is not unique, so its values for the frames with several components can differ from other programs'.

### `log_coin`

```python
def log_coin(data: 'pd.DataFrame', variables: 'Sequence[str] | None' = None, exogenous: 'str | Sequence[str] | None' = None, no_first_cat: 'str | Sequence[str] | None' = None, weight: 'str | Sequence[float] | None' = None, order: 'int' = 2, pairwise: 'bool' = False, twotails: 'bool' = False, pmax: 'float' = 0.05, frequency: 'bool' = False, percentage: 'bool' = False, directed: 'bool' = False, igraph: 'bool' = False, **netcoin_args)
```

A network of the associations among categorical variables from a log-linear model.

Port of `logCoin`. A model with all the interactions up to `order` is fitted to the
table of the variables (Poisson regression on the counts of the cells, with contrasts that
make the effects deviations from the mean). Every category is a node, its estimate being the
main effect; the interactions between two categories whose estimate is significant are links,
and those of three or more categories are joined through a node `u01`, `u02`...

**Parameters**

- `data` (DataFrame): One row per case. Spaces in the names of columns are replaced by `_`.
- `variables` (list of str, optional): The variables (all the columns by default).
- `exogenous` (str or list of str, optional): Variables that are only used to explain the others: their own interactions are not in the model.
- `no_first_cat` (str or list of str, optional): Variables whose first category is not shown as a node (the pivot of the contrasts). At most five variables can have their first category included.
- `weight` (str or sequence, optional): Weight of each case or the column that holds it. The counts are rounded.
- `order` (int): Order of the interactions of the model.
- `twotails` (bool): Two-sided contrast: significant negative links are kept too and the links are bipolar, coloured by their estimate. Otherwise only the positive ones are kept, with the one-sided p-value `Pr(>z)`.
- `pmax` (float): Significance level of the links.
- `frequency, percentage` (bool): Put the number and/or the proportion of cases of each category in the nodes (and size them by it).
- `pairwise, directed`: Accepted for compatibility with R; they have no effect.
- `igraph` (bool): Return an `igraph.Graph`.
- `**netcoin_args`: The rest goes to `netCoin` (`language`, `color`, `lwidth`, `main`...).

**Returns**

- `netCoin or igraph.Graph`

### `glm_coin`

```python
def glm_coin(formulas: 'str | Sequence[str]', data: 'pd.DataFrame', weights: 'str | Sequence[float] | None' = None, pmax: 'float' = 0.05, robust: 'bool' = False, show_arrows: 'bool' = True, color: 'str' = 'variable', lwidth: 'str' = 's.value', lcolor: 'str' = 'estimate', language: 'str' = 'en', igraph: 'bool' = False, **netcoin_args)
```

A network of the effects in one or several regression equations (port of `glmCoin`).

Each line of `formulas` is an equation, `"response ~ term + term, family"`, where the family
is `gaussian` (the default, a linear model), `binomial` (fitted as quasibinomial, with
estimated dispersion), `poisson`, `quasipoisson`, `Gamma`, `inverse.gaussian` or
`quasibinomial`, or `multinomial` (a factor response with several categories: one group of effects per category, the first being the reference). The categories of every factor are nodes, and so are the numeric variables
and the response; the links go from every explanatory category to the response, with the
effect of the category against the mean of its variable weighted by the sample (`estimate`)
and against its first category (`estimate.base`), and from every numeric variable with its
average slope. Effects in binomial and multinomial equations are in percentage points. The standard errors
come from the delta method with the covariance of the model, or the robust HC1 one.

**Parameters**

- `formulas` (str or list of str): The equations (variable names with spaces and other signs are accepted).
- `data` (DataFrame): Text and categorical columns are the factors, numeric columns the quantitative variables.
- `weights` (str or sequence, optional): The weights of the cases, or the column that holds them.
- `pmax` (float): Links whose p-value is not below it are hidden.
- `robust` (bool): Use the HC1 covariance.
- `show_arrows, color, lwidth, lcolor, language`: As in `netCoin`.
- `igraph` (bool): Return an `igraph.Graph`.
- `**netcoin_args`: The rest goes to `netCoin`.

**Returns**

- `netCoin`: With the attribute `equations` (fit of every equation, see `eq_summary`) and `coeffs` (their coefficients, see `eq_coefs`).

### `caring`

```python
def caring(filepath: 'str', arguments: 'Mapping', show: 'bool' = False) -> 'str'
```

Read `filepath` and write the page of networks that `arguments` ask for (see `caring_create_graphs`).

**Parameters**

- `filepath` (str): Data file (`.sav`, `.dta`, `.xlsx` or delimited text); see `caring_read_file`.
- `arguments` (mapping): What to draw and how; see `caring_create_graphs`.
- `show` (bool): Open the page in the browser when it is written.

**Returns**

- `str`: The folder that holds the page.

**Examples**

```python
>>> caring("survey.sav", {"variables": ["gender", "age", "works"], "plot": ["surCoin", "Associations"],
...                       "dir": "out"})  # doctest: +SKIP
```

### `caring_create_graphs`

```python
def caring_create_graphs(data: 'pd.DataFrame', arguments: 'Mapping', show: 'bool' = False) -> 'str'
```

The networks of `arguments["plot"]` in one page; returns the folder.

`arguments` holds the names of R (`variables`, `weight`, `maxL`, `criteria`,
`dichotomies`, `valueDicho`, `metric`, `exogenous`, `chaine`, `family`, `nclustersMin`...) or
python's; `plot` lists, one after the other, the kind of network (`"surCoin"`, `"surScat"`,
`"logCoin2"` to `"logCoin5"`, `"allNet"`, `"glmCoin"`) and the title it gets in the page;
`dir` (or `directory`) is the folder to write.

**Parameters**

- `data` (DataFrame): One row per case.
- `arguments` (mapping): The keys of the R function (`variables`, `weight`, `maxL`, `criteria`, `dichotomies`, `valueDicho`, `metric`, `exogenous`, `chaine`, `family`, `nclustersMin`...) or their python names (`max_l`...). `plot` is a list that alternates the kind of network and its title; `dir` (or `directory`) is the folder to write.
- `show` (bool): Open the page in the browser when it is written.

**Returns**

- `str`: The folder that holds the page.

### `caring_read_file`

```python
def caring_read_file(filepath: 'str') -> 'pd.DataFrame'
```

Read a `.sav`, `.dta`, `.xlsx` or delimited text file as a table of cases.

Text columns become categories, and those that hold only numbers become numeric. A Stata file
takes the label of each variable as its name and labelled values become categories.

**Parameters**

- `filepath` (str): `.sav` (needs `pyreadstat`), `.dta`, `.xlsx` (needs `openpyxl`), or a delimited text file. Both are in the `files` extra.

**Returns**

- `DataFrame`: One row per case.

## Scattergram and clusters

### `sur_scat`

```python
def sur_scat(data: 'pd.DataFrame', variables: 'Sequence[str] | None' = None, active: 'Sequence[str] | None' = None, columns: 'Sequence[str] | None' = None, patterns: 'bool' = False, v_patterns: 'Sequence[str] | None' = None, weight: 'str | Sequence[float] | None' = None, seed: 'int | None' = 2020, max_n: 'int' = 2000, jitter: 'float | bool' = 0, type: 'str | None' = None, xaxis=None, yaxis=None, scale_axes: 'bool' = True, nfactors: 'int' = 2, crit_factors: 'float' = 0, layouts=None, suffix: 'str | None' = None, align_to=None, nclusters=None, cluster_on: 'str' = 'factors', scale_clusters: 'bool | None' = None, nstart: 'int' = 25, cluster_index: 'str' = 'ch', sort_clusters: 'bool' = True, clusters=None, force_execution: 'bool' = False, **netcoin_args)
```

A scattergram of the cases of a survey, coloured by groups (port of `surScat`).

The cases (or the patterns of answers, with `patterns`) are drawn on the plane of a
multiple correspondence analysis (qualitative active variables) or of a principal components
analysis (quantitative ones), or on observed variables, and grouped with k-means. The result
is a `netCoin` object without links, whose nodes hold the variables, the groups and
the coordinates; its attribute `scat` keeps what the later steps use (the case-to-pattern
map, the shares of every group in every pattern...), and `clusters` and `cluster_means` the
centroids of the groups.

**Parameters**

- `data` (DataFrame): One row per case.
- `variables, active` (list of str, optional): The variables shown in the nodes and those the analysis is made of (`active` defaults to `variables`). A variable used on an axis joins `variables`.
- `columns` (list of str, optional): Columns of the node table to put first (`"_vpattern"`, `"_actives"`, `"_others"` and `"_nclusters"` stand for those blocks).
- `patterns, v_patterns` (bool, list of str): Merge the cases that give the same answers to `v_patterns` into one node (a node of patterns is sized by its number of cases). Stating `v_patterns` implies `patterns`.
- `weight` (str or sequence, optional): Weight of each case.
- `seed` (int, optional): Seed of the random steps (k-means, samples, jitter).
- `max_n` (int): At most this many nodes are drawn, a random sample of them if there are more.
- `jitter` (float or bool): A little random spread of the coordinates, as a fraction of their range (`True` is 0.02).
- `type` ({"mca", "pca"}, optional): Kind of analysis; by default `"pca"` if all the active variables are numbers.
- `xaxis, yaxis` (int, str or (axis, bool), optional): What is drawn on each axis: a factorial axis (`1`, `"F1"`, `"PC2"`) or a numeric (or ordered) variable, optionally followed by whether to standardize it. With a variable on both axes there are no groups by default.
- `scale_axes` (bool): Standardize the factorial axes.
- `nfactors, crit_factors` (int, float): Factors used to cluster: those of the first `nfactors` whose eigenvalue exceeds `crit_factors` times the mean eigenvalue.
- `nclusters` (int, "auto" or list): Number(s) of k-means groups, one column each (the last stated colours the graph). 0 draws no groups. `"auto"` chooses among the numbers stated next to it (2 to 10) with `cluster_index`.
- `cluster_on` ({"factors", "variables"}): What the k-means works on.
- `scale_clusters` (bool, optional): Standardize the clustering space (by default yes, except indicators of categories).
- `nstart` (int): Random starts of the k-means (the best is kept).
- `cluster_index` ({"ch", "silhouette", "bic"}): Index that chooses the number of groups: Calinski-Harabasz, silhouette, or BIC of a latent model (LPA for `type="pca"`, LCA for `"mca"`, lower is better); with `"auto"` next to a 1 the BIC also weighs the absence of groups.
- `sort_clusters` (bool): Number the groups from left to right on the horizontal axis.
- `layouts` (bool or list of str, optional): Planes added to the main one: the reductions of `looking4clusters` (`"pca"`, `"tsne"`, `"cmds"`, `"nmf"`, `"umap"`) and the rotations of the components (`"varimax"`, `"quartimax"`, `"equamax"`, only with `type="pca"`); `True` or `"all"` asks for all of them.
- `clusters` (bool or list of str, optional): Groups added to the main ones, by number of groups in `nclusters`: `"kmeans"`, `"pam"`, `"hclust"` (of `looking4clusters`, on the active variables as they come) and the latent models `"lpa"` / `"lca"`; `True` asks for all of them (the latent model that suits `type`).
- `suffix, align_to` (str, str or int, optional): With `suffix`, every plane added is also shown aligned by Procrustes to the main plane (or to `align_to`, a plane's name or position), under its name followed by the suffix.
- `force_execution` (bool): Compute `pam`, `hclust` and `cmds` even beyond 5000 cases.
- `**netcoin_args`: The rest goes to `netCoin` (`language`, `main`, `color`...).

**Returns**

- `netCoin`

### `add_axes`

```python
def add_axes(net, axes, name: 'str | None' = None, which=None, weight=None, axes_labels=None)
```

Add one or several planes to a scattergram (`addAxes`).

`axes` is a table with an even number of columns (a pair per plane), a list of them, a dict or a
`Looking4Clusters`, whose reductions become planes. `name` names the plane (a prefix if
several); `weight` weighs the cases when they are collapsed into patterns. The coordinates already
drawn become the first plane. `axes_labels`: `None` (the default) lets every plane keep its own
axis titles, a pair puts the same titles on every plane, and `"keep"` leaves what the object had.

**Parameters**

- `net` (netCoin): A scattergram made by `sur_scat`.
- `axes` (DataFrame, list, dict or Looking4Clusters): The coordinates, one pair of columns per plane.
- `name` (str, optional): Name of the plane (a prefix when there are several).
- `which` (optional): For a `Looking4Clusters`, which reductions to add.
- `weight` (sequence, optional): Weight of each case, used when the cases are collapsed into patterns.
- `axes_labels` (None, "keep" or pair of str): Titles of the axes; see above.

**Returns**

- `netCoin`: A modified copy; the planes are in `net.layouts`.

### `add_clusters`

```python
def add_clusters(net, clusters, name: 'str | None' = None, sort: 'bool' = True, weight=None, which=None, max_groups: 'int | None' = 25, rest: 'int | None' = None)
```

Add one or several cluster columns to a scattergram (`addClusters`).

`clusters` is a vector (one group per node, or per case when the nodes are patterns), a table of
one clusterization per column, a fitted scikit-learn clusterer or a `Looking4Clusters`
(every clusterization it holds, or the one `which` chooses). The groups are renumbered from left
to right on the first axis when `sort`. `max_groups` skips clusterizations with more groups
(a continuous variable passed by mistake); `rest` is the number of a residual group.

**Parameters**

- `net` (netCoin): A scattergram made by `sur_scat`.
- `clusters` (sequence, DataFrame, clusterer or Looking4Clusters): The groups.
- `name` (str, optional): Name of the column (a prefix when there are several).
- `sort` (bool): Renumber the groups from left to right on the first axis.
- `weight` (sequence, optional): Weight of each node or case, for the centroids and the order of the groups.
- `which` (optional): For a `Looking4Clusters`, which clusterizations to add.
- `max_groups` (int, optional): Skip clusterizations with more groups than this.
- `rest` (int, optional): Number of the group that gathers what is left over.

**Returns**

- `netCoin`: A modified copy; the columns are added to `net.nodes`, and the centroids of every group on every plane to `net.clusters`.

### `replace_clusters`

```python
def replace_clusters(net, clusters, name: 'str | None' = None, sort: 'bool' = True, weight=None, which=None, max_groups: 'int | None' = 25)
```

Replace the cluster columns of a scattergram by new ones (`replaceClusters`).

**Parameters**

- `net` (netCoin): A scattergram made by `sur_scat`.
- `clusters, name, sort, weight, which, max_groups`: As in `add_clusters`; the cluster columns the object held are removed first.

**Returns**

- `netCoin`: A modified copy.

### `layouts_align`

```python
def layouts_align(net, align: 'str' = 'procrustes', to=1, gpa: 'bool' = False, weights=None, which=None, suffix: 'str | None' = None, tol: 'float' = 1e-10, maxit: 'int' = 100, who: 'str' = 'layoutsAlign')
```

Lay the planes of a scattergram over one another (`layoutsAlign`).

`align` is `"procrustes"` (centre, equal size and best rotation or reflection), `"cloud"`
(centre and size only), `"canvas"` (each axis to [0, 1]) or `"none"`. `to` is the plane the
rest are brought to (name or 1-based index); with `gpa` every plane is rotated towards the mean of
them all. With `suffix` the originals stay and the aligned planes are added, right after their
own, under their name and the suffix; otherwise they take the place of the originals.

**Parameters**

- `net` (netCoin): A scattergram with more than one plane.
- `align` ({"procrustes", "cloud", "canvas", "none"}): How the planes are laid over one another.
- `to` (int or str): The plane the rest are brought to (1-based position or name).
- `gpa` (bool): Rotate every plane towards the mean of them all (generalized Procrustes).
- `weights` (sequence or str, optional): Weight of each node in the fit.
- `which` (sequence, optional): The planes to align (all by default).
- `suffix` (str, optional): Keep the originals and add the aligned planes under their name followed by `suffix`.
- `tol, maxit` (float, int): Convergence of the generalized Procrustes iteration.

**Returns**

- `netCoin`: A modified copy.

### `add_communities`

```python
def add_communities(net, data: 'pd.DataFrame', variables=None, community='lo', scale: 'bool' = True, method: 'str' = 'pearson', criteria: 'str' = 'value', min_link: 'float' = 0.95, max_link: 'float' = inf, pairwise: 'bool' = False, min_size: 'float' = 0.01, max_communities: 'int | None' = None, name: 'str | None' = None, sort: 'bool' = True, weight=None, max_groups: 'int | None' = 25, max_nodes: 'float | None' = None, seed: 'int | None' = 2020, weighted: 'bool' = True)
```

Add a cluster column per community-detection method, from the network of correlations between the cases.

Port of `addCommunities`. The cases are correlated with one another over `variables` (the transposed
reading of the table); the communities of that network are a grouping of the cases, and those smaller than
`min_size` (a proportion below 1, a count otherwise) are gathered into a residual group.
`community` takes codes (`"lo"`) or names (`"louvain"`): edgeBetweenness, fastGreedy, labelProp,
leadingEigen, louvain, optimal, spinglass, walktrap. Needs `python-igraph`.

The links of the network are weighed by their correlation (what the R code means to do). R's `toIgraph`
does not carry that weight over from the network of `netCorr`, so R finds the communities of the
unweighted network; `weighted=False` does the same here.

**Parameters**

- `net` (netCoin): A scattergram made by `sur_scat`.
- `data` (DataFrame): The table the scattergram was drawn from.
- `variables` (list of str, optional): The variables the cases are compared over (all by default); numeric or ordered categorical.
- `community` (str or list of str): The methods, by code or name.
- `scale` (bool): Standardize each variable first.
- `method, criteria, min_link, max_link, pairwise`: Passed to `net_corr`, which builds the network of the cases (links need `|r| >= min_link` by default 0.95).
- `min_size` (float): Communities smaller than this (a proportion if below 1, a count otherwise) go into a residual group.
- `max_communities` (int, optional): At most this many groups, the residual one among them.
- `name` (str, optional): Name of the column (`name.method(k)` with several methods).
- `sort, weight, max_groups`: As in `add_clusters`.
- `max_nodes` (float, optional): Override the size each method can bear (`optimal` 100, `edgeBetweenness` and `spinglass` 2000).
- `seed` (int, optional): Seed of the methods that draw at random (louvain, labelProp, spinglass).
- `weighted` (bool): Weigh the links by their correlation; `False` reproduces R.

**Returns**

- `netCoin`: A modified copy with one cluster column per method, named `method(k)`.

### `cluster_match`

```python
def cluster_match(x, model, columns=None, weights=None, align: 'bool' = True, weighted: 'bool' = True, modal: 'bool' = True, sort=None, digits: 'int' = 3, model_name: 'str | None' = None) -> 'pd.DataFrame'
```

Compare a model classification with the others of a scattergram, category by category.

Port of `clusterMatch`. `model` is the name of a node column or a vector (one per node, or per case
when the nodes are patterns); `columns` are the classifications compared (the groups of the object by
default; only those with as many categories as the model); the categories are matched one to one
(`align`) so that the diagonal of the table is as large as possible. `weighted` adds the means of
the specificity and F1 weighted by the prevalence; `modal` compares, where a pattern spans several
groups, its modal group. The result has one row per classification; its `attrs` hold `model`,
`unit`, `weights`, `spanning`, `tables`, `byClass`, `mapping` and `ties` (the
classifications for which several matchings were equally good: R may pick another of them).

**Parameters**

- `x` (netCoin): A scattergram made by `sur_scat`.
- `model` (str or sequence): A node column, or a vector (per node, or per case when the nodes are patterns).
- `columns` (list of str, optional): The classifications compared.
- `weights` (str or sequence, optional): Weight of each node (a column or a vector).
- `align, weighted, modal` (bool): See above.
- `sort` (str or list of str, optional): Column(s) to sort by, from the largest down; a `"-"` prefix sorts from the smallest up.
- `digits` (int): Decimals of the result.
- `model_name` (str, optional): Name of the model when it is a vector.

**Returns**

- `DataFrame`: One row per classification: `Classification`, `k`, `N`, `Accuracy`, `Kappa`, `Sensitivity`, `Precision`, `F1` and `Specificity` (macro means; `.w` are the weighted ones).

### `cluster_agree`

```python
def cluster_agree(x, model, columns=None, weights=None, modal: 'bool' = True, sort=None, digits: 'int' = 3, model_name: 'str | None' = None) -> 'pd.DataFrame'
```

The same comparison as agreement between partitions (ARI, NMI, VI): no matching, any number of categories.

Port of `clusterAgree`. The `k` column is that of each classification (the model's is in
`attrs['kModel']`).

**Parameters**

- `x` (netCoin): A scattergram made by `sur_scat`.
- `model` (str or sequence): A node column, or a vector (per node, or per case when the nodes are patterns).
- `columns, weights, modal, sort, digits, model_name`: As in `cluster_match`.

**Returns**

- `DataFrame`: One row per classification: `Classification`, `k`, `N`, `ARI`, `NMI`, `VI` and `VI.n` (the variation of information, raw and normalized).

### `looking4clusters`

```python
def looking4clusters(data, groups=None, components: 'bool' = False, running_all: 'bool' = True, distance: 'str' = 'euclidean', agglomeration: 'str' = 'complete', selectedk: 'int | None' = None, perplex: 'float' = 30, max_iter: 'int' = 1000, force_execution: 'bool' = False, seed: 'int | None' = None) -> 'Looking4Clusters'
```

Look for groups of cases and for two-dimensional views of them.

Port of `looking4clusters`. With `running_all` it runs k-means, PAM and hierarchical
clustering for 2 to 10 groups (or `selectedk` plus or minus five) and the reductions PCA,
t-SNE, classical MDS (on Canberra distances), NMF and UMAP; each needs its own library and
is skipped with a warning if it is missing. The best number of groups of every method is
stored in `optim` (Calinski-Harabasz).

**Parameters**

- `data` (DataFrame or 2-d array): One row per case, one numeric column per variable.
- `groups` (sequence, optional): Groups of your own, one per case.
- `components` (bool): Work on the principal components of `log2(data + 1)` that explain the data instead of on the variables (at least five cases and five variables).
- `running_all` (bool): Run the methods; otherwise only your `groups` are stored.
- `distance` (str): `"euclidean"`, `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"` for PAM and hierarchical clustering.
- `agglomeration` (str): `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`, `"mcquitty"`, `"median"` or `"centroid"`.
- `selectedk` (int, optional): The number of groups you prefer.
- `perplex, max_iter`: Of t-SNE.
- `force_execution` (bool): Run the slow methods also with more than 5000 cases (or variables).
- `seed` (int, optional): Seed of the methods that start at random (k-means, t-SNE, NMF, UMAP).

**Returns**

- `Looking4Clusters`

### `Looking4Clusters`

```python
class Looking4Clusters(data: 'np.ndarray', samples: 'list[str]', variables: 'list[str]', options: 'dict' = <factory>, clusters: 'dict' = <factory>, reductions: 'dict' = <factory>, optim: 'dict' = <factory>, pctvar: 'dict' = <factory>) -> None
```

Result of `looking4clusters`: the data, its clusterings and its reductions.

`clusters` maps the name of a method (`kmeans`, `pam`, `hclust`...) to a categorical
Series, or to a DataFrame with one column per number of groups (`levels_2`, `levels_3`...);
`optim` holds the number of groups the method found best. `reductions` maps a name
(`pca`, `tsne`, `cmds`...) to a DataFrame of two coordinates, and `pctvar` the
percentage of the variance of those that have it.

**Attributes**

- `data` (ndarray): The values clustered (samples by variables).
- `samples, variables` (list of str): Their names.
- `options` (dict): Options kept for the page (`groupStatsBy`, `myGroups`).
- `clusters` (dict): Name of the method to a categorical Series, or a DataFrame with a column per number of groups.
- `reductions` (dict): Name of the reduction to a DataFrame of two coordinates.
- `optim` (dict): The number of groups each method found best.
- `pctvar` (dict): Percentage of variance of the first two axes of the reductions that have it.

## Layouts

### `layout_mca`

```python
def layout_mca(incidences: 'pd.DataFrame') -> 'pd.DataFrame'
```

Coordinates of the variables from a simple correspondence analysis.

Port of `layoutMCA`: classical correspondence analysis of the table of
incidences completed with its complement. Returns principal coordinates
`F1` and `F2`, one row per variable.

**Parameters**

- `incidences` (DataFrame): One row per case, one 0/1 column per variable.

**Returns**

- `DataFrame`: `F1` and `F2`, one row per variable.

### `layout_pca`

```python
def layout_pca(coincidences: 'Coin | pd.DataFrame') -> 'pd.DataFrame'
```

Coordinates of the variables from the eigenvectors of their Pearson correlations.

Port of `layoutPCA`. Takes a `Coin` (or a correlation matrix) and
returns `F1` and `F2`, one row per variable.

**Parameters**

- `coincidences` (Coin or DataFrame): The coincidences, or a correlation matrix.

**Returns**

- `DataFrame`: `F1` and `F2`, one row per variable.

### `layout_fact`

```python
def layout_fact(data: 'pd.DataFrame', method='pc')
```

Coordinates of the variables from their loadings on two components or factors.

Port of `layoutFact`. `method` is any of `"pc"` (principal components of the
correlations), `"nf"` (maximum likelihood factors, not rotated), `"vf"` (varimax) and
`"of"` (oblimin), or `"null"` for no layout. The cases with a missing value are left
out. One method gives a DataFrame (`Comp.1`, `Comp.2` or `Factor1`, `Factor2`, one
row per variable); several give a dict of them named as in R. The sign of each axis is
arbitrary, as it is in R.

**Parameters**

- `data` (DataFrame): One numeric column per variable.
- `method` (str or list of str): `"pc"`, `"nf"`, `"vf"`, `"of"` or `"null"`.

**Returns**

- `DataFrame or dict of DataFrame`

### `layout_circle`

```python
def layout_circle(nodes: 'pd.DataFrame', items: 'Sequence | None' = None, deg: 'float' = 0, name: 'str | None' = None) -> 'pd.DataFrame'
```

Place nodes on a circle (rD3plot's `layoutCircle`).

**Parameters**

- `nodes` (DataFrame): The table of nodes (one row each).
- `items` (names or 0-based positions, optional): The nodes to place (all by default); the rest get `NaN`. They are spread evenly in the order given, the first at angle `deg`.
- `deg` (float): Angle in degrees of the first node.
- `name` (str, optional): Column with the names of the nodes (otherwise `"1"`, `"2"`...).

**Returns**

- `DataFrame`: `x` and `y` (rounded to 3 decimals), one row per node.

### `layout_grid`

```python
def layout_grid(nodes: 'pd.DataFrame', string: 'str', name: 'str | None' = None, byrow: 'bool' = False) -> 'pd.DataFrame'
```

Place nodes on a grid described by a text (rD3plot's `layoutGrid`).

`string` lists the columns of the grid separated by `"."` and the nodes of each
column from top to bottom separated by `","`; `";"` leaves an empty cell and a
column that ends with `","` leaves one at the bottom. A column with fewer nodes
is centred. For example `"a,b.c.d,e,f"`. With `byrow` the grid is read by rows.

Returns a DataFrame `x`, `y` in [0, 1] (`NaN` for the nodes that are not in the
text).

**Parameters**

- `nodes` (DataFrame): The nodes (their `name` column or `name` is matched against the text).
- `string` (str): The description of the grid.
- `name` (str, optional): Column of `nodes` with the names.
- `byrow` (bool): Read the grid by rows.

## igraph and files

### `to_igraph`

```python
def to_igraph(net)
```

The network of a `netCoin` object as an `igraph.Graph`.

Port of `toIgraph`. The vertices carry the node table (`name`, `label`,
`size`, `color`, `shape`... taken from the columns the object uses for them;
`fx`/`fy` become `x`/`y`) and the edges the link table (`width`,
`weight`, `color`, `label` likewise). The graph is directed if
`show_arrows` was set.

**Parameters**

- `net` (netCoin): The network.

**Returns**

- `igraph.Graph`

### `from_igraph`

```python
def from_igraph(graph, **kwargs)
```

A `netCoin` object from an `igraph.Graph`.

Port of `fromIgraph`. The vertex attributes become node columns and the edge
attributes link columns; those called `label`, `label.cex`, `size`, `color`
and `shape` (vertices) or `width`, `weight`, `color` and `label` (edges)
are used for the corresponding argument unless you give it. Other keyword
arguments go to `netCoin`.

**Parameters**

- `graph` (igraph.Graph): The graph.
- `**kwargs`: Arguments of `netCoin`.

**Returns**

- `netCoin`

### `cal_centr`

```python
def cal_centr(graph, measures: 'Sequence[str] | str' = ('degree', 'wdegree', 'closeness', 'betweenness', 'eigen'), order: 'str' = '') -> 'pd.DataFrame'
```

Centrality of every node (port of `calCentr`).

**Parameters**

- `graph` (netCoin or igraph.Graph)
- `measures` (str or list of str): `"degree"`, `"wdegree"` (strength; weights in the edge attribute `N`, as in R), `"closeness"`, `"betweenness"`, `"eigen"` or `"all"` (only the first letter counts). In directed graphs `indegree` and `outdegree` come with `degree`.
- `order` (str): Sort the result by this column, decreasing.

**Returns**

- `DataFrame`: One row per node: `nodes` and a column per measure. (R also returns the centralization of each measure; that is not computed here.)

### `save_pajek`

```python
def save_pajek(net, file: 'str' = 'file.net', arcs: 'Sequence[str] | None' = None, edges: 'Sequence[str] | None' = None, partitions: 'Sequence[str] | None' = None, vectors: 'Sequence[str] | None' = None, encoding: 'str' = 'utf-8') -> 'str'
```

Write a netCoin object as a Pajek file.

Port of `savePajek`.

**Parameters**

- `net` (netCoin)
- `file` (str): Name of the file; `.net` is added if it has no extension, and a project file `.paj` is written instead when there are partitions or vectors.
- `arcs, edges` (list of str): Link columns written as weights of directed (`*Arcs`) and undirected (`*Edges`) links, one block for each.
- `partitions, vectors` (list of str): Node columns written as partitions (categories coded 1, 2...) and as vectors (missing values as 0).

**Returns**

- `str`: The path written.

### `save_graphml`

```python
def save_graphml(net, file: 'str' = 'netCoin.graphml') -> 'str'
```

Write a netCoin object as GraphML (port of `saveGhml`). Needs python-igraph.

Categorical columns are written as text, since GraphML only has numbers, text and
logical values.

**Parameters**

- `net` (netCoin): The network.
- `file` (str): Path of the file.

**Returns**

- `str`: The path written.
