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
