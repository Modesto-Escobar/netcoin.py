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
| `dichotomize(..., sep = "C")` | Returns meaningless columns (e.g. `V3`, `azul;rojo;verde`). |
| `edgeList(..., criteria = "Fisher")` | The p-value cap (`max = .5`) is applied for `"Z"` and `"Hyp"` but not for `"Fisher"`. |
| `sim(coin, "L")` alone | Returns a single matrix instead of `conf.L`, `confidence`, `conf.U`. |
