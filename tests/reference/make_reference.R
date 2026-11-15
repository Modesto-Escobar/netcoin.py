# Generates the reference values that the Python tests compare against.
#
# Run from the repository root with the R installation that has netCoin:
#
#   Rscript tests/reference/make_reference.R
#
# Inputs  -> tests/reference/inputs/*.csv     (datasets shipped with netCoin
#                                              plus small synthetic ones)
# Outputs -> tests/reference/expected/*.csv   (what netCoin returns)
#         -> tests/reference/manifest.json    (versions and the R call of each case)
#
# The CSV files are committed, so running the Python tests does not need R.

suppressMessages({
  library(netCoin)
  library(jsonlite)
})

root <- file.path("tests", "reference")
dir_in  <- file.path(root, "inputs")
dir_out <- file.path(root, "expected")
unlink(c(dir_in, dir_out), recursive = TRUE)
dir.create(dir_in,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

manifest <- list(
  R = R.version.string,
  netCoin = as.character(packageVersion("netCoin")),
  rD3plot = as.character(packageVersion("rD3plot")),
  generated = format(Sys.time(), "%Y-%m-%d"),
  collate = Sys.getlocale("LC_COLLATE"),
  cases = list()
)

write_csv <- function(x, name, dir, row_names = FALSE) {
  utils::write.csv(x, file.path(dir, paste0(name, ".csv")),
                   row.names = row_names, fileEncoding = "UTF-8", na = "NA")
}
add_case <- function(id, call, files) {
  manifest$cases[[id]] <<- list(call = call, files = files)
}
# a square matrix keeps its row names in the first column
write_matrix <- function(m, name) {
  m <- as.matrix(m)
  df <- as.data.frame(m, check.names = FALSE, stringsAsFactors = FALSE)
  df <- cbind(.row = rownames(m), df)
  write_csv(df, name, dir_out)
}

# ---------------------------------------------------------------- inputs ----
data(Galapagos, dice, ess, package = "netCoin", envir = environment())
galapagos <- Galapagos
inc_dice  <- dice[, 2:11]

# incidences with missing values, for pairwise = TRUE
galapagos_na <- galapagos
galapagos_na[cbind(c(1, 3, 5, 8, 12, 15), c(2, 5, 7, 7, 1, 10))] <- NA

ess_chr <- data.frame(lapply(ess[, 1:4], as.character), stringsAsFactors = FALSE)

txt <- data.frame(
  id = 1:8,
  keywords = c("red;blue", "blue; green", "red", "Green;red;blue", NA, "", "red;red", "azul;rojo;verde"),
  kw2      = c("red;blue", "blue;green", "red", "Green;red;blue", NA, "", "red;red", "azul;rojo;verde"),
  colour   = c("rojo", "azul", NA, "verde", "rojo", "azul", "rojo", "rojo"),
  grp      = c("A", "B", NA, "A", "B", "A", NA, "B"),
  sentence = c("The cat and the dog", "a Dog, a cat!", "the bird", "Cat dog bird", "", NA, "The Dog", "dog dog cat"),
  stringsAsFactors = FALSE
)

weights <- c(1, 2, 1, 0.5, 1, 3, 1, 1, 2, 1, 1, 0.5, 1, 1, 2, 1, 1)

write_csv(galapagos,    "galapagos",    dir_in)
write_csv(galapagos_na, "galapagos_na", dir_in)
write_csv(inc_dice,     "dice",         dir_in)
write_csv(ess_chr,      "ess",          dir_in)
write_csv(txt,          "txt",          dir_in)
write_csv(data.frame(weight = weights), "weights", dir_in)

C_gal  <- coin(galapagos)
C_dice <- coin(inc_dice)

# ------------------------------------------------------------------ coin ----
coin_case <- function(id, call_txt, obj) {
  write_matrix(obj[, ], paste0("coin_", id))
  write_csv(data.frame(n = attr(obj, "n"),
                       m = if (!is.null(attr(obj, "m"))) 1 else 0), paste0("coin_", id, "_attrs"), dir_out)
  if (!is.null(attr(obj, "m"))) {
    write_matrix(attr(obj, "m"), paste0("coin_", id, "_m"))
    write_matrix(attr(obj, "x"), paste0("coin_", id, "_x"))
  }
  add_case(paste0("coin_", id), call_txt, paste0("coin_", id))
}
coin_case("galapagos",         "coin(galapagos)",                          C_gal)
coin_case("dice",              "coin(dice[, 2:11])",                       C_dice)
coin_case("galapagos_min3",    "coin(galapagos, minimum = 3)",             coin(galapagos, minimum = 3))
coin_case("galapagos_sort",    "coin(galapagos, sort = TRUE)",             coin(galapagos, sort = TRUE))
coin_case("galapagos_sort_inc","coin(galapagos, sort = TRUE, decreasing = FALSE)",
          coin(galapagos, sort = TRUE, decreasing = FALSE))
coin_case("galapagos_weight",  "coin(galapagos, weight = weights)",        coin(galapagos, weight = weights))
coin_case("galapagos_total",   "coin(galapagos, total = TRUE)",            coin(galapagos, total = TRUE))
coin_case("galapagos_na_pairwise", "coin(galapagos_na, pairwise = TRUE)",  coin(galapagos_na, pairwise = TRUE))
coin_case("galapagos_na_listwise", "coin(galapagos_na)",                   coin(galapagos_na))

# ------------------------------------------------------------------- sim ----
sim_procs <- c("matching", "Rogers", "Gower", "Sneath", "Anderberg", "Jaccard", "dice",
               "antiDice", "Ochiai", "Kulczynski", "Hamann", "Yule", "Pearson",
               "tetrachoric", "odds", "Russell", "expected", "confidence", "Haberman", "Z",
               "Fisher", "coincidences", "relative", "sConditional", "tConditional",
               "c.conditional", "c.probable")
# method codes understood by netCoin::sim (the words above are not all valid)
sim_codes <- c(matching = "M", Rogers = "T", Gower = "G", Sneath = "S", Anderberg = "B",
               Jaccard = "J", dice = "D", antiDice = "A", Ochiai = "O", Kulczynski = "K",
               Hamann = "N", Yule = "Y", Pearson = "P", tetrachoric = "V", odds = "C",
               Russell = "R", expected = "E", confidence = "L", Haberman = "H", Z = "Z",
               Fisher = "W", coincidences = "F", relative = "X", sConditional = "I",
               tConditional = "0", c.conditional = "U", c.probable = "Q")

# runner(codes) calls sim() with those method codes
sim_case <- function(ds, proc, runner, call_txt, suffix = "") {
  code <- sim_codes[[proc]]
  # sim(C, "L") alone returns a single matrix, so it is asked together with "J"
  codes <- if (code == "L") c("L", "J") else code
  r <- runner(codes)
  if (code == "L") r <- r[c("conf.L", "confidence", "conf.U")]
  if (!is.list(r)) r <- setNames(list(r), proc)
  files <- character()
  for (nm in names(r)) {
    f <- paste0("sim_", ds, "_", gsub("[^A-Za-z0-9.]", "_", nm), suffix)
    write_matrix(r[[nm]], f)
    files <- c(files, f)
  }
  add_case(paste0("sim_", ds, "_", proc, suffix), sprintf(call_txt, code), files)
}
for (p in sim_procs) {
  sim_case("galapagos", p, function(cd) sim(C_gal, cd),   "sim(coin(galapagos), \"%s\")")
  sim_case("dice",      p, function(cd) sim(C_dice, cd),  "sim(coin(dice), \"%s\")")
}
sim_case("galapagos", "confidence", function(cd) sim(C_gal, cd, level = .90),
         "sim(coin(galapagos), \"%s\", level = .90)", suffix = "_l90")
# sim() on the incidence data frame (computes coin internally)
r <- sim(galapagos, "J", minimum = 3)
write_matrix(r, "sim_galapagos_Jaccard_min3")
add_case("sim_galapagos_Jaccard_min3", "sim(galapagos, \"J\", minimum = 3)", "sim_galapagos_Jaccard_min3")
# distances
sim_case("galapagos", "Jaccard", function(cd) sim(C_gal, cd, distance = TRUE) , "sim(coin(galapagos), \"%s\", distance = TRUE)", suffix = "_dist")
# missing values treated pairwise, and weighted cases
C_pw <- coin(galapagos_na, pairwise = TRUE)
for (p in c("Jaccard", "Haberman", "Z", "expected", "confidence", "Fisher", "c.conditional", "Pearson")) {
  sim_case("galapagos_na", p, function(cd) sim(C_pw, cd, pairwise = TRUE),
           "sim(coin(galapagos_na, pairwise = TRUE), \"%s\", pairwise = TRUE)", suffix = "_pairwise")
}
for (p in c("Jaccard", "Haberman", "Z", "expected", "confidence", "Fisher")) {
  sim_case("galapagos", p, function(cd) sim(galapagos, cd, weight = weights),
           "sim(galapagos, \"%s\", weight = weights)", suffix = "_weight")
}

# -------------------------------------------------------------- edgeList ----
edge_case <- function(id, call_txt, E) {
  if (is.null(E)) {
    write_csv(data.frame(), paste0("edge_", id), dir_out)
  } else {
    E[] <- lapply(E, function(x) if (is.factor(x)) as.character(x) else x)
    write_csv(E, paste0("edge_", id), dir_out)
  }
  add_case(paste0("edge_", id), call_txt, paste0("edge_", id))
}
edge_case("default_galapagos", "edgeList(coin(galapagos))", edgeList(C_gal))
edge_case("jaccard_dice_min", "edgeList(coin(galapagos), c(\"Jaccard\",\"dice\"), criteria = \"Jaccard\", min = .3)",
          edgeList(C_gal, c("Jaccard", "dice"), criteria = "Jaccard", min = .3))
edge_case("dice_sorted", "edgeList(coin(dice), c(\"Haberman\",\"Z\"), criteria = \"Haberman\", min = 2, sort = \"Haberman\")",
          edgeList(C_dice, c("Haberman", "Z"), criteria = "Haberman", min = 2, sort = "Haberman"))
edge_case("bonferroni", "edgeList(coin(galapagos), \"Haberman\", Bonferroni = TRUE)",
          edgeList(C_gal, "Haberman", Bonferroni = TRUE))
edge_case("support", "edgeList(coin(galapagos), c(\"Haberman\",\"coincidences\",\"expected\",\"confidence\"), support = 3)",
          edgeList(C_gal, c("Haberman", "coincidences", "expected", "confidence"), support = 3))
edge_case("directed", "edgeList(coin(galapagos), \"Haberman\", directed = TRUE, diagonal = TRUE)",
          edgeList(C_gal, "Haberman", directed = TRUE, diagonal = TRUE))
edge_case("conditional", "edgeList(coin(galapagos), c(\"CC\",\"CP\"), criteria = \"Z\")",
          edgeList(C_gal, c("CC", "CP"), criteria = "Z"))
# "Hyp" (not "Fisher") so that the p-value cap max = .5 applies, as for "Z"
edge_case("fisher", "edgeList(coin(galapagos), \"Hyp\", criteria = \"Hyp\")",
          edgeList(C_gal, "Hyp", criteria = "Hyp"))
edge_case("level90", "edgeList(coin(dice), \"confidence\", criteria = \"Haberman\", min = 3, level = .90)",
          edgeList(C_dice, "confidence", criteria = "Haberman", min = 3, level = .90))
M <- as.matrix(C_gal[, ])
edge_case("shape", "edgeList(as.matrix(coin(galapagos)), \"shape\")", edgeList(M, "shape"))
adj <- data.frame(a = c("odd", "even", "odd"), b = c("1", "2", "3"), c = c("5", "4", NA),
                  stringsAsFactors = FALSE)
write_csv(adj, "adjacency", dir_in)
# in R an adjacency list is requested with "shape" on a data.frame (no "tree" procedure)
edge_case("adjacency", "edgeList(adjacency, \"shape\")", edgeList(adj, "shape"))

# ---------------------------------------------------------- expectedList ----
edge_case("expected_list", "expectedList(coin(galapagos))", expectedList(C_gal))
edge_case("expected_list_min2", "expectedList(coin(galapagos), min = 2)", expectedList(C_gal, min = 2))
# NOTE: expectedList(C_gal, confidence = TRUE) fails in netCoin 2.1.49: it calls
# edgeList(<matrix>, "shape", min = -Inf), which rewrites min to 1 and drops the
# confidences below 1, so the final cbind() has 78 and 64 rows. Not used as reference.

# ------------------------------------------------------------ dichotomize ----
di_case <- function(id, call_txt, expr) {
  r <- as.data.frame(expr, check.names = FALSE, stringsAsFactors = FALSE)
  r[] <- lapply(r, function(x) if (is.factor(x)) as.character(x) else x)
  write_csv(r, paste0("dicho_", id), dir_out)
  add_case(paste0("dicho_", id), call_txt, paste0("dicho_", id))
}
di_case("ess_age",        "dichotomize(ess, \"Age\")",                      dichotomize(ess_chr, "Age"))
di_case("ess_two",        "dichotomize(ess, c(\"Gender\",\"Age\"))",        dichotomize(ess_chr, c("Gender", "Age")))
di_case("kw",             "dichotomize(txt, \"kw2\", sep = \";\")",    dichotomize(txt, "kw2", sep = ";"))
di_case("kw_min2",        "dichotomize(txt, \"kw2\", sep = \";\", min = 2)", dichotomize(txt, "kw2", sep = ";", min = 2))
di_case("kw_length2",     "dichotomize(txt, \"kw2\", sep = \";\", length = 2)", dichotomize(txt, "kw2", sep = ";", length = 2))
di_case("kw_nosort",      "dichotomize(txt, \"kw2\", sep = \";\", sort = FALSE)", dichotomize(txt, "kw2", sep = ";", sort = FALSE))
di_case("kw_nas_null",    "dichotomize(txt, \"kw2\", sep = \";\", nas = NULL)", dichotomize(txt, "kw2", sep = ";", nas = NULL))
di_case("kw_nas_label",   "dichotomize(txt, \"kw2\", sep = \";\", nas = \"Missing\")", dichotomize(txt, "kw2", sep = ";", nas = "Missing"))
di_case("kw_values",      "dichotomize(txt, \"kw2\", sep = \";\", values = c(\"red\",\"blue\"))", dichotomize(txt, "kw2", sep = ";", values = c("red", "blue")))
di_case("kw_add",         "dichotomize(txt, \"kw2\", sep = \";\", add = TRUE)", dichotomize(txt, "kw2", sep = ";", add = TRUE))
di_case("kw_spaces",       "dichotomize(txt, \"keywords\", sep = \";\")", dichotomize(txt, "keywords", sep = ";"))
di_case("words_lower",    "dichotomize(txt, \"sentence\", sep = \" \", text = \"lower\", stopwords = c(\"the\",\"a\",\"and\"))",
        dichotomize(txt, "sentence", sep = " ", text = "lower", stopwords = c("the", "a", "and")))
di_case("two_seps",       "dichotomize(txt, c(\"kw2\",\"colour\"), sep = c(\";\",\"\"))", dichotomize(txt, c("kw2", "colour"), sep = c(";", "")))
di_case("regex_pool",     "dichotomize(txt, c(\"kw2\",\"sentence\"), sep = \"[;, ]+\", regex = TRUE, text = \"lower\")",
        dichotomize(txt, c("kw2", "sentence"), sep = "[;, ]+", regex = TRUE, text = "lower"))
di_case("sep_c",          "dichotomize(txt, c(\"grp\",\"colour\"), sep = \"C\")", dichotomize(txt, c("grp", "colour"), sep = "C"))
di_case("sep_c_nas",      "dichotomize(txt, c(\"grp\",\"colour\"), sep = \"C\", nas = \"empty\")", dichotomize(txt, c("grp", "colour"), sep = "C", nas = "empty"))
di_case("min_fraction",   "dichotomize(txt, \"kw2\", sep = \";\", min = .3)", dichotomize(txt, "kw2", sep = ";", min = .3))
di_case("values_nomin",   "dichotomize(txt, \"kw2\", sep = \";\", values = c(\"verde\",\"red\",\"zzz\"), min = 0, sort = FALSE)",
        dichotomize(txt, "kw2", sep = ";", values = c("verde", "red", "zzz"), min = 0, sort = FALSE))
di_case("regex_ws",       "dichotomize(txt, \"sentence\", sep = \"\\s+\", regex = TRUE, text = \"upper\", length = 3)",
        dichotomize(txt, "sentence", sep = "\\s+", regex = TRUE, text = "upper", length = 3))
# R keeps an unnamed column for empty tokens (e.g. "a;;b" or sep = "C" with a ""
# value); python drops them, so no reference is generated for those inputs.
# R's dichotomize(sep = "C") on c("keywords","colour") shows it: column "V3".


# ------------------------------------------------------ coocur / propCoin ----
counts <- data.frame(a = c(2, 0, 1, 3, 0, 1), b = c(1, 1, 0, 2, 2, 0), c = c(0, 2, 2, 1, 1, 1),
                     d = c(1, 0, 0, 0, 1, 0))
write_csv(counts, "counts", dir_in)
coocur_case <- function(id, call_txt, obj) {
  write_matrix(unclass(obj)[, ], paste0("coocur_", id))
  write_csv(data.frame(n = attr(obj, "n"), m = attr(obj, "m")), paste0("coocur_", id, "_attrs"), dir_out)
  add_case(paste0("coocur_", id), call_txt, paste0("coocur_", id))
}
coocur_case("default", "coocur(counts)", coocur(counts))
coocur_case("min3_sort", "coocur(counts, minimum = 3, sort = TRUE)", coocur(counts, minimum = 3, sort = TRUE))

prop_case <- function(id, call_txt, m) {
  write_matrix(m, paste0("prop_", id))
  add_case(paste0("prop_", id), call_txt, paste0("prop_", id))
}
for (mg in 0:2) prop_case(paste0("margin", mg), sprintf("propCoin(coin(galapagos), margin = %d)", mg),
                          propCoin(C_gal, margin = mg))
prop_case("margin0_dec3", "propCoin(coin(galapagos), decimals = 3)", propCoin(C_gal, decimals = 3))
prop_case("pairwise", "propCoin(coin(galapagos_na, pairwise = TRUE))", propCoin(C_pw))

# --------------------------------------------------------------- asNodes ----
nodes_case <- function(id, call_txt, x) {
  write_csv(x, paste0("nodes_", id), dir_out)
  add_case(paste0("nodes_", id), call_txt, paste0("nodes_", id))
}
nodes_case("default", "asNodes(coin(galapagos))", asNodes(C_gal))
nodes_case("both", "asNodes(coin(galapagos), TRUE, TRUE)", asNodes(C_gal, TRUE, TRUE))
nodes_case("pct_es", "asNodes(coin(galapagos), FALSE, TRUE, \"es\")", asNodes(C_gal, FALSE, TRUE, "es"))
nodes_case("freq_ca", "asNodes(coin(galapagos), TRUE, FALSE, \"ca\")", asNodes(C_gal, TRUE, FALSE, "ca"))
nodes_case("pairwise", "asNodes(coin(galapagos_na, pairwise = TRUE), TRUE, TRUE)", asNodes(C_pw, TRUE, TRUE))
nodes_case("edges", "asNodes(edgeList(coin(galapagos)))", asNodes(edgeList(C_gal)))

# ---------------------------------------------------------------- allNet ----
user_nodes <- data.frame(name = c("Geospiza fortis", "Certhidea olivacea", "Geospiza magnirostris", "Nonexistent"),
                         group = c("g1", "c", "g1", "x"), frequency = 999, stringsAsFactors = FALSE)
write_csv(user_nodes, "user_nodes", dir_in)
allnet_case <- function(id, call_txt, N) {
  nd <- N$nodes; rownames(nd) <- NULL
  write_csv(nd, paste0("allnet_", id, "_nodes"), dir_out)
  lk <- N$links; if (!is.null(lk)) { rownames(lk) <- NULL; lk[] <- lapply(lk, function(x) if (is.factor(x)) as.character(x) else x) }
  write_csv(lk, paste0("allnet_", id, "_links"), dir_out)
  keep <- c("nodeName", "nodeSize", "nodeLabel", "linkWidth", "linkColor", "scenarios", "language")
  o <- N$options[intersect(keep, names(N$options))]
  write_json(lapply(o, function(x) unname(x)), file.path(dir_out, paste0("allnet_", id, "_options.json")), auto_unbox = TRUE)
  add_case(paste0("allnet_", id), call_txt, paste0("allnet_", id))
}
allnet_case("default", "allNet(galapagos)", allNet(galapagos))
allnet_case("links", "allNet(galapagos, procedures = c(\"Jaccard\",\"Haberman\"), lwidth = \"Jaccard\", lcolor = \"Z\", criteria = \"Jaccard\", minL = .3, language = \"es\", frequency = TRUE, percentages = FALSE, layout = \"pc\")",
            allNet(galapagos, procedures = c("Jaccard", "Haberman"), lwidth = "Jaccard", lcolor = "Z",
                   criteria = "Jaccard", minL = .3, language = "es", frequency = TRUE, percentages = FALSE, layout = "pc"))
allnet_case("nodes", "allNet(galapagos, nodes = user_nodes, name = \"name\")", allNet(galapagos, nodes = user_nodes))
allnet_case("weight", "allNet(galapagos, weight = weights, support = 3, sort = TRUE)",
            allNet(galapagos, weight = weights, support = 3, sort = TRUE))
allnet_case("mca", "allNet(dice[,2:11], layout = \"mc\", procedures = \"Jaccard\", criteria = \"Jaccard\", minL = .1)",
            allNet(inc_dice, layout = "mc", procedures = "Jaccard", criteria = "Jaccard", minL = .1))
# NOTE: allNet(galapagos_na, pairwise = TRUE) fails in netCoin 2.1.49 (colSums() is NA with
# missing values and is used as a column index), so it is not used as reference.
allnet_case("bonferroni", "allNet(dice[,2:11], Bonferroni = TRUE, maxL = .2, directed = TRUE, sortL = \"Haberman\")",
            allNet(inc_dice, Bonferroni = TRUE, maxL = .2, directed = TRUE, sortL = "Haberman"))


# ------------------------------------------------------- evolving network ----
evo_nodes <- list(
  first  = data.frame(name = c("a", "b", "c"),      v = c(1, 2, 3),     g = c("x", "x", "y"), stringsAsFactors = FALSE),
  second = data.frame(name = c("a", "b", "c", "d"), v = c(1, 5, 3, 4),  g = c("x", "x", "y", "y"), stringsAsFactors = FALSE),
  third  = data.frame(name = c("b", "c", "d"),      v = c(5, 3, 4),     g = c("x", "z", "y"), stringsAsFactors = FALSE))
evo_links <- list(
  first  = data.frame(Source = c("a", "b"),           Target = c("b", "c"),           w = c(1, 2), stringsAsFactors = FALSE),
  second = data.frame(Source = c("a", "c"),           Target = c("b", "d"),           w = c(3, 4), stringsAsFactors = FALSE),
  third  = data.frame(Source = c("b", "c", "b"),      Target = c("c", "d", "d"),      w = c(5, 6, 7), stringsAsFactors = FALSE))
evo_nets <- lapply(names(evo_nodes), function(k) {
  write_csv(evo_nodes[[k]], paste0("evo_nodes_", k), dir_in)
  write_csv(evo_links[[k]], paste0("evo_links_", k), dir_in)
  netCoin(evo_nodes[[k]], evo_links[[k]], main = paste("Frame", k), repulsion = if (k == "third") 50 else 25)
})
names(evo_nets) <- names(evo_nodes)
E <- do.call(multigraphCreate, c(evo_nets, list(mode = "frame", frame = 1, speed = 30, loop = TRUE,
                                                lineplots = c("v", "degree", "bogus"))))
nd <- E$nodes; rownames(nd) <- NULL
write_csv(nd, "evo_result_nodes", dir_out)
lk <- E$links; rownames(lk) <- NULL
write_csv(lk, "evo_result_links", dir_out)
keep <- c("frames", "frame", "speed", "loop", "lineplotsKeys", "main", "note", "repulsion", "distance", "zoom", "nodeName")
write_json(lapply(E$options[intersect(keep, names(E$options))], function(x) unname(x)),
           file.path(dir_out, "evo_result_options.json"), auto_unbox = FALSE, na = "null")
add_case("evo_result", "multigraphCreate(first=, second=, third=, mode = \"frame\", frame = 1, speed = 30, loop = TRUE, lineplots = c(\"v\",\"degree\",\"bogus\"))",
         c("evo_result_nodes", "evo_result_links"))


# ------------------------------------------------- bar / time / pie charts ----
jsonw <- function(x, file) write_json(x, file.path(dir_out, file), auto_unbox = TRUE, na = "null", digits = NA, null = "null")
plain <- function(o) lapply(o, function(x) if (is.factor(x)) as.character(x) else unname(x))
write_table <- function(df, name) {
  df <- as.data.frame(df); rownames(df) <- NULL
  df[] <- lapply(df, function(x) if (is.factor(x)) as.character(x) else x)
  write_csv(df, name, dir_out)
}
bar_case <- function(id, call_txt, B) {
  write_table(B$nodes, paste0("bar_", id, "_nodes"))
  write_table(B$links, paste0("bar_", id, "_links"))
  jsonw(plain(B$options), paste0("bar_", id, "_options.json"))
  add_case(paste0("bar_", id), call_txt, paste0("bar_", id, c("_nodes", "_links")))
}
ess_f <- ess
bar_nodes_user <- data.frame(name = c("Geospiza fortis", "Certhidea olivacea", "Geospiza magnirostris", "Missing species"),
                             label = c("Fortis", "Olivacea", "Magnirostris", "?"), group = c("g", "c", "g", "x"),
                             stringsAsFactors = FALSE)
write_csv(bar_nodes_user, "bar_nodes_user", dir_in)
ess_w <- data.frame(lapply(ess[, 1:4], as.character), check.names = FALSE, stringsAsFactors = FALSE)
ess_w$cweight <- ess$cweight
write_csv(ess_w, "ess_cweight", dir_in)
bar_case("categorical", 'barCoin(ess, c("Gender","Age"), dichotomies = "_none")',
         barCoin(ess_f, c("Gender", "Age"), dichotomies = "_none"))
bar_case("full", 'barCoin(ess, 4 variables, dichotomies = c("Social participation","Political participation"), valueDicho = "Yes", expected = TRUE, confidence = TRUE, significance = TRUE, level = .9, percentages = TRUE, language = "es", sort = "incidences", select = "Gender:Male")',
         barCoin(ess_f, c("Gender", "Age", "Social participation", "Political participation"),
                 dichotomies = c("Social participation", "Political participation"), valueDicho = "Yes",
                 expected = TRUE, confidence = TRUE, significance = TRUE, level = .9, percentages = TRUE,
                 language = "es", sort = "incidences", select = "Gender:Male"))
bar_case("weight", 'barCoin(ess, c("Gender","Age","cweight"), dichotomies = "_none", weight = "cweight", criteria = "hyp", support = 5, maxL = .6, confidence = TRUE)',
         barCoin(ess_f, c("Gender", "Age", "cweight"), dichotomies = "_none", weight = "cweight", criteria = "hyp",
                 support = 5, maxL = .6, confidence = TRUE))
bar_case("dichotomies", "barCoin(galapagos)", barCoin(galapagos))
bar_case("nodes", 'barCoin(galapagos, nodes = bar_nodes_user, label = "label", color = "group", expected = TRUE, decreasing = FALSE, sort = "label", scalebar = TRUE, note = "A note", cex = 1.2, defaultColor = "firebrick")',
         barCoin(galapagos, nodes = bar_nodes_user, label = "label", color = "group", expected = TRUE,
                 decreasing = FALSE, sort = "label", scalebar = TRUE, note = "A note", cex = 1.2, defaultColor = "firebrick"))
bar_case("commonlabel", 'barCoin(ess, c("Gender","Age"), dichotomies = "_none", commonlabel = "Gender")',
         barCoin(ess_f, c("Gender", "Age"), dichotomies = "_none", commonlabel = "Gender"))

# timeline
data(sociologists, works, package = "netCoin")
soc <- as.data.frame(sociologists); wrk <- as.data.frame(works)
write_csv(soc, "soc", dir_in); write_csv(wrk, "wrk", dir_in)
time_case <- function(id, call_txt, T) {
  write_table(T$periods, paste0("time_", id, "_periods"))
  files <- paste0("time_", id, "_periods")
  if (!is.null(T$events)) { write_table(T$events, paste0("time_", id, "_events")); files <- c(files, paste0("time_", id, "_events")) }
  jsonw(plain(T$options), paste0("time_", id, "_options.json"))
  add_case(paste0("time_", id), call_txt, files)
}
time_case("periods", 'timeCoin(sociologists, start = "birth", end = "death", group = "school", text = "label", main = "Sociologists", info = "birth_country", language = "es", cex = 1.1)',
          timeCoin(soc, start = "birth", end = "death", group = "school", text = "label", main = "Sociologists",
                   info = "birth_country", language = "es", cex = 1.1))
time_case("events", 'timeCoin(sociologists, start = "birth", end = "death", events = works, eventNames = "works", eventPeriod = "name", eventTime = "date", eventColor = "label", eventShape = "name")',
          timeCoin(soc, start = "birth", end = "death", events = wrk, eventNames = "works", eventPeriod = "name",
                   eventTime = "date", eventColor = "label", eventShape = "name"))
posix <- data.frame(name = c("a", "b", "c"),
                    start = as.POSIXct(c("2020-01-01 10:00:00", "2020-01-02 12:30:00", "2020-01-03 00:00:00"), tz = "UTC"),
                    end = as.POSIXct(c("2020-01-05 10:00:00", "2020-01-04 08:15:00", "2020-01-06 18:00:00"), tz = "UTC"))
write_csv(data.frame(name = posix$name, start = format(posix$start, "%Y-%m-%d %H:%M:%S"), end = format(posix$end, "%Y-%m-%d %H:%M:%S")), "posix", dir_in)
time_case("posix", "timeCoin(posix)", timeCoin(posix))

# incTime / mobileEdges
inc <- incTime(soc, "name", "birth", "death")
write_csv(cbind(year = rownames(inc), inc), "inctime", dir_out)
add_case("inctime", 'incTime(sociologists, "name", "birth", "death")', "inctime")
for (d in c(0, 10)) {
  me <- mobileEdges(soc, "name", "birth", d)
  write_table(me, paste0("mobile_", d)); add_case(paste0("mobile_", d), sprintf('mobileEdges(sociologists, "name", "birth", %d)', d), paste0("mobile_", d))
}

# pie
pie_case <- function(id, call_txt, P) {
  o <- list(data = P$data, dim = as.integer(P$dim), colnames = P$colnames, rownames = P$rownames,
            dataw = P$dataw, nodes = P$nodes, links = P$links, options = plain(P$options))
  for (k in c("nodes", "links")) if (!is.null(o[[k]])) { o[[k]] <- as.data.frame(o[[k]]); rownames(o[[k]]) <- NULL; o[[k]][] <- lapply(o[[k]], function(x) if (is.factor(x)) as.character(x) else x) }
  write_json(o, file.path(dir_out, paste0("pie_", id, ".json")), auto_unbox = TRUE, na = "null", digits = NA, dataframe = "columns")
  add_case(paste0("pie_", id), call_txt, paste0("pie_", id, ".json"))
}
pie_nodes <- data.frame(name = c("Geospiza fortis", "Certhidea olivacea", "Unknown"), size = c(1, 2, 3), stringsAsFactors = FALSE)
write_csv(pie_nodes, "pie_nodes", dir_in)
pie_links <- edgeList(C_gal, c("Haberman", "Jaccard"), criteria = "Jaccard", min = .8)
write_csv(pie_links, "pie_links", dir_in)
pie_case("default", "pieCoin(coin(galapagos))", pieCoin(C_gal))
pie_case("options", 'pieCoin(coin(galapagos), expected = FALSE, abline = 3, main = "Pies", note = "n", language = "es", colors = c("red","#00F","lightblue","gray 50"), cex = 1.2, showLegend = FALSE, helpOn = TRUE, help = "hi")',
         pieCoin(C_gal, expected = FALSE, abline = 3, main = "Pies", note = "n", language = "es",
                 colors = c("red", "#00F", "lightblue", "gray 50"), cex = 1.2, showLegend = FALSE, helpOn = TRUE, help = "hi"))
pie_case("links_column", 'pieCoin(coin(galapagos), nodes = pie_nodes, links = pie_links, lcolor = "Haberman")',
         pieCoin(C_gal, nodes = pie_nodes, links = pie_links, lcolor = "Haberman"))
pie_case("links_color", 'pieCoin(coin(galapagos), links = pie_links, lcolor = "firebrick")',
         pieCoin(C_gal, links = pie_links, lcolor = "firebrick"))
pie_case("links_numeric", 'pieCoin(coin(galapagos), links = pie_links, lcolor = seq_len(nrow(pie_links)))',
         pieCoin(C_gal, links = pie_links, lcolor = seq_len(nrow(pie_links))))

write_json(manifest, file.path(root, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Reference written:", length(manifest$cases), "cases\n")
