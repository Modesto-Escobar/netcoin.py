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


# ------------------------------------------------------------ netCoin() ----
net_nodes <- data.frame(
  name  = letters[1:10],
  label = paste("Node", LETTERS[1:10]),
  group = c("g1", "g2", "g1", "g3", "g2", "g1", "g3", "g2", "g1", "g3"),
  kind  = c("x", "y", "x", "y", "x", "y", "x", "y", "x", "y"),
  size  = c(5, 3, 8, 1, 9, 4, 7, 2, 6, 10),
  value = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5),
  text  = paste0("t", 1:10),
  stringsAsFactors = FALSE)
net_links <- data.frame(
  Source = c("a", "a", "b", "c", "d", "e", "f", "g", "h", "i", "a", "c", "e", "g"),
  Target = c("b", "c", "c", "d", "e", "f", "g", "h", "i", "j", "d", "e", "g", "i"),
  weight = c(.9, .5, .2, .7, .4, .8, .1, .6, .3, .95, .55, .25, .45, .75),
  type   = rep(c("p", "q", "r"), length.out = 14),
  note   = paste0("l", 1:14),
  stringsAsFactors = FALSE)
net_tree <- data.frame(parent = c("a", "a", "b", "c"), child = c("b", "c", "d", "e"), stringsAsFactors = FALSE)
write_csv(net_nodes, "net_nodes", dir_in); write_csv(net_links, "net_links", dir_in); write_csv(net_tree, "net_tree", dir_in)

net_case <- function(id, call_txt, N) {
  nd <- as.data.frame(N$nodes); rownames(nd) <- NULL
  nd[] <- lapply(nd, function(x) if (is.factor(x)) as.character(x) else x)
  write_csv(nd, paste0("net_", id, "_nodes"), dir_out)
  if (!is.null(N$links)) { lk <- as.data.frame(N$links); rownames(lk) <- NULL; lk[] <- lapply(lk, function(x) if (is.factor(x)) as.character(x) else x); write_csv(lk, paste0("net_", id, "_links"), dir_out) }
  o <- list(options = plain(N$options), nlinks = if (is.null(N$links)) 0 else nrow(N$links))
  if (!is.null(N$tree)) o$tree <- lapply(as.list(N$tree), unname)
  if (!is.null(N$layouts)) o$layouts <- lapply(N$layouts, function(m) unname(split(unname(m), row(m))))
  jsonw(o, paste0("net_", id, ".json"))
  add_case(paste0("net_", id), call_txt, paste0("net_", id, c("_nodes", "_links")))
}

net_case("basic", 'netCoin(nodes, links, label = "label", size = "size", color = "group", shape = "kind", group = "group", lwidth = "weight", lcolor = "type", ltext = "note", ntext = "text", info = "text", main = "A network", note = "A note", scenarios = 10)',
         netCoin(net_nodes, net_links, label = "label", size = "size", color = "group", shape = "kind", group = "group",
                 lwidth = "weight", lcolor = "type", ltext = "note", ntext = "text", info = "text",
                 main = "A network", note = "A note", scenarios = 10))
net_case("single", 'netCoin(nodes, links, color = "firebrick", shape = "Square", lcolor = "#00f", size = 12, label = "", legend = "kind", border = "group", labelSize = "size")',
         netCoin(net_nodes, net_links, color = "firebrick", shape = "Square", lcolor = "#00f", size = 12, label = "",
                 legend = "kind", border = "group", labelSize = "size"))
net_case("numeric", 'netCoin(nodes, links, color = "value", nodeScaleLimits = c(0, 5, 14), lcolor = "weight", linkScaleLimits = c(0, 1), intensity = "weight", nodeBipolar = TRUE, linkBipolar = TRUE)',
         netCoin(net_nodes, net_links, color = "value", nodeScaleLimits = c(0, 5, 14), lcolor = "weight",
                 linkScaleLimits = c(0, 1), intensity = "weight", nodeBipolar = TRUE, linkBipolar = TRUE))
net_case("lists", 'netCoin(nodes, links, color = list(group = c(g1 = "red", g2 = "blue", g3 = "green")), shape = list(kind = c(x = "Circle", y = "star")), lcolor = list(type = c(p = "orange", q = "grey", r = "black")))',
         netCoin(net_nodes, net_links, color = list(group = c(g1 = "red", g2 = "blue", g3 = "green")),
                 shape = list(kind = c(x = "Circle", y = "star")),
                 lcolor = list(type = c(p = "orange", q = "grey", r = "black"))))
net_case("vectors", 'netCoin(nodes, links, color = <10 colours>, shape = 1:10, lcolor = 1:14)',
         netCoin(net_nodes, net_links, color = rep(c("red", "#0000FF", "light blue", "gray 50", "Forest Green"), 2),
                 shape = 1:10, lcolor = 1:14))
net_case("categories", 'netCoin(nodes, links, color = factor(group), shape = factor(kind), lcolor = factor(type))',
         netCoin(net_nodes, net_links, color = factor(net_nodes$group), shape = factor(net_nodes$kind),
                 lcolor = factor(net_links$type)))
net_case("filters", 'netCoin(nodes, links, nodeFilter = "size > 2", linkFilter = "weight > 0.3")',
         netCoin(net_nodes, net_links, nodeFilter = "size > 2", linkFilter = "weight > 0.3"))
net_case("degree", 'netCoin(nodes, links, degreeFilter = c(2, 4))', netCoin(net_nodes, net_links, degreeFilter = c(2, 4)))
net_case("degree_one", 'netCoin(nodes, links, degreeFilter = 3, nodeFilter = "kind == \'x\'")',
         netCoin(net_nodes, net_links, degreeFilter = 3, nodeFilter = "kind == 'x'"))
net_case("options", 'netCoin(nodes, links, controls = c(1,2,3), showLegend = FALSE, mode = "network", repulsion = 40, distance = 15, zoom = 2, limits = c(0,0,10,10), background = "#ABC", defaultColor = "firebrick", helpOn = TRUE, help = "Help text", showCoordinates = TRUE, showArrows = TRUE, showAxes = TRUE, axesLabels = c("X","Y"), nodesLabel = "N", linksLabel = "L", sort = "size", decreasing = TRUE, frequencies = TRUE, statistics = TRUE, roundedItems = TRUE, cex = 1.3, groupText = TRUE, group = "group")',
         netCoin(net_nodes, net_links, controls = c(1, 2, 3), showLegend = FALSE, mode = "network", repulsion = 40,
                 distance = 15, zoom = 2, limits = c(0, 0, 10, 10), background = "#ABC", defaultColor = "firebrick",
                 helpOn = TRUE, help = "Help text", showCoordinates = TRUE, showArrows = TRUE, showAxes = TRUE,
                 axesLabels = c("X", "Y"), nodesLabel = "N", linksLabel = "L", sort = "size", decreasing = TRUE,
                 frequencies = TRUE, statistics = TRUE, roundedItems = TRUE, cex = 1.3, groupText = TRUE, group = "group"))
net_case("controls", 'netCoin(nodes, links, controls = c(1,5,8), showLegend = FALSE)', netCoin(net_nodes, net_links, controls = c(1, 5, 8), showLegend = FALSE))
net_case("controls0", 'netCoin(nodes, links, controls = 0)', netCoin(net_nodes, net_links, controls = 0))
net_case("spanish", 'netCoin(nodes, links, name = "name", language = "es")', netCoin(net_nodes, net_links, language = "es", name = "name"))
net_case("tree", 'netCoin(nodes, links, tree = tree)', netCoin(net_nodes, net_links, tree = net_tree))
net_case("linksonly", 'netCoin(links = links)', netCoin(links = net_links))
net_case("nolinks", 'netCoin(nodes, size = "size", name = "label")', netCoin(net_nodes, size = "size", name = "label"))
net_case("images", 'netCoin(nodes, links, image = c("label","text"), imageNames = c("group","kind"))',
         netCoin(net_nodes, net_links, image = c("label", "text"), imageNames = c("group", "kind")))
net_case("sparse", 'netCoin(nodes[c(1,2,4)], links, name = "name")', netCoin(net_nodes[c(1, 2, 4), c("name", "size")], net_links))

# layouts: given matrices and deterministic names
lay <- cbind(x = seq(0, 9), y = rep(c(0, 1), 5))
net_case("layout_matrix", 'netCoin(nodes, links, layout = cbind(x=0:9, y=rep(0:1,5)))', netCoin(net_nodes, net_links, layout = lay))
lay_nocols <- unname(lay)
net_case("layout_nocols", 'netCoin(nodes, links, layout = unname(cbind(...)))', netCoin(net_nodes, net_links, layout = lay_nocols))
net_case("layout_planes", 'netCoin(nodes, links, layout = list(first = cbind(x,y), second = cbind(u,v)))',
         netCoin(net_nodes, net_links, layout = list(first = lay, second = cbind(u = lay[, 2], v = lay[, 1]))))
for (lname in c("circle", "grid", "star", "reingold_tilford", "kamada_kawai", "sugiyama")) {
  net_case(paste0("layout_", lname), sprintf('netCoin(nodes, links, layout = "%s")', lname), netCoin(net_nodes, net_links, layout = lname))
}
net_case("layout_two_names", 'netCoin(nodes, links, layout = c("circle", "grid"))', netCoin(net_nodes, net_links, layout = c("circle", "grid")))
net_case("layout_mixed", 'netCoin(nodes, links, layout = list(mine = cbind(x,y), "circle"))', netCoin(net_nodes, net_links, layout = list(mine = lay, "circle")))

# communities and centralities
for (cm in c("walktrap", "fast_greedy", "leading_eigen", "edge_betweenness")) {
  net_case(paste0("community_", cm), sprintf('netCoin(nodes, links, community = "%s")', cm), netCoin(net_nodes, net_links, community = cm))
}
net_case("community_color", 'netCoin(nodes, links, community = "walktrap", color = "community")', netCoin(net_nodes, net_links, community = "walktrap", color = "group"))
net_case("centrality", 'netCoin(nodes, links, centrality = c("wdegree","closeness","betweenness","eigen"))',
         netCoin(net_nodes, net_links, centrality = c("wdegree", "closeness", "betweenness", "eigen")))

# conversions and exports
G <- toIgraph(netCoin(net_nodes, net_links, label = "label", size = "size", color = "group", lwidth = "weight", lcolor = "type", showArrows = TRUE))
vtab <- igraph::as_data_frame(G, what = "vertices"); rownames(vtab) <- NULL
etab <- igraph::as_data_frame(G, what = "edges"); rownames(etab) <- NULL
write_table(vtab, "igraph_vertices"); write_table(etab, "igraph_edges")
write_json(list(directed = igraph::is_directed(G)), file.path(dir_out, "igraph_meta.json"), auto_unbox = TRUE)
add_case("igraph_to", 'toIgraph(netCoin(nodes, links, label, size, color, lwidth, lcolor, showArrows = TRUE))', c("igraph_vertices", "igraph_edges"))
cent <- netCoin:::calCentr(netCoin(net_nodes, net_links), "all")$nodes
write_table(cent, "centrality_all"); add_case("centrality_all", 'calCentr(netCoin(nodes, links), "all")', "centrality_all")
write_table(netCoin:::calCentr(netCoin(net_nodes, net_links, showArrows = TRUE), c("degree", "closeness"))$nodes, "centrality_directed")
add_case("centrality_directed", 'calCentr(netCoin(nodes, links, showArrows = TRUE), c("degree","closeness"))', "centrality_directed")
pj <- netCoin(net_nodes, net_links, main = "Pajek test", showArrows = TRUE)
pjf <- file.path(dir_out, "pajek_basic.net"); savePajek(pj, pjf, arcs = "weight", edges = c("weight"))
pjp <- file.path(dir_out, "pajek_project.paj"); savePajek(pj, file.path(dir_out, "pajek_project.net"), arcs = "weight", partitions = c("group", "kind"), vectors = "size")
file.rename(file.path(dir_out, "pajek_basic.net"), file.path(dir_out, "pajek_basic.txt"))
file.rename(pjp, file.path(dir_out, "pajek_project.txt"))
add_case("pajek", "savePajek(...)", c("pajek_basic", "pajek_project"))

# circle and grid layouts
write_matrix_plain <- function(m, name) { df <- as.data.frame(m); df <- cbind(.row = rownames(m), df); write_csv(df, name, dir_out) }
write_matrix_plain(layoutCircle(net_nodes, name = "name"), "layoutcircle_all")
write_matrix_plain(layoutCircle(net_nodes, name = "name", deg = 90), "layoutcircle_deg")
write_matrix_plain(layoutCircle(net_nodes, nodes = c("c", "a", "e", "zz"), name = "name"), "layoutcircle_subset")
write_matrix_plain(layoutCircle(net_nodes, nodes = "b", name = "name"), "layoutcircle_one")
gstrings <- c(simple = "a,b.c.d,e,f", gaps = "a;b,c.d,e.f", bottom = "a,b,.c.d,e", one = "a,b,c,d", byrow = "a,b,c.d,e.f,g")
write_json(as.list(gstrings), file.path(dir_out, "layoutgrid_strings.json"), auto_unbox = TRUE)
for (g in names(gstrings)) write_matrix_plain(layoutGrid(net_nodes, gstrings[[g]], name = "name", byrow = (g == "byrow")), paste0("layoutgrid_", g))

# ---------------------------------------------------- surCoin / cobCoin ----
# The survey is `ess` (one text column per variable, missing values included) with the
# weights; `ess_metric` adds two quantitative variables made up for the occasion.
set.seed(1)
ess_m <- ess_w
ess_m$income <- round(rnorm(nrow(ess_m), 50, 10) + as.integer(factor(ess_m$Age)), 1)
ess_m$edu    <- round(runif(nrow(ess_m)) * 10 + 2 * as.integer(factor(ess_m$Gender)), 1)
ess_m$income[c(3, 7, 400)] <- NA
write_csv(ess_m, "ess_metric", dir_in)
ess_cc <- ess_w[complete.cases(ess_w), ]
write_csv(ess_cc, "ess_complete", dir_in)
sur_nodes <- data.frame(name = c("Gender:Male", "Age:15-29", "Gender:Female", "Age:30-39", "Social participation:Yes", "Zzz:none"),
                        grp = c("a", "b", "a", "b", "c", "d"), stringsAsFactors = FALSE)
write_csv(sur_nodes, "sur_nodes", dir_in)
sur_nodes_metric <- data.frame(name = c("income", "Gender:Male", "Social participation:Yes", "Gender:Female", "edu"),
                          grp = c("m", "a", "c", "a", "m"), stringsAsFactors = FALSE)
write_csv(sur_nodes_metric, "sur_nodes_metric", dir_in)

vars3 <- c("Gender", "Age", "Social participation")
sur_case <- function(id, call_txt, expr) {
  tryCatch(net_case(paste0("sur_", id), call_txt, expr),
           error = function(e) message("surCoin case ", id, " FAILED: ", conditionMessage(e)))
}
sur_case("default", 'surCoin(ess, c("Gender","Age","Social participation"))', surCoin(ess_w, vars3))
sur_case("one_dicho", 'surCoin(ess, "Gender", dichotomies = "Social participation", valueDicho = "Yes", maxL = 1)',
         surCoin(ess_w, "Gender", dichotomies = "Social participation", valueDicho = "Yes", maxL = 1))
sur_case("dichos", 'surCoin(ess, c("Gender","Age"), dichotomies = c("Social participation","Political participation"), valueDicho = "Yes", maxL = 1)',
         surCoin(ess_w, c("Gender", "Age"), dichotomies = c("Social participation", "Political participation"),
                 valueDicho = "Yes", maxL = 1))
sur_case("dicho_inside", 'surCoin(ess, c("Social participation","Gender","Age"), dichotomies = "Social participation", valueDicho = "Yes", maxL = 1)',
         surCoin(ess_w, c("Social participation", "Gender", "Age"), dichotomies = "Social participation",
                 valueDicho = "Yes", maxL = 1))
sur_case("links", 'surCoin(ess, vars3, frequency = TRUE, percentages = FALSE, procedures = c("Haberman","Jaccard"), lwidth = "Jaccard", lcolor = "Haberman", criteria = "Jaccard", minL = .2, language = "es", sort = TRUE, decreasing = FALSE, sortL = "Jaccard", main = "Survey")',
         surCoin(ess_w, vars3, frequency = TRUE, percentages = FALSE, procedures = c("Haberman", "Jaccard"),
                 lwidth = "Jaccard", lcolor = "Haberman", criteria = "Jaccard", minL = .2, language = "es",
                 sort = TRUE, decreasing = FALSE, sortL = "Jaccard", main = "Survey"))
sur_case("minimum", 'surCoin(ess, vars3, minimum = 150, maxL = 1)', surCoin(ess_w, vars3, minimum = 150, maxL = 1))
sur_case("weight", 'surCoin(ess, c("Gender","Age","cweight"), weight = "cweight", maxL = 1)',
         surCoin(ess_w, c("Gender", "Age", "cweight"), weight = "cweight", maxL = 1))
sur_case("weight_vector", 'surCoin(ess[complete.cases(ess), ], vars3, weight = cweight, maxL = 1)',
         surCoin(ess_cc, vars3, weight = ess_cc$cweight, maxL = 1))
sur_case("subsample", 'surCoin(ess, vars3, subsample = TRUE, maxL = 1)', surCoin(ess_w, vars3, subsample = TRUE, maxL = 1))
sur_case("pairwise", 'surCoin(ess, vars3, pairwise = TRUE, maxL = 1)', surCoin(ess_w, vars3, pairwise = TRUE, maxL = 1))
sur_case("nodes", 'surCoin(ess, vars3, nodes = sur_nodes, color = "grp", maxL = 1)',
         surCoin(ess_w, vars3, nodes = sur_nodes, color = "grp", maxL = 1))
sur_case("layout_mc", 'surCoin(ess, vars3, layout = "mc", maxL = 1)', surCoin(ess_w, vars3, layout = "mc", maxL = 1))
sur_case("layout_pc", 'surCoin(ess, vars3, layout = "pc")', surCoin(ess_w, vars3, layout = "pc"))
sur_case("layout_two", 'surCoin(ess, vars3, layout = c("mc","pc"))', surCoin(ess_w, vars3, layout = c("mc", "pc")))
sur_case("nodes_mc", 'surCoin(ess, vars3, nodes = sur_nodes, layout = "mc", maxL = 1)',
         surCoin(ess_w, vars3, nodes = sur_nodes, layout = "mc", maxL = 1))
sur_lay <- cbind(x = c(0, 1, 2, 3, 4, 5), y = c(5, 4, 3, 2, 1, 0))
sur_case("layout_matrix", 'surCoin(ess, vars3, nodes = sur_nodes, layout = cbind(x = 0:5, y = 5:0), maxL = 1)',
         surCoin(ess_w, vars3, nodes = sur_nodes, layout = sur_lay, maxL = 1))
sur_case("commonlabel", 'surCoin(ess, c("Gender","Age"), commonlabel = "Gender", maxL = 1)',
         surCoin(ess_w, c("Gender", "Age"), commonlabel = "Gender", maxL = 1))
sur_case("exogenous", 'surCoin(ess, vars3, exogenous = c("Gender","Age"), maxL = 1)',
         surCoin(ess_w, vars3, exogenous = c("Gender", "Age"), maxL = 1))
sur_case("exogenous_filter", 'surCoin(ess, vars3, exogenous = "Gender", linkFilter = "Haberman > 0", maxL = 1)',
         surCoin(ess_w, vars3, exogenous = "Gender", linkFilter = "Haberman > 0", maxL = 1))
sur_case("spanish", 'surCoin(ess, vars3, language = "es", lwidth = "Haberman", ltext = "Z")',
         surCoin(ess_w, vars3, language = "es", lwidth = "Haberman", ltext = "Z"))
sur_case("metric", 'surCoin(ess_metric, c("Gender","Social participation"), metric = c("income","edu"), maxL = 1)',
         surCoin(ess_m, c("Gender", "Social participation"), metric = c("income", "edu"), maxL = 1))
sur_case("metric_pearson", 'surCoin(ess_metric, vars3, metric = "income", procedures = c("Pearson","Haberman"), criteria = "Pearson", minL = .05, maxL = Inf)',
         surCoin(ess_m, vars3, metric = "income", procedures = c("Pearson", "Haberman"), criteria = "Pearson",
                 minL = .05, maxL = Inf))
sur_case("metric_pairwise", 'surCoin(ess_metric, c("Gender","Age"), metric = c("income","edu"), pairwise = TRUE, maxL = 1)',
         surCoin(ess_m, c("Gender", "Age"), metric = c("income", "edu"), pairwise = TRUE, maxL = 1))
sur_case("metric_weight", 'surCoin(ess_metric, c("Gender","Age","cweight"), metric = "edu", weight = "cweight", maxL = 1)',
         surCoin(ess_m, c("Gender", "Age", "cweight"), metric = "edu", weight = "cweight", maxL = 1))
sur_case("metric_nodes", 'surCoin(ess_metric, vars3, metric = c("income","edu"), nodes = sur_nodes_metric, color = "grp", maxL = 1)',
         surCoin(ess_m, vars3, metric = c("income", "edu"), nodes = sur_nodes_metric, color = "grp", maxL = 1))
sur_case("metric_nopct", 'surCoin(ess_metric, vars3, metric = "edu", percentages = FALSE, frequency = TRUE, maxL = 1)',
         surCoin(ess_m, vars3, metric = "edu", percentages = FALSE, frequency = TRUE, maxL = 1))
sur_case("metric_exogenous", 'surCoin(ess_metric, vars3, metric = "income", exogenous = c("Gender","income"), maxL = 1)',
         surCoin(ess_m, vars3, metric = "income", exogenous = c("Gender", "income"), maxL = 1))
coin_case("sur", 'surCoin(ess, c("Gender","Age"), coin = TRUE)', surCoin(ess_w, c("Gender", "Age"), coin = TRUE))
coin_case("sur_weight", 'surCoin(ess, c("Gender","Age","cweight"), weight = "cweight", coin = TRUE)',
          surCoin(ess_w, c("Gender", "Age", "cweight"), weight = "cweight", coin = TRUE))

vars4 <- c("Gender", "Age", "Social participation", "Political participation")
sur_case("cob_default", 'cobCoin(ess, c("Gender","Age","Social participation","Political participation"))', cobCoin(ess_w, vars4))
sur_case("cob_args", 'cobCoin(ess, vars4, degree = 45, significance = .1, language = "es", main = "CoWeb")',
         cobCoin(ess_w, vars4, degree = 45, significance = .1, language = "es", main = "CoWeb"))

# --------------------------------------------------- netCorr / d_netCorr ----
suppressMessages(library(GPArotation))       # layoutFact(..., "of") needs it
mt <- data.frame(mtcars[, c("mpg", "disp", "hp", "drat", "wt", "qsec")], row.names = NULL)
set.seed(2)
mt_w <- cbind(mt, w = round(runif(nrow(mt), 0.3, 2), 2))
mt_na <- mt
mt_na[cbind(c(1, 4, 9, 15, 22, 30, 31, 5), c(1, 2, 3, 3, 5, 6, 4, 4))] <- NA
small <- data.frame(a = round(rnorm(7), 3), b = round(rnorm(7), 3), c = round(rnorm(7), 3), d = round(rnorm(7), 3))
small_ties <- data.frame(a = c(1, 2, 2, 3, 5, 5, 6, 8), b = c(2, 1, 4, 3, 3, 6, 8, 7),
                         c = c(5, 3, 2, 8, 1, 4, 4, 9), d = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5))
small_na <- small; small_na[c(2, 5), 1] <- NA; small_na[3, 3] <- NA; small_na[7, 4] <- NA
write_csv(mt, "corr_mt", dir_in); write_csv(mt_w, "corr_mt_w", dir_in); write_csv(mt_na, "corr_mt_na", dir_in)
write_csv(small, "corr_small", dir_in); write_csv(small_ties, "corr_small_ties", dir_in)
write_csv(small_na, "corr_small_na", dir_in)
nc_nodes <- data.frame(name = c("qsec", "mpg", "wt", "extra"), grp = c("x", "y", "x", "z"), stringsAsFactors = FALSE)
write_csv(nc_nodes, "corr_nodes", dir_in)
nc_nodes_es <- nc_nodes; names(nc_nodes_es)[1] <- "nombre"
write_csv(nc_nodes_es, "corr_nodes_es", dir_in)
# variables in blocks, so that the components split as the threshold rises
set.seed(5)
lat <- matrix(rnorm(60 * 3), 60, 3)
blocks <- data.frame(
  a1 = lat[, 1] + rnorm(60, sd = .5), a2 = lat[, 1] + rnorm(60, sd = .8), a3 = lat[, 1] + rnorm(60, sd = 1.2),
  b1 = lat[, 2] + rnorm(60, sd = .5), b2 = lat[, 2] + rnorm(60, sd = .9), b3 = lat[, 2] + .5 * lat[, 1] + rnorm(60, sd = 1),
  c1 = lat[, 3] + rnorm(60, sd = .6), c2 = lat[, 3] + .3 * lat[, 2] + rnorm(60, sd = 1.1))
blocks <- round(blocks, 3)
write_csv(blocks, "corr_blocks", dir_in)

nc_case <- function(id, call_txt, expr) {
  tryCatch(net_case(paste0("nc_", id), call_txt, expr),
           error = function(e) message("netCorr case ", id, " FAILED: ", conditionMessage(e)))
}
nc_case("default", "netCorr(mt)", netCorr(mt))
for (mth in c("spearman", "kendall", "partial", "dcor"))
  nc_case(mth, sprintf('netCorr(mt, method = "%s")', mth), netCorr(mt, method = mth))
nc_case("multi", 'netCorr(mt, method = c("pearson","spearman","kendall"), criteria = "p(kendall)")',
        netCorr(mt, method = c("pearson", "spearman", "kendall"), criteria = "p(kendall)"))
nc_case("small", 'netCorr(small, method = c("spearman","kendall"), criteria = "p")',
        netCorr(small, method = c("spearman", "kendall"), criteria = "p", maxL = 1))
nc_case("small_one", 'netCorr(small, method = c("spearman","kendall"), linkBipolar = FALSE, maxL = 1)',
        netCorr(small, method = c("spearman", "kendall"), linkBipolar = FALSE, maxL = 1))
nc_case("small_ties", 'netCorr(small_ties, method = c("spearman","kendall"), maxL = 1)',
        netCorr(small_ties, method = c("spearman", "kendall"), maxL = 1))
nc_case("small_na", 'netCorr(small_na, pairwise = TRUE, method = c("spearman","kendall","pearson"), maxL = 1)',
        netCorr(small_na, pairwise = TRUE, method = c("spearman", "kendall", "pearson"), maxL = 1))
nc_case("weight", 'netCorr(mt_w, weight = "w", method = c("pearson","spearman","partial"))',
        netCorr(mt_w, weight = "w", method = c("pearson", "spearman", "partial")))
nc_case("weight_vector", 'netCorr(mt, weight = mt_w$w, method = "spearman")', netCorr(mt, weight = mt_w$w, method = "spearman"))
nc_case("na_listwise", "netCorr(mt_na)", netCorr(mt_na))
nc_case("na_pairwise", 'netCorr(mt_na, pairwise = TRUE, method = c("pearson","spearman","kendall","dcor"))',
        netCorr(mt_na, pairwise = TRUE, method = c("pearson", "spearman", "kendall", "dcor")))
mt_naw <- cbind(mt_na, w = mt_w$w)
nc_case("na_pairwise_weight", 'netCorr(mt_na + w, weight = "w", pairwise = TRUE, method = c("pearson","spearman"))',
        netCorr(mt_naw, weight = "w", pairwise = TRUE, method = c("pearson", "spearman")))
nc_case("na_partial", 'netCorr(mt_na, pairwise = TRUE, method = "partial")', suppressWarnings(netCorr(mt_na, pairwise = TRUE, method = "partial")))
nc_case("value", 'netCorr(mt, criteria = "value", minL = .5)', netCorr(mt, criteria = "value", minL = .5))
nc_case("value_one", 'netCorr(mt, criteria = "value", minL = .5, maxL = .9, linkBipolar = FALSE)',
        netCorr(mt, criteria = "value", minL = .5, maxL = .9, linkBipolar = FALSE))
nc_case("one_sided", 'netCorr(mt, linkBipolar = FALSE)', netCorr(mt, linkBipolar = FALSE))
nc_case("bonferroni", 'netCorr(mt, Bonferroni = TRUE, maxL = .05)', netCorr(mt, Bonferroni = TRUE, maxL = .05))
nc_case("sort", 'netCorr(mt_na, pairwise = TRUE, sort = TRUE, decreasing = FALSE, sortL = "value", decreasingL = FALSE)',
        netCorr(mt_na, pairwise = TRUE, sort = TRUE, decreasing = FALSE, sortL = "value", decreasingL = FALSE))
nc_case("sort_decreasing", 'netCorr(mt_na, pairwise = TRUE, sort = TRUE, sortL = "z")', netCorr(mt_na, pairwise = TRUE, sort = TRUE, sortL = "z"))
nc_case("nomeans", "netCorr(mt, means = FALSE, frequency = TRUE)", netCorr(mt, means = FALSE, frequency = TRUE))
nc_case("nodes", 'netCorr(mt, nodes = nc_nodes_es, language = "es", color = "grp")', netCorr(mt, nodes = nc_nodes_es, language = "es", color = "grp"))
nc_case("layout_pc", 'netCorr(mt, layout = "pc")', netCorr(mt, layout = "pc"))
nc_case("layout_pc_nodes", 'netCorr(mt, nodes = nc_nodes, layout = "pc")', netCorr(mt, nodes = nc_nodes, layout = "pc"))
nc_case("options", 'netCorr(mt, lwidth = "z", lcolor = "p", size = "mean", color = "std", main = "Cars")',
        netCorr(mt, lwidth = "z", lcolor = "p", size = "mean", color = "std", main = "Cars"))

# distance correlation and factor layouts
for (f in c("pc", "nf", "vf", "of")) {
  L <- netCoin:::layoutFact(blocks, f)
  write_matrix_plain(L, paste0("layoutfact_", f))
}
Lall <- netCoin:::layoutFact(blocks, c("pc", "nf", "vf", "of"))
write_json(list(names = names(Lall)), file.path(dir_out, "layoutfact_names.json"), auto_unbox = TRUE)

dn_seq <- seq(.2, 1, .05)
nc_case("d_default", "d_netCorr(blocks, sequence = seq(.2, 1, .05))", d_netCorr(blocks, sequence = dn_seq))
nc_case("d_pc", 'd_netCorr(blocks, factorial = "pc", sequence = seq(.2, 1, .05), method = "spearman")',
        d_netCorr(blocks, factorial = "pc", sequence = dn_seq, method = "spearman"))
nc_case("d_fact", 'd_netCorr(blocks, factorial = c("pc","nf","vf"), sequence = seq(.2, 1, .05), weight = NULL)',
        d_netCorr(blocks, factorial = c("pc", "nf", "vf"), sequence = dn_seq))
nc_case("d_all", "d_netCorr(blocks, components = FALSE, sequence = seq(.3, .9, .1), means = FALSE)",
        d_netCorr(blocks, components = FALSE, sequence = seq(.3, .9, .1), means = FALSE))
nc_case("d_back", "d_netCorr(blocks, backcomponents = TRUE, sequence = seq(.2, 1, .05), textFilter = c(.5, .1))",
        d_netCorr(blocks, backcomponents = TRUE, sequence = dn_seq, textFilter = c(.5, .1)))

# ------------------------------------------------------------- logCoin ----
set.seed(7)
nl <- 320
lg <- data.frame(
  A = sample(c("a1", "a2", "a3"), nl, TRUE, c(.5, .3, .2)),
  B = sample(c("b1", "b2"), nl, TRUE),
  C = sample(c("c1", "c2"), nl, TRUE, c(.7, .3)),
  D = sample(c("d1", "d2", "d3"), nl, TRUE),
  E = sample(c("e1", "e2"), nl, TRUE),
  N = sample(1:3, nl, TRUE),
  W = round(runif(nl, .5, 2), 2),
  stringsAsFactors = FALSE)
lg$B[lg$A == "a1" & runif(nl) < .5] <- "b1"                        # some association
lg$C[lg$B == "b2" & runif(nl) < .4] <- "c2"
lg$E[lg$D == "d3" & runif(nl) < .6] <- "e2"
lg$A[c(3, 50)] <- NA; lg$E[c(10, 11, 200)] <- NA
names(lg)[names(lg) == "D"] <- "D var"
write_csv(lg, "log_cat5", dir_in)

lc_case <- function(id, call_txt, expr) {
  tryCatch(net_case(paste0("lc_", id), call_txt, expr),
           error = function(e) message("logCoin case ", id, " FAILED: ", conditionMessage(e)))
}
vars3l <- c("Gender", "Age", "Social participation")
lc_case("default", 'logCoin(ess, c("Gender","Age","Social participation"))', logCoin(ess_w, vars3l))
lc_case("twotails", 'logCoin(ess, c("Gender","Age","Political participation"), frequency = TRUE, percentage = TRUE, twotails = TRUE, pmax = .1)',
        logCoin(ess_w, c("Gender", "Age", "Political participation"), frequency = TRUE, percentage = TRUE, twotails = TRUE, pmax = .1))
lc_case("order3", 'logCoin(ess, c("Age","Social participation","Political participation","Gender"), exogenous = "Gender", order = 3, pmax = .2)',
        logCoin(ess_w, c("Age", "Social participation", "Political participation", "Gender"), exogenous = "Gender", order = 3, pmax = .2))
lc_case("order3_freq", 'logCoin(ess, 4 variables, exogenous = "Gender", order = 3, pmax = .2, frequency = TRUE, percentage = TRUE)',
        logCoin(ess_w, c("Age", "Social participation", "Political participation", "Gender"), exogenous = "Gender", order = 3, pmax = .2,
                frequency = TRUE, percentage = TRUE))
lc_case("nofirst", 'logCoin(ess, vars3, noFirstCat = "Gender")', logCoin(ess_w, vars3l, noFirstCat = "Gender"))
lc_case("nofirst_two", 'logCoin(ess, 4 variables, noFirstCat = c("Gender","Age"), twotails = TRUE)',
        logCoin(ess_w, c("Gender", "Age", "Social participation", "Political participation"), noFirstCat = c("Gender", "Age"), twotails = TRUE))
lc_case("weight", 'logCoin(ess, c("Gender","Age","Social participation"), weight = "cweight")',
        logCoin(ess_w, vars3l, weight = "cweight"))
lc_case("spanish", 'logCoin(ess, vars3, language = "es", frequency = TRUE)', logCoin(ess_w, vars3l, language = "es", frequency = TRUE))
lc_case("four", 'logCoin(ess, 4 variables, pmax = .1)',
        logCoin(ess_w, c("Gender", "Age", "Social participation", "Political participation"), pmax = .1))
lc_case("five", 'logCoin(lg, c("A","B","C","D var","E"), pmax = .1)', logCoin(lg, c("A", "B", "C", "D var", "E"), pmax = .1))
lc_case("five_two", 'logCoin(lg, c("A","B","C","D var","E"), twotails = TRUE, pmax = .2, frequency = TRUE)',
        logCoin(lg, c("A", "B", "C", "D var", "E"), twotails = TRUE, pmax = .2, frequency = TRUE))
lc_case("numeric", 'logCoin(lg, c("N","B","C"), weight = "W", pmax = .3)', logCoin(lg, c("N", "B", "C"), weight = "W", pmax = .3))
lc_case("exo_two", 'logCoin(lg, c("A","B","C","E"), exogenous = c("A","E"), order = 3, pmax = .3, noFirstCat = "B")',
        logCoin(lg, c("A", "B", "C", "E"), exogenous = c("A", "E"), order = 3, pmax = .3, noFirstCat = "B"))
lc_case("options", 'logCoin(ess, vars3, color = "cat", lwidth = "z value", main = "Log-linear", degreeFilter = 0)',
        logCoin(ess_w, vars3l, color = "cat", lwidth = "z value", main = "Log-linear", degreeFilter = 0))

write_json(manifest, file.path(root, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Reference written:", length(manifest$cases), "cases\n")
