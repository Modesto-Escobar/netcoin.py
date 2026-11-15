# -*- coding: utf-8 -*-
"""Guided tour of what Phase 2 added. Run it from any folder:

    python python_examples/review_phase2.py          # prints each step
    python python_examples/review_phase2.py --open   # also opens the networks in the browser

The same steps in R are in review_phase2.R, to compare the numbers.
"""
import os
import sys

import pandas as pd

# use the netcoin of this repository, even if another copy is installed in the environment
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from netcoin import (all_net, as_nodes, coin, coocur, dichotomize, edge_list,
                     expected_list, prop_coin, sim)

here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "tmp")                      # ignored by git
show = "--open" in sys.argv
pd.set_option("display.width", 160, "display.max_columns", 12)


def step(title):
    print("\n" + "=" * 8, title, "=" * 8)


# Galapagos: 17 islands (rows) x 13 finch species (columns), 1 = the species lives there
galapagos = pd.read_csv(os.path.join(here, "..", "tests", "reference", "inputs", "galapagos.csv"))

# ------------------------------------------------------------------ 1 ----
step("1. dichotomize: text -> indicator columns (R-style)")
books = pd.DataFrame({"keywords": ["red;blue", "blue;green", "red", "Green;red;blue", None, "red;red"]})
print(dichotomize(books, "keywords", sep=";"))        # note the 'keywords:None' column, as in R
print(dichotomize(books, "keywords", sep=";", nas=None, length=2))

words = pd.DataFrame({"s": ["The cat and the dog", "a Dog, a cat!", "the bird"]})
print(dichotomize(words, "s", sep=" ", text="lower", stopwords=["the", "a", "and"], nas=None))

# ------------------------------------------------------------------ 2 ----
step("2. coin: a Coin object carries the matrix AND the number of cases")
c = coin(galapagos)
print(c.n, "cases,", len(c), "species")
print(c)                                              # same print as R: n= 17 + lower triangle
print(c.frequencies.head(3))
print("sorted, minimum 10:", coin(galapagos, minimum=10, sort=True).labels)

# ------------------------------------------------------------------ 3 ----
step("3. sim: similarities, with R's names and abbreviations")
print(sim(c, "Jaccard").iloc[:3, :3])
print(sim(c, "J").equals(sim(c, "Jaccard")), "<- 'J' and 'Jaccard' are the same")
many = sim(c, ["Jaccard", "Haberman", "confidence"])  # several -> dict of matrices
print(list(many))

# ------------------------------------------------------------------ 4 ----
step("4. edge_list: ordered as R does (Source = earlier variable)")
edges = edge_list(c, ["Jaccard", "Haberman"], criteria="Jaccard", min=0.7, sort="Jaccard")
print(edges.head(6))
print("p-value criterion defaults to max=0.5:\n", edge_list(c).head(3))
print(expected_list(c).head(3))

# ------------------------------------------------------------------ 5 ----
step("5. coocur, prop_coin, as_nodes")
counts = pd.DataFrame({"a": [2, 0, 1, 3, 0, 1], "b": [1, 1, 0, 2, 2, 0], "c": [0, 2, 2, 1, 1, 1]})
print(coocur(counts))
print(prop_coin(c, margin=1).iloc[:3, :3], "<- % of the row species")
print(as_nodes(c, frequency=True, percentages=True, language="es").head(3))

# ------------------------------------------------------------------ 6 ----
step("6. all_net: the whole pipeline in one call")
finches = pd.read_csv(os.path.join(here, "sanderson", "finches.tsv"), sep="\t")
net = all_net(galapagos, procedures=["Jaccard", "Haberman"], lwidth="Jaccard", criteria="Jaccard",
              min_link=0.5, nodes=finches[["name", "type"]], color="type", layout="pc",
              main="Galapagos finches (Jaccard >= 0.5)", language="en")
print(net)                                            # first rows of nodes and links
net.summary()
net.plot(directory=os.path.join(out, "review_net")) if show else None

# ------------------------------------------------------------------ 7 ----
step("7. netCoin.add: a modified copy, the original is untouched")
other = net.add(main="Same network, other title", size="frequency", color=None)
print(net.options["main"], "|", other.options["main"])
print(other.options.get("nodeColor"), "vs", net.options.get("nodeColor"))
os.makedirs(out, exist_ok=True)
other.plot(directory=os.path.join(out, "review_net2")) if show else None
print("\nHTML files are written to", out, "when --open is given")
