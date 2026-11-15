# -*- coding: utf-8 -*-
"""Bar, timeline and pie charts with the data of the netCoin package.

    python python_examples/review_charts.py           # writes the three pages
    python python_examples/review_charts.py --open    # and opens them in the browser

The same calls in R are in review_charts.R.
"""
import os
import re
import sys

import pandas as pd

# use the netcoin of this repository, even if another copy is installed in the environment
_f = globals().get("__file__", os.path.join(os.getcwd(), "python_examples", "x.py"))
_f = _f[1:] if re.match(r"[/\\][A-Za-z]:", _f) else _f
here = os.path.dirname(os.path.abspath(_f))
sys.path.insert(0, os.path.join(here, ".."))

from netcoin import bar_coin, coin, inc_time, mobile_edges, pie_coin, time_coin

show = "--open" in sys.argv
out = os.path.join(here, "tmp")
data = os.path.join(here, "..", "tests", "reference", "inputs")
pd.set_option("display.width", 160, "display.max_columns", 12)


def finish(graph, name):
    if show:
        graph.plot(os.path.join(out, name))        # writes the page and opens it
    else:
        graph._create(os.path.join(out, name))


# ------------------------------------------------------------------ bars ----
print("== bar_coin: survey of 1000 people (ESS), weighted")
ess = pd.read_csv(os.path.join(data, "ess_cweight.csv"))
bar = bar_coin(ess, ["Gender", "Age", "Social participation", "Political participation"],
               dichotomies=["Social participation", "Political participation"], value_dicho="Yes",
               weight="cweight", expected=True, confidence=True, significance=True,
               select="Gender:Male")
print(bar.nodes.head(4), "\n", bar.links.head(3), "\n", bar.options)
finish(bar, "review_bar")

# ------------------------------------------------------------- timeline ----
print("\n== time_coin: sociologists and their works")
soc = pd.read_csv(os.path.join(data, "soc.csv"))
works = pd.read_csv(os.path.join(data, "wrk.csv"))
timeline = time_coin(soc, start="birth", end="death", group="school", text="label",
                     main="Sociologists", events=works, event_names="works", event_period="name",
                     event_time="date", event_color="label")
print(timeline.periods[["name", "birth", "death", "school"]].head(3))
timeline.summary()
finish(timeline, "review_time")

print("\n== inc_time: who is alive in each year (first rows and columns)")
print(inc_time(soc, "name", "birth", "death").iloc[:3, :4])
print("== mobile_edges: pairs born within 10 years of each other")
print(mobile_edges(soc, "name", "birth", 10).head(3))

# ------------------------------------------------------------------ pie ----
print("\n== pie_coin: Galapagos finches")
galapagos = pd.read_csv(os.path.join(data, "galapagos.csv"))
pie = pie_coin(coin(galapagos), main="Galapagos finches")
print("pies:", pie.data.shape, "(row, column, part: both / X only / neither / Y only)")
print(pie.data[0, 1], "<- finches 1 and 2: cases of both, only 1, neither, only 2")
finish(pie, "review_pie")
print("\nPages written to", out, "(add --open to see them)")
