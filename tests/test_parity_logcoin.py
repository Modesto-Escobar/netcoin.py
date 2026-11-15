"""``log_coin`` against netCoin 2.1.49."""
import json
import warnings

import pytest

from netcoin import log_coin
from conftest import EXPECTED, compare_frame, compare_options, read_expected, read_input

ESS = read_input("ess_cweight", strings=True).astype({"cweight": float})
LG = read_input("log_cat5", strings=False)
LG["N"] = LG["N"].astype(int)
VARS3 = ["Gender", "Age", "Social participation"]
VARS4 = VARS3 + ["Political participation"]
EXO4 = ["Age", "Social participation", "Political participation", "Gender"]

# python call of each R call (the R calls are in tests/reference/make_reference.R)
CASES = {
    "default": lambda: log_coin(ESS, VARS3),
    "twotails": lambda: log_coin(ESS, ["Gender", "Age", "Political participation"], frequency=True,
                                 percentage=True, twotails=True, pmax=0.1),
    "order3": lambda: log_coin(ESS, EXO4, exogenous="Gender", order=3, pmax=0.2),
    "order3_freq": lambda: log_coin(ESS, EXO4, exogenous="Gender", order=3, pmax=0.2, frequency=True,
                                    percentage=True),
    "nofirst": lambda: log_coin(ESS, VARS3, no_first_cat="Gender"),
    "nofirst_two": lambda: log_coin(ESS, VARS4, no_first_cat=["Gender", "Age"], twotails=True),
    "weight": lambda: log_coin(ESS, VARS3, weight="cweight"),
    "spanish": lambda: log_coin(ESS, VARS3, language="es", frequency=True),
    "four": lambda: log_coin(ESS, VARS4, pmax=0.1),
    "five": lambda: log_coin(LG, ["A", "B", "C", "D var", "E"], pmax=0.1),
    "five_two": lambda: log_coin(LG, ["A", "B", "C", "D var", "E"], twotails=True, pmax=0.2, frequency=True),
    "numeric": lambda: log_coin(LG, ["N", "B", "C"], weight="W", pmax=0.3),
    "exo_two": lambda: log_coin(LG, ["A", "B", "C", "E"], exogenous=["A", "E"], order=3, pmax=0.3,
                                no_first_cat="B"),
    "options": lambda: log_coin(ESS, VARS3, color="cat", lwidth="z value", main="Log-linear", degree_filter=0),
}


def _ref(case):
    return json.loads((EXPECTED / f"net_lc_{case}.json").read_text(encoding="utf-8"))


def _run(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return CASES[case]()


@pytest.mark.parametrize("case", CASES)
def test_nodes_and_options_match_r(case):
    net, ref = _run(case), _ref(case)
    compare_frame(net.nodes, read_expected(f"net_lc_{case}_nodes"), f"{case}/nodes", rtol=1e-6)
    compare_options(net.options, ref["options"], case)


@pytest.mark.parametrize("case", CASES)
def test_links_match_r(case):
    net, ref = _run(case), _ref(case)
    assert net.lcount() == ref["nlinks"], f"{case}: links"
    if ref["nlinks"]:
        compare_frame(net.links, read_expected(f"net_lc_{case}_links"), f"{case}/links", rtol=1e-6)
