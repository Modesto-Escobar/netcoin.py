"""bar_coin, time_coin, inc_time, mobile_edges and pie_coin against netCoin 2.1.49."""
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from netcoin import coin, multigraphCreate
from netcoin.barplot import BarCoin, bar_coin
from netcoin.pie import PieCoin, pie_coin
from netcoin.timeline import TimeCoin, inc_time, mobile_edges, time_coin
from netcoin.webgraph import category_colors, col2hex, is_color, r_json, records
from conftest import EXPECTED, read_expected, read_input
from test_html_output import _browser, _payload


def _options(name):
    return json.loads((EXPECTED / name).read_text(encoding="utf-8"))


def _same_table(actual: pd.DataFrame, expected: pd.DataFrame, label: str, rtol=1e-8):
    assert list(actual.columns) == list(expected.columns), f"{label}: columns"
    assert len(actual) == len(expected), f"{label}: {len(actual)} rows, R has {len(expected)}"
    for col in expected.columns:
        a, e = actual[col].reset_index(drop=True), expected[col].reset_index(drop=True)
        if pd.api.types.is_numeric_dtype(e) and pd.api.types.is_numeric_dtype(a):
            np.testing.assert_allclose(a.astype(float), e.astype(float), rtol=rtol, atol=1e-10,
                                       equal_nan=True, err_msg=f"{label}/{col}")
        else:
            assert a.astype(str).tolist() == e.astype(str).tolist(), f"{label}/{col}"


def _same_options(actual: dict, expected: dict, label: str):
    assert set(actual) == set(expected), f"{label}: {sorted(set(actual) ^ set(expected))}"
    for key, value in expected.items():
        if isinstance(value, float):
            assert actual[key] == pytest.approx(value), f"{label}/{key}"
        else:
            assert actual[key] == value, f"{label}/{key}"


# ------------------------------------------------------------------- bars ----

def _ess():
    return read_input("ess_cweight")


BAR_CASES = {
    "categorical": lambda: bar_coin(_ess(), ["Gender", "Age"], dichotomies="_none"),
    "full": lambda: bar_coin(
        _ess(), ["Gender", "Age", "Social participation", "Political participation"],
        dichotomies=["Social participation", "Political participation"], value_dicho="Yes",
        expected=True, confidence=True, significance=True, level=0.9, percentages=True,
        language="es", sort="incidences", select="Gender:Male"),
    "weight": lambda: bar_coin(_ess(), ["Gender", "Age", "cweight"], dichotomies="_none",
                               weight="cweight", criteria="hyp", support=5, max_link=0.6,
                               confidence=True),
    "dichotomies": lambda: bar_coin(read_input("galapagos")),
    "nodes": lambda: bar_coin(read_input("galapagos"), nodes=read_input("bar_nodes_user", strings=True),
                              label="label", color="group", expected=True, decreasing=False,
                              sort="label", scalebar=True, note="A note", cex=1.2,
                              default_color="firebrick"),
    "commonlabel": lambda: bar_coin(_ess(), ["Gender", "Age"], dichotomies="_none",
                                    commonlabel="Gender"),
}


@pytest.mark.parametrize("case", BAR_CASES)
def test_bar_coin_matches_r(case):
    if case == "nodes":
        with pytest.warns(UserWarning, match="not present among incidences"):
            bar = BAR_CASES[case]()
    else:
        bar = BAR_CASES[case]()
    assert isinstance(bar, BarCoin)
    _same_table(bar.nodes, read_expected(f"bar_{case}_nodes"), f"{case}/nodes")
    _same_table(bar.links, read_expected(f"bar_{case}_links"), f"{case}/links")
    _same_options(bar.options, _options(f"bar_{case}_options.json"), case)


def test_bar_coin_links_must_be_nodes():
    nodes = pd.DataFrame({"name": ["a"], "incidences": [1]})
    links = pd.DataFrame({"Source": ["a"], "Target": ["zzz"], "coincidences": [1]})
    with pytest.raises(ValueError, match="not defined"):
        BarCoin(nodes, links)


def test_bar_coin_sort_must_be_a_column():
    nodes = pd.DataFrame({"name": ["a"], "incidences": [1]})
    links = pd.DataFrame({"Source": [], "Target": [], "coincidences": []})
    with pytest.warns(UserWarning, match="sort"):
        bar = BarCoin(nodes, links, sort="nope")
    assert "order" not in bar.options


def test_bar_coin_without_dichotomies_or_variables_fails():
    with pytest.raises(ValueError):
        bar_coin(pd.DataFrame({"a": [1, 2]}), [], dichotomies="_none")


# ------------------------------------------------------------- timelines ----

def test_time_coin_periods_match_r():
    soc = read_input("soc", strings=True)
    soc[["birth", "death"]] = soc[["birth", "death"]].astype(float)
    t = time_coin(soc, start="birth", end="death", group="school", text="label", main="Sociologists",
                  info="birth_country", language="es", cex=1.1)
    assert isinstance(t, TimeCoin) and t.events is None
    _same_table(t.periods, read_expected("time_periods_periods"), "periods")
    _same_options(t.options, _options("time_periods_options.json"), "periods")


def test_time_coin_events_match_r():
    soc = read_input("soc", strings=True)
    soc[["birth", "death"]] = soc[["birth", "death"]].astype(float)
    wrk = read_input("wrk", strings=True)
    wrk["date"] = wrk["date"].astype(float)
    t = time_coin(soc, start="birth", end="death", events=wrk, event_names="works",
                  event_period="name", event_time="date", event_color="label", event_shape="name")
    _same_table(t.events, read_expected("time_events_events"), "events")
    _same_options(t.options, _options("time_events_options.json"), "events")


def test_time_coin_datetimes_are_seconds_and_flag_posixct():
    posix = read_input("posix", strings=True)
    posix["start"] = pd.to_datetime(posix["start"])
    posix["end"] = pd.to_datetime(posix["end"])
    t = time_coin(posix)
    _same_table(t.periods, read_expected("time_posix_periods"), "posix")
    _same_options(t.options, _options("time_posix_options.json"), "posix")


def test_time_coin_needs_its_columns():
    with pytest.raises(ValueError, match="present"):
        time_coin(pd.DataFrame({"name": ["a"], "start": [1]}))


def test_time_coin_events_are_ordered_by_time():
    periods = pd.DataFrame({"name": ["a"], "start": [1], "end": [9]})
    events = pd.DataFrame({"name": ["e2", "e1"], "period": ["a", "a"], "date": [5, 2]})
    assert time_coin(periods, events=events).events["name"].tolist() == ["e1", "e2"]


def test_inc_time_matches_r():
    soc = read_input("soc", strings=True)
    got = inc_time(soc, "name", "birth", "death")
    expected = read_expected("inctime", strings=True).set_index("year")
    assert [str(y) for y in got.index] == list(expected.index)
    assert list(got.columns) == list(expected.columns)
    np.testing.assert_array_equal(got.to_numpy(dtype=float), expected.to_numpy(dtype=float))


def test_inc_time_open_end_goes_on_and_missing_beginning_is_missing():
    data = pd.DataFrame({"name": ["a", "b", "c"], "birth": [2000, 2001, np.nan], "death": [2002, np.nan, 2003]})
    got = inc_time(data)
    assert got.index.tolist() == [2000, 2001, 2002, 2003]
    assert got["a"].tolist() == [1, 1, 1, 0]
    assert got["b"].tolist() == [0, 1, 1, 1]
    assert got["c"].isna().all()


@pytest.mark.parametrize("difference", [0, 10])
def test_mobile_edges_match_r(difference):
    soc = read_input("soc", strings=True)
    got = mobile_edges(soc, "name", "birth", difference)
    _same_table(got, read_expected(f"mobile_{difference}"), f"mobile_{difference}")


# ------------------------------------------------------------------ pies ----

def _pie_dict(pie: PieCoin) -> dict:
    out = pie._data()
    out["data"] = np.asarray(out["data"]).tolist()
    if "dataw" in out:
        out["dataw"] = np.asarray(out["dataw"]).tolist()
    for key, names in (("nodes", "nodenames"), ("links", "linknames")):
        if key in out:
            out[key] = dict(zip(out.pop(names), out[key]))
    return out


def _clean(x):
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_clean(v) for v in x]
    if isinstance(x, float) and np.isnan(x):
        return None
    return x


def _compare(a, b, path=""):
    if isinstance(b, dict):
        assert isinstance(a, dict) and set(a) == set(b), f"{path}: {sorted(set(a) ^ set(b)) if isinstance(a, dict) else a}"
        for k in b:
            _compare(a[k], b[k], f"{path}/{k}")
    elif isinstance(b, list):
        assert isinstance(a, list) and len(a) == len(b), f"{path}: length"
        for i, (x, y) in enumerate(zip(a, b)):
            _compare(x, y, f"{path}[{i}]")
    elif isinstance(b, float) and a is not None:
        assert a == pytest.approx(b, rel=1e-8, abs=1e-10), path
    else:
        assert a == b, f"{path}: {a!r} != {b!r}"


PIE_CASES = {
    "default": lambda c: pie_coin(c),
    "options": lambda c: pie_coin(c, expected=False, abline=3, main="Pies", note="n", language="es",
                                  colors=["red", "#00F", "lightblue", "gray 50"], cex=1.2,
                                  show_legend=False, help_on=True, help="hi"),
    "links_column": lambda c: pie_coin(c, nodes=read_input("pie_nodes", strings=True).assign(size=[1, 2, 3]),
                                       links=read_input("pie_links"), lcolor="Haberman"),
    "links_color": lambda c: pie_coin(c, links=read_input("pie_links"), lcolor="firebrick"),
    "links_numeric": lambda c: pie_coin(c, links=read_input("pie_links"),
                                        lcolor=list(range(1, len(read_input("pie_links")) + 1))),
}


@pytest.mark.parametrize("case", PIE_CASES)
def test_pie_coin_matches_r(case):
    expected = json.loads((EXPECTED / f"pie_{case}.json").read_text(encoding="utf-8"))
    pie = PIE_CASES[case](coin(read_input("galapagos")))
    assert isinstance(pie, PieCoin)
    got = _clean(_pie_dict(pie))
    # R writes length-one vectors as scalars
    for key in ("dim", "colnames", "rownames"):
        assert got[key] == expected[key], key
    # R writes a missing element as {}
    expected = {k: v for k, v in expected.items() if v != {}}
    _compare(got["options"], expected["options"], "options")
    _compare(got["data"], expected["data"], "data")
    if "dataw" in expected:
        _compare(got["dataw"], expected["dataw"], "dataw")
    else:
        assert "dataw" not in got
    for key in ("nodes", "links"):
        assert (key in expected) == (key in got), key
        if key in expected:
            exp = {k: v if isinstance(v, list) else [v] for k, v in expected[key].items()}
            _compare(got[key], exp, key)


def test_pie_coin_needs_a_coin():
    with pytest.raises(TypeError):
        pie_coin(read_input("galapagos"))


def test_pie_parts_add_up_to_n():
    c = coin(read_input("galapagos"))
    pie = pie_coin(c)
    assert np.allclose(pie.data.sum(axis=2), c.n)
    assert np.allclose(pie.dataw[~np.isnan(pie.dataw[:, :, 0])].sum(axis=1), c.n)


def test_pie_coin_warns_about_bad_arguments():
    c = coin(read_input("galapagos"))
    with pytest.warns(UserWarning, match="ablineX"):
        pie_coin(c, abline=99)
    with pytest.warns(UserWarning, match="colors"):
        pie_coin(c, colors=["red", "blue"])
    links = pd.DataFrame({"Source": ["zz"], "Target": ["yy"]})
    with pytest.warns(UserWarning, match="links"):
        pie_coin(c, links=links)


# ------------------------------------------------- colours and the JSON ----

@pytest.mark.parametrize("color,hexa", [
    ("red", "#ff0000"), ("Red", "#ff0000"), ("light blue", "#add8e6"), ("LightBlue", "#add8e6"),
    ("firebrick3", "#cd2626"), ("#1F77B4", "#1f77b4"), ("#abc", "#aabbcc"), ("#11223344", "#112233"),
    ("transparent", "#ffffff"), ("gray50", "#7f7f7f"), ("grey 50", "#7f7f7f"),
])
def test_col2hex_matches_r(color, hexa):
    assert col2hex(color) == hexa


def test_invalid_colours():
    assert not is_color("notacolour") and not is_color(["red", "nope"]) and not is_color([])
    assert is_color(["red", "#00F"])
    with pytest.raises(ValueError):
        col2hex("#12")


def test_category_colors():
    assert category_colors(["a", "b", "a"]) == ["#1f77b4", "#ff7f0e", "#1f77b4"]
    assert category_colors([1, 2, 21]) == ["#1f77b4", "#ff7f0e", "#1f77b4"]


def test_r_json_follows_rd3plot():
    assert r_json(10.705882352) == 10.71
    assert r_json(1234.5678) == 1235
    assert r_json(0.000123456) == 0.0001235
    assert r_json(3.0) == 3 and isinstance(r_json(3.0), int)
    assert r_json(float("nan")) is None and r_json(float("inf")) is None
    assert r_json("1990") == 1990 and r_json(" 2.5 ") == 2.5 and r_json("1e3") == 1000
    assert r_json("0x1A") == 26
    assert r_json("abc") == "abc" and r_json("12abc") == "12abc" and r_json("NA") == "NA"
    assert r_json("Inf") is None
    assert r_json({"a": np.int64(3), "b": [np.float64(0.5), None]}) == {"a": 3, "b": [0.5, None]}
    frame = pd.DataFrame({"x": [1.5, np.nan], "y": ["a", None]})
    assert records(frame) == [{"x": 1.5, "y": "a"}, {"x": None, "y": None}]


# ------------------------------------------------------------- the pages ----

def _bar():
    return bar_coin(_ess(), ["Gender", "Age"], dichotomies="_none", expected=True, confidence=True)


def test_pages_are_written_with_their_script(tmp_path):
    graphs = {"bar": (_bar(), "barplot.js"),
              "time": (time_coin(pd.DataFrame({"name": ["a", "b"], "start": [1, 2], "end": [5, 6]})), "timeline.js"),
              "pie": (pie_coin(coin(read_input("galapagos"))), "pie.js")}
    for key, (graph, script) in graphs.items():
        out = tmp_path / key
        graph._create(str(out))
        assert (out / "scripts" / script).is_file() and (out / "scripts" / "images.js").is_file()
        data = _payload(out / "index.html")
        assert data["options"]["language"] == "en"


def test_the_bar_page_holds_records_as_r_writes_them(tmp_path):
    _bar()._create(str(tmp_path / "b"))
    data = _payload(tmp_path / "b" / "index.html")
    assert data["nodes"][0] == {"name": "Gender:Female", "incidences": 534}
    assert set(data["links"][0]) == {"Source", "Target", "coincidences", "expected", "conf.L", "conf.U"}
    assert data["options"]["confidence"] == ["conf.L", "conf.U"]


def test_directory_argument_writes_the_page(tmp_path):
    bar_coin(_ess(), ["Gender", "Age"], dichotomies="_none", directory=str(tmp_path / "d"))
    assert (tmp_path / "d" / "index.html").is_file()


def test_charts_go_into_a_multigraph(tmp_path):
    bar = _bar()
    time = time_coin(pd.DataFrame({"name": ["a", "b"], "start": [1, 2], "end": [5, 6]}))
    pie = pie_coin(coin(read_input("galapagos")))
    out = Path(multigraphCreate([bar, time, pie], names=["bars", "time", "pies"],
                                directory=str(tmp_path / "m"), show=False))
    for i in range(3):
        assert _payload(out / "data" / str(i) / "index.html")["options"]["multigraph"]["idx"] == i
    assert "multigraph" not in bar.options          # the caller's object is untouched


def test_summaries(capsys):
    _bar().summary()
    assert "nodes and" in capsys.readouterr().out
    time_coin(pd.DataFrame({"name": ["a", "b"], "start": [1, 2], "end": [5, 6]})).summary()
    assert "2 categories" in capsys.readouterr().out


# -------------------------------------------------------- real browser ----

def _dom(page: Path) -> str:
    return subprocess.run(
        [_browser(), "--headless=new", "--disable-gpu", "--virtual-time-budget=6000", "--dump-dom",
         page.as_uri()], capture_output=True, text=True, encoding="utf-8", timeout=120).stdout


@pytest.mark.skipif(_browser() is None, reason="no Chromium-based browser found")
@pytest.mark.parametrize("kind", ["bar", "time", "pie"])
def test_browser_draws_the_chart(kind, tmp_path):
    if kind == "bar":
        graph = bar_coin(_ess(), ["Gender", "Age"], dichotomies="_none", expected=True,
                         select="Gender:Male")
        marker = "Gender:Male"
    elif kind == "time":
        soc = read_input("soc", strings=True)
        soc[["birth", "death"]] = soc[["birth", "death"]].astype(float)
        graph = time_coin(soc, start="birth", end="death", main="Sociologists timeline")
        marker = "Sociologists timeline"
    else:
        graph = pie_coin(coin(read_input("galapagos")), main="Finches pies")
        marker = "Finches pies"
    graph._create(str(tmp_path / kind))
    dom = _dom(tmp_path / kind / "index.html")
    assert marker in dom
    assert dom.count("<rect") + dom.count("<path") + dom.count("<circle") > 10, "nothing was drawn"
