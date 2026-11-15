"""``multigraphCreate``: default, parallel and frame modes."""
import json
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from netcoin import multigraphCreate, netCoin
from netcoin.multigraph import evolvingNetwork
from netcoin.utils import PROJECT_MARK
from conftest import EXPECTED, read_expected, read_input
from test_html_output import _browser, _payload

NAMES = ["first", "second", "third"]


def _nets():
    nets = {}
    for k in NAMES:
        nodes = read_input(f"evo_nodes_{k}", strings=True)
        nodes["v"] = nodes["v"].astype(float)
        links = read_input(f"evo_links_{k}", strings=True)
        links["w"] = links["w"].astype(float)
        nets[k] = netCoin(nodes, links, main=f"Frame {k}", repulsion=50 if k == "third" else 25)
    return nets


# ---------------------------------------------------------------- frame ----

def test_evolving_network_matches_r():
    net = evolvingNetwork(_nets(), frame=1, speed=30, loop=True, lineplots=["v", "degree", "bogus"])
    nodes = read_expected("evo_result_nodes", strings=True)
    assert net.nodes["name"].tolist() == ["a", "b", "c", "d"]
    for col in ["v", "g"]:
        got = [x if isinstance(x, str) else ("" if pd.isna(x) else str(int(x))) for x in net.nodes[col]]
        want = ["" if x == "NA" else x for x in nodes[col]]
        assert got == want, col

    links = read_expected("evo_result_links")
    for col in ["w", "_frame_"]:
        assert net.links[col].tolist() == links[col].tolist(), col
    assert net.links["Source"].tolist() == links["Source"].tolist()
    assert net.links["Target"].tolist() == links["Target"].tolist()

    ref = json.loads((EXPECTED / "evo_result_options.json").read_text(encoding="utf-8"))
    options = net.options
    assert options["frames"] == ref["frames"]
    assert options["frame"] == ref["frame"][0] and options["speed"] == ref["speed"][0]
    assert options["loop"] is True
    assert options["lineplotsKeys"] == ref["lineplotsKeys"]
    assert options["main"] == ref["main"]
    assert options["repulsion"] == ref["repulsion"]


def test_evolving_network_needs_two_netcoin_graphs_with_the_same_columns():
    nets = _nets()
    with pytest.raises(ValueError, match="only one"):
        evolvingNetwork({"a": nets["first"]})
    with pytest.raises(ValueError, match="netCoin"):
        evolvingNetwork({"a": nets["first"], "b": "folder"})
    odd = netCoin({"name": ["a", "b"], "other": [1, 2]}, {"Source": ["a"], "Target": ["b"]})
    with pytest.raises(ValueError, match="node columns"):
        evolvingNetwork({"a": nets["first"], "b": odd})


def test_frame_mode_writes_one_page_with_the_frames(tmp_path):
    out = multigraphCreate(list(_nets().values()), names=NAMES, mode="frame",
                           directory=str(tmp_path / "evo"), show=False, frame=2, speed=10, loop=True)
    data = _payload(Path(out) / "index.html")
    assert data["options"]["frames"] == NAMES and data["options"]["frame"] == 2
    assert "_frame_" in data["linknames"]
    assert len(data["nodes"][0]) == 4


def test_frame_mode_with_a_single_graph_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="only one"):
        multigraphCreate([_nets()["first"]], mode="frame", directory=str(tmp_path / "x"), show=False)


# -------------------------------------------------------------- default ----

def test_default_mode_writes_a_selector_and_one_folder_per_graph(tmp_path):
    nets = _nets()
    out = Path(multigraphCreate(list(nets.values()), names=NAMES, directory=str(tmp_path / "multi"),
                                show=False))
    assert _payload(out / "index.html") == {"options": {"names": NAMES}}
    assert (out / "scripts" / "multigraph.js").is_file()
    for i, name in enumerate(NAMES):
        index = out / "data" / str(i) / "index.html"
        html = index.read_text(encoding="utf-8")
        assert PROJECT_MARK in html and f"<title>{name}</title>" in html
        assert _payload(index)["options"]["multigraph"] == {"idx": i, "names": NAMES}
    for net in nets.values():                      # the caller's objects are not modified
        assert "multigraph" not in net.options


def test_graphs_can_be_given_as_a_dict_and_are_named_automatically(tmp_path):
    nets = list(_nets().values())
    out = Path(multigraphCreate({"A": nets[0], "B": nets[1]}, directory=str(tmp_path / "d"), show=False))
    assert _payload(out / "index.html")["options"]["names"] == ["A", "B"]
    out = Path(multigraphCreate(nets, directory=str(tmp_path / "auto"), show=False))
    assert _payload(out / "index.html")["options"]["names"] == ["graph1", "graph2", "graph3"]


def test_names_of_the_wrong_length_are_ignored_with_a_warning(tmp_path):
    with pytest.warns(UserWarning, match="names"):
        out = Path(multigraphCreate(list(_nets().values()), names=["only"],
                                    directory=str(tmp_path / "w"), show=False))
    assert _payload(out / "index.html")["options"]["names"][0] == "graph1"


def test_a_folder_that_already_holds_a_graph_is_included(tmp_path):
    inner = tmp_path / "inner"
    nets = _nets()
    nets["first"]._create(str(inner))
    out = Path(multigraphCreate([nets["second"], str(inner)], names=["fresh", "old"],
                                directory=str(tmp_path / "m"), show=False))
    assert _payload(out / "index.html")["options"] == {"names": ["fresh", "old"], "external": [False, True]}
    assert (out / "data" / "1" / "index.html").is_file()


def test_unsupported_objects_are_skipped_with_a_warning(tmp_path):
    nets = list(_nets().values())
    with pytest.warns(UserWarning, match="Not supported"):
        out = Path(multigraphCreate([nets[0], 42, nets[1]], directory=str(tmp_path / "s"), show=False))
    assert _payload(out / "index.html")["options"]["names"] == ["graph1", "graph2"]


def test_wrong_arguments(tmp_path):
    with pytest.raises(ValueError, match="non-empty"):
        multigraphCreate([], directory=str(tmp_path / "e"), show=False)
    with pytest.raises(ValueError, match="mode"):
        multigraphCreate(list(_nets().values()), mode="sideways", directory=str(tmp_path / "e"), show=False)


def test_an_unrelated_folder_is_never_overwritten(tmp_path):
    keep = tmp_path / "keep"
    keep.mkdir()
    (keep / "file.txt").write_text("mine")
    with pytest.raises(FileExistsError):
        multigraphCreate(list(_nets().values()), directory=str(keep), show=False)
    assert (keep / "file.txt").read_text() == "mine"


# ------------------------------------------------------------- parallel ----

def test_parallel_mode_has_a_grid_and_removes_controls_when_it_is_crowded(tmp_path):
    nets = list(_nets().values())
    out = Path(multigraphCreate(nets, names=NAMES, mode="parallel", mfrow=(1, 3),
                                directory=str(tmp_path / "p"), show=False))
    assert _payload(out / "index.html")["options"] == {"mfrow": [1, 3], "n": 3}
    inner = _payload(out / "multiGraph" / "data" / "0" / "index.html")
    assert "controls" not in inner["options"]
    assert all("controls" in n.options for n in nets)          # the originals keep them

    out2 = Path(multigraphCreate(nets[:2], mode="parallel", directory=str(tmp_path / "p2"), show=False))
    assert _payload(out2 / "index.html")["options"]["mfrow"] == [1, 2]
    assert "controls" in _payload(out2 / "multiGraph" / "data" / "0" / "index.html")["options"]


def test_parallel_mode_validates_mfrow(tmp_path):
    with pytest.raises(ValueError, match="mfrow"):
        multigraphCreate(list(_nets().values()), mode="parallel", mfrow=(0, 2),
                         directory=str(tmp_path / "x"), show=False)


# -------------------------------------------------------- real browser ----

def _dom(page: Path) -> str:
    return subprocess.run(
        [_browser(), "--headless=new", "--disable-gpu", "--virtual-time-budget=8000", "--dump-dom",
         page.as_uri()], capture_output=True, text=True, encoding="utf-8", timeout=120).stdout


@pytest.mark.skipif(_browser() is None, reason="no Chromium-based browser found")
def test_browser_default_mode_draws_the_first_graph(tmp_path):
    out = Path(multigraphCreate(list(_nets().values()), names=NAMES,
                                directory=str(tmp_path / "m"), show=False))
    dom = _dom(out / "data" / "0" / "index.html")
    assert "Frame first" in dom and dom.count("<line") >= 2


@pytest.mark.skipif(_browser() is None, reason="no Chromium-based browser found")
def test_browser_frame_mode_shows_the_frame_slider(tmp_path):
    out = Path(multigraphCreate(list(_nets().values()), names=NAMES, mode="frame", speed=0,
                                directory=str(tmp_path / "e"), show=False))
    assert "Frame: 1/3" in _dom(out / "index.html")      # speed=0: it does not play
