"""``netcoin/www`` is rD3plot's engine: check it is intact and complete."""
import hashlib
import json
from pathlib import Path

import pytest

import netcoin
from netcoin import network

WWW = Path(netcoin.__file__).parent / "www"
SOURCE = json.loads((WWW / "SOURCE.json").read_text(encoding="utf-8"))


def test_source_is_recorded():
    assert SOURCE["package"] == "rD3plot" and SOURCE["version"]
    assert SOURCE["commit"], "run tools/sync_www.py from a git checkout of rD3plot"


@pytest.mark.parametrize("name", sorted(SOURCE["files"]))
def test_file_matches_its_recorded_hash(name):
    digest = hashlib.sha256((WWW / name).read_bytes()).hexdigest()
    assert digest == SOURCE["files"][name], (
        f"{name} was edited by hand; change it in rD3plot and run tools/sync_www.py")


def test_no_unlisted_files():
    present = {p.relative_to(WWW).as_posix() for p in WWW.rglob("*") if p.is_file()}
    assert present - set(SOURCE["files"]) == {"include.html", "SOURCE.json"}


def test_galleries_are_not_shipped():
    names = [p.name for p in WWW.rglob("*")]
    assert not [n for n in names if n.startswith(("gallery", "tutorial", "multipages"))]


@pytest.mark.parametrize("language", ["en", "es", "ca"])
def test_every_script_of_a_network_exists(language, tmp_path, monkeypatch):
    written = {}
    monkeypatch.setattr(network, "createHTML",
                        lambda directory, styles, scripts, data, title=None: written.update(styles=styles, scripts=scripts))
    net = {"options": {"language": language, "statistics": True}, "nodes": {}, "links": {}}
    network.netCreate(net, str(tmp_path / "x"))
    assert written["scripts"][-2:] == [f"{language}.js", "network.js"]
    for name in written["styles"] + written["scripts"]:
        assert (WWW / name).is_file(), name
