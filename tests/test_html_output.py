"""``createHTML`` and the pages it writes."""
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from netcoin import all_net, netCoin
from netcoin.utils import PROJECT_MARK, createHTML
from conftest import read_input


def _payload(index: Path) -> dict:
    html = index.read_text(encoding="utf-8")
    raw = re.search(r'<script type="application/json" id="data">(.*?)</script>', html, re.S).group(1)
    return json.loads(raw)


def test_writes_a_marked_utf8_page(tmp_path):
    out = tmp_path / "net"
    createHTML(str(out), ["reset.css"], ["d3.min.js"], {"text": "Contemporáneos ñ"})
    html = (out / "index.html").read_text(encoding="utf-8")
    assert PROJECT_MARK in html
    assert 'href="styles/reset.css"' in html and 'src="scripts/d3.min.js"' in html
    assert (out / "styles" / "reset.css").is_file() and (out / "scripts" / "d3.min.js").is_file()
    assert _payload(out / "index.html") == {"text": "Contemporáneos ñ"}


def test_fonts_travel_with_styles_css(tmp_path):
    createHTML(str(tmp_path / "n"), ["styles.css"], [], None)
    assert (tmp_path / "n" / "styles" / "Roboto-Regular-latin.woff2").is_file()


def test_title_is_the_folder_name_with_either_separator(tmp_path):
    createHTML(str(tmp_path / "a" / "my graph"), [], [], None)
    assert "<title>my graph</title>" in (tmp_path / "a" / "my graph" / "index.html").read_text(encoding="utf-8")
    createHTML(str(tmp_path / "b") + "/sub/other", [], [], None)
    assert "<title>other</title>" in (tmp_path / "b" / "sub" / "other" / "index.html").read_text(encoding="utf-8")


def test_data_cannot_close_the_script_element(tmp_path):
    createHTML(str(tmp_path / "n"), [], [], {"x": "</script><b>boom</b>"})
    html = (tmp_path / "n" / "index.html").read_text(encoding="utf-8")
    assert html.count("</script>") == 1
    assert _payload(tmp_path / "n" / "index.html")["x"] == "</script><b>boom</b>"


def test_an_empty_directory_is_used(tmp_path):
    createHTML(str(tmp_path), [], [], {"a": 1})
    assert (tmp_path / "index.html").is_file()


def test_a_netcoin_page_is_replaced(tmp_path):
    createHTML(str(tmp_path / "n"), [], [], {"a": 1})
    (tmp_path / "n" / "stale.txt").write_text("old")
    createHTML(str(tmp_path / "n"), [], [], {"a": 2})
    assert not (tmp_path / "n" / "stale.txt").exists()
    assert _payload(tmp_path / "n" / "index.html") == {"a": 2}


def test_foreign_directories_are_never_deleted(tmp_path):
    mine = tmp_path / "documents"
    mine.mkdir()
    (mine / "thesis.docx").write_text("precious")
    with pytest.raises(FileExistsError):
        createHTML(str(mine), [], [], {"a": 1})
    assert (mine / "thesis.docx").read_text() == "precious"

    other = tmp_path / "site"
    other.mkdir()
    (other / "index.html").write_text("<html>not ours</html>")
    with pytest.raises(FileExistsError):
        createHTML(str(other), [], [], {"a": 1})
    assert (other / "index.html").read_text() == "<html>not ours</html>"


def test_language_selects_the_language_script(tmp_path):
    netCoin({"nombre": ["a", "b"]}, {"Source": ["a"], "Target": ["b"]},
            language="es", directory=str(tmp_path / "n"))
    assert (tmp_path / "n" / "scripts" / "es.js").is_file()
    assert not (tmp_path / "n" / "scripts" / "en.js").exists()


# ----------------------------------------------------------- real browser ----
def _browser():
    candidates = [
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
        "/usr/bin/google-chrome", "/usr/bin/chromium", "/usr/bin/chromium-browser",
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    ]
    for c in candidates:
        if Path(c).is_file():
            return c
    return shutil.which("chromium") or shutil.which("google-chrome") or shutil.which("msedge")


@pytest.mark.skipif(_browser() is None, reason="no Chromium-based browser found")
def test_the_network_is_drawn_by_the_browser(tmp_path):
    """Smoke test: the page loads rD3plot's scripts, draws every link and shows the title."""
    net = all_net(read_input("galapagos"), procedures="Jaccard", criteria="Jaccard", min_link=0.5,
                  main="Pinzones de Galápagos", language="es", directory=str(tmp_path / "net"))
    dom = subprocess.run(
        [_browser(), "--headless=new", "--disable-gpu", "--virtual-time-budget=5000", "--dump-dom",
         (tmp_path / "net" / "index.html").as_uri()],
        capture_output=True, text=True, encoding="utf-8", timeout=120).stdout
    assert "Pinzones de Galápagos" in dom
    assert dom.count("<line") >= net.lcount()
    assert "Número de escenarios: 17" in dom
