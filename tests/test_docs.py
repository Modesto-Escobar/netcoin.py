"""The code of the documentation runs, and the public API is documented."""
import inspect
import re
import warnings
from pathlib import Path

import pytest

import netcoin

DOCS = Path(__file__).resolve().parents[1] / "docs"
GUIDES = ["guide.md", "scattergram.md"]


def blocks(path):
    text = path.read_text(encoding="utf-8")
    return [b for b in re.findall(r"```python\n(.*?)```", text, re.S) if not b.lstrip().startswith("# no-run")]


@pytest.mark.parametrize("guide", GUIDES)
def test_the_examples_of_the_guides_run(guide, tmp_path, monkeypatch):
    pytest.importorskip("igraph")
    pytest.importorskip("statsmodels")
    pytest.importorskip("sklearn")
    monkeypatch.chdir(tmp_path)
    (tmp_path / "out").mkdir()
    namespace: dict = {}
    for i, code in enumerate(blocks(DOCS / guide)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                exec(compile(code, f"{guide}[{i}]", "exec"), namespace)
            except Exception as error:                                     # say which block
                raise AssertionError(f"{guide}, block {i}:\n{code}\n-> {type(error).__name__}: {error}") from error


def public():
    for name in sorted(dir(netcoin)):
        obj = getattr(netcoin, name)
        if name.startswith("_") or not (inspect.isfunction(obj) or inspect.isclass(obj)):
            continue
        if getattr(obj, "__module__", "").split(".")[0] == "netcoin":
            yield name, obj


@pytest.mark.parametrize("name,obj", list(public()))
def test_every_public_function_says_what_it_takes_and_gives(name, obj):
    doc = inspect.getdoc(obj) or ""
    assert doc, f"{name} has no docstring"
    assert "Parameters" in doc or "Attributes" in doc, f"{name} does not describe its parameters"


def test_the_reference_is_the_one_the_docstrings_make():
    import importlib.util

    spec = importlib.util.spec_from_file_location("make_reference", DOCS.parent / "tools" / "make_reference.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    reference = (DOCS / "reference.md").read_text(encoding="utf-8")
    missing = [name for name, _ in public() if f"### `{name}" not in reference]
    assert not missing, f"docs/reference.md lacks {missing}"
    assert module.build() == reference, "docs/reference.md is out of date: python tools/make_reference.py"
