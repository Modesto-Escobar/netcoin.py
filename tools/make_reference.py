"""Write docs/reference.md from the docstrings of the public API.

    python tools/make_reference.py
"""
from __future__ import annotations

import inspect
import re
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")
import netcoin  # noqa: E402

AREAS = {
    "The network and its page": ["netCoin", "multigraphCreate", "bar_coin", "BarCoin", "time_coin", "TimeCoin",
                                 "pie_coin", "PieCoin", "inc_time", "mobile_edges"],
    "Coincidences": ["dichotomize", "coin", "Coin", "coocur", "Cooc", "sim", "edge_list", "expected_list",
                     "prop_coin", "as_nodes", "all_net"],
    "Surveys": ["sur_coin", "cob_coin", "net_corr", "d_net_corr", "log_coin", "glm_coin", "caring",
                "caring_create_graphs", "caring_read_file"],
    "Scattergram and clusters": ["sur_scat", "add_axes", "add_clusters", "replace_clusters", "layouts_align",
                                 "add_communities", "cluster_match", "cluster_agree", "looking4clusters",
                                 "Looking4Clusters"],
    "Layouts": ["layout_mca", "layout_pca", "layout_fact", "layout_circle", "layout_grid"],
    "igraph and files": ["to_igraph", "from_igraph", "cal_centr", "save_pajek", "save_graphml"],
}

SECTION = re.compile(r"^(Parameters|Returns|Attributes|Examples|Raises|Notes|See Also)\n-{3,}$", re.M)


def entries(text: str) -> list[tuple[str, str]]:
    """``name : type`` lines followed by indented descriptions."""
    out: list[tuple[str, list[str]]] = []
    for line in text.strip("\n").split("\n"):
        if line and not line.startswith(" "):
            out.append((line, []))
        elif out:
            out[-1][1].append(line.strip())
    return [(head, " ".join(x for x in body if x)) for head, body in out]


def markdown(text: str) -> str:
    """reStructuredText of the docstrings (roles, double backticks) as markdown."""
    text = re.sub(r":(?:func|class|meth|mod|attr):`~?([^`]+)`", r"`\1`", text)
    return re.sub(r"``([^`]+)``", r"`\1`", text)


def convert(doc: str) -> str:
    parts = SECTION.split(markdown(doc))
    lines = [parts[0].strip()]
    for title, body in zip(parts[1::2], parts[2::2]):
        lines.append(f"**{title}**")
        if title == "Notes":
            lines.append(" ".join(body.split()))
            continue
        if title == "Examples":
            lines.append("```python\n" + body.strip("\n") + "\n```")
            continue
        items = []
        for head, text in entries(body):
            name, _, kind = head.partition(" : ")
            label = f"`{name}`" + (f" ({kind})" if kind else "")
            items.append(f"- {label}" + (f": {text}" if text else ""))
        lines.append("\n".join(items))
    return "\n\n".join(lines)


def public():
    found = {}
    for name in dir(netcoin):
        obj = getattr(netcoin, name)
        if name.startswith("_") or not (inspect.isfunction(obj) or inspect.isclass(obj)):
            continue
        if getattr(obj, "__module__", "").split(".")[0] == "netcoin":
            found[name] = obj
    return found


def build() -> str:
    found = public()
    listed = [n for names in AREAS.values() for n in names]
    left = sorted(set(found) - set(listed))
    areas = dict(AREAS)
    if left:
        areas["Other"] = left
    out = ["# Reference", "",
           "Generated from the docstrings by `python tools/make_reference.py`; do not edit it by hand.", ""]
    for area, names in areas.items():
        out += [f"## {area}", ""]
        for name in names:
            obj = found.get(name)
            if obj is None:
                continue
            try:
                signature = str(inspect.signature(obj))
            except (TypeError, ValueError):
                signature = "(...)"
            kind = "class" if inspect.isclass(obj) else "def"
            out += [f"### `{name}`", "", f"```python\n{kind} {name}{signature}\n```", "",
                    convert(inspect.getdoc(obj) or ""), ""]
    return "\n".join(out)


def main() -> None:
    text = build()
    (ROOT / "docs" / "reference.md").write_text(text, encoding="utf-8")
    print(f"docs/reference.md: {text.count('### `')} entries")


if __name__ == "__main__":
    main()
