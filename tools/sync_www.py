"""Copy the web assets of the R package rD3plot into ``netcoin/www``.

The JavaScript/CSS engine is not written here: it is rD3plot's. Run this script
to bring ``netcoin/www`` in line with a given checkout of rD3plot (galleries,
tutorials and other pages that netcoin does not use are left out)::

    python tools/sync_www.py --source D:/Documentos/Git/rD3plot
    python tools/sync_www.py --source D:/Documentos/Git/rD3plot --check   # only compare

It writes ``netcoin/www/SOURCE.json`` with the rD3plot version and commit and the
SHA-256 of every file, which ``tests/test_www_sync.py`` verifies.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WWW = ROOT / "netcoin" / "www"

# What netcoin uses, by purpose (paths relative to rD3plot/inst/www)
GROUPS = {
    "network": [
        "reset.css", "styles.css", "Roboto-Regular-latin.woff2", "Roboto-Regular-latin-ext.woff2",
        "d3.min.js", "jspdf.min.js", "jszip.min.js", "iro.min.js", "images.js", "colorScales.js",
        "functions.js", "en.js", "es.js", "ca.js", "network.js", "help_es.html",
    ],
    "barplot": ["barplot.js"],
    "timeline": ["timeline.js"],
    "pie": ["pie.js"],
    "multigraph": ["multigraph.js", "polygraph.js"],
    "statistics": [
        "statistics.css", "stdlib.bundle.js", "statistics_en.js", "statistics_es.js",
        "statistics_ca.js", "stats-icons",
    ],
    "page": ["template.html"],
}
# Files of netcoin/www that are not rD3plot's and must survive a sync
OWN_FILES = {"include.html", "SOURCE.json"}


def wanted_files() -> list[str]:
    names: list[str] = []
    for group in GROUPS.values():
        names.extend(group)
    return names


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expand(base: Path, name: str) -> list[Path]:
    path = base / name
    if path.is_dir():
        return sorted(p for p in path.rglob("*") if p.is_file())
    return [path]


def git(source: Path, *args: str) -> str | None:
    try:
        out = subprocess.run(["git", "-C", str(source), *args], capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def source_info(source: Path) -> dict:
    description = (source / "DESCRIPTION").read_text(encoding="utf-8")
    field = lambda key: (re.search(rf"^{key}:\s*(.+)$", description, re.M) or [None, None])[1]
    status = git(source, "status", "--porcelain", "--", "inst/www")
    return {
        "package": field("Package"),
        "version": field("Version"),
        "date": field("Date"),
        "commit": git(source, "rev-parse", "--short", "HEAD"),
        "dirty": bool(status),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default=os.environ.get("RD3PLOT_DIR", str(ROOT.parent / "rD3plot")),
                        help="rD3plot checkout (default: $RD3PLOT_DIR or ../rD3plot)")
    parser.add_argument("--check", action="store_true", help="do not copy, only report the differences")
    args = parser.parse_args(argv)

    source = Path(args.source)
    base = source / "inst" / "www"
    if not base.is_dir():
        print(f"error: {base} not found", file=sys.stderr)
        return 2

    files: dict[str, Path] = {}
    for name in wanted_files():
        for path in expand(base, name):
            if not path.exists():
                print(f"error: {path} is not in rD3plot", file=sys.stderr)
                return 2
            files[path.relative_to(base).as_posix()] = path

    stale = [p.relative_to(WWW).as_posix() for p in sorted(WWW.rglob("*")) if p.is_file()
             and p.relative_to(WWW).as_posix() not in files and p.name not in OWN_FILES]
    changed = [n for n, p in files.items() if not (WWW / n).is_file() or sha256(WWW / n) != sha256(p)]

    info = source_info(source)
    print(f"rD3plot {info['version']} ({info['commit']}{', modified' if info['dirty'] else ''})")
    print(f"{len(files)} files; {len(changed)} differ; {len(stale)} not used by netcoin")
    if args.check:
        for n in changed:
            print("  differs:", n)
        for n in stale:
            print("  not used:", n)
        return 1 if changed or stale else 0

    for n in stale:
        (WWW / n).unlink()
    for n in changed:
        (WWW / n).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(files[n], WWW / n)
    for d in sorted((p for p in WWW.rglob("*") if p.is_dir()), reverse=True):
        if not any(d.iterdir()):
            d.rmdir()

    manifest = {**info, "files": {n: sha256(WWW / n) for n in sorted(files)}}
    (WWW / "SOURCE.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("updated netcoin/www and SOURCE.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
