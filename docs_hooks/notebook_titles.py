"""Keep the notebook's own H1 in the browser tab while the sidebar shows nav keys.

``ignore_h1_titles`` stops mkdocs-jupyter from replacing the configured nav
title ("Python:" / "R:") with the notebook's first heading, but MkDocs Material
then also renders ``page.title`` (the bare nav key) as the page ``<title>``.
This hook recovers the notebook's first H1 into ``page.meta.title``, which
Material prefers for the ``<title>`` tag, so tabs stay descriptive while the
sidebar keeps the configured labels.
"""
import json

_H1_CACHE = {}


def _first_h1(src_path):
    notebook = json.load(open(src_path, encoding="utf-8"))
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "markdown":
            continue
        for line in "".join(cell.get("source", [])).splitlines():
            if line.startswith("# "):
                return line[2:].strip()
    return None


def on_page_content(html, page, **kwargs):
    if page.file.src_path.endswith(".ipynb"):
        path = page.file.abs_src_path
        if path not in _H1_CACHE:
            try:
                _H1_CACHE[path] = _first_h1(path)
            except (OSError, ValueError):
                _H1_CACHE[path] = None
        title = _H1_CACHE[path]
        if title:
            page.meta["title"] = title
