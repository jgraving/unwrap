"""Keep links portable between the source notebooks and the published guide."""

import logging
import posixpath
import re
from html import escape, unescape
from urllib.parse import unquote, urlsplit, urlunsplit

from mkdocs.utils import get_relative_url


# nbconvert emits quoted href attributes. Restrict substitutions to anchors
# so the original code, images, and saved output markup survive the build.
_ANCHOR = re.compile(r"""(<a\b[^<>]*?\shref=)(["'])(.*?)\2""", re.IGNORECASE)
_LOG = logging.getLogger("mkdocs.hooks.notebook_links")
_FIGURE_NOTE = "notebooks/BEETLE_FIGURE.md"
_REPO_NOTEBOOK = re.compile(r'(\]\(|href=")\.\./(notebooks(?:_r)?/[\w.-]+\.ipynb)(?=[)#?"])')


def on_page_markdown(markdown, page, files, **kwargs):
    # Top-level docs pages link to notebooks in the repository. sync_docs_notebooks
    # copies them inside docs/ for MkDocs. Markdown links go through MkDocs'
    # resolver; raw HTML buttons need the rendered relative URL directly.
    if not page.file.src_uri.endswith(".md") or "/" in page.file.src_uri:
        return markdown

    def replace(match):
        prefix, path = match.groups()
        if prefix == 'href="':
            target = files.get_file_from_path(path)
            if target is None:
                _LOG.warning("%s links to missing notebook %s", page.file.src_uri, path)
                return match[0]
            path = get_relative_url(target.url, page.file.url)
        return prefix + path

    return _REPO_NOTEBOOK.sub(replace, markdown)


def on_page_content(html, page, config, files):
    if not page.file.src_uri.endswith(".ipynb"):
        return html

    notebook_dir = posixpath.dirname(page.file.src_uri)

    def resolve(match):
        href = urlsplit(unescape(match[3]))
        if (href.scheme or href.netloc or href.path.startswith("/")
                or not href.path.endswith((".md", ".ipynb"))):
            return match[0]

        path = posixpath.normpath(posixpath.join(notebook_dir, unquote(href.path)))
        target = files.get_file_from_path(path)
        if target is not None:
            url = get_relative_url(target.url, page.file.url)
        elif path == _FIGURE_NOTE:
            # This illustrated note lives alongside the source notebooks.
            url = f"{config.repo_url}/blob/main/{path}"
        else:
            _LOG.warning("%s links to missing notebook or documentation source %s",
                         page.file.src_uri, href.path)
            return match[0]

        result = urlunsplit(("", "", url, href.query, href.fragment))
        return f"{match[1]}{match[2]}{escape(result, quote=True)}{match[2]}"

    return _ANCHOR.sub(resolve, html)
