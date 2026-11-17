"""Remove the broken CSS footer appended by mkdocs-jupyter 0.27."""

import re

from mkdocs.plugins import event_priority

# mkdocs-jupyter 0.27's nbconvert pass appends a <style> block at the end
# of each notebook page whose content is not CSS but the Python repr of
# the page's list of styles, e.g. "['pre { line-height: 125%; }', ...]".
# The block is dead markup: browsers ignore it, and the highlight styles
# the page actually needs already come from the custom template. Under
# markdown_cell_processor: mkdocs the same footer goes through the
# markdown pipeline instead and renders as a visible code block.

# The strip has to run on the rendered page: with
# markdown_cell_processor: jupyter the plugin renders notebook pages
# through nbconvert's HTML (it replaces page.render) and never passes
# them through on_page_markdown, which is why an earlier
# on_page_markdown strip was inert. Priority 0 runs after the
# plugin's own default-priority on_page_content. The rendered page has
# scripts appended after the footer, so the match cannot anchor at end
# of file; the list-repr content is the signature. The match is at most
# one per page, so stripping is idempotent.
_CSS_FOOTER = re.compile(
    r"\n *<style>\n *\[['\"]pre \{ line-height: 125%; \}[\s\S]*?\]\n"
    r" *</style>\n ?"
)


@event_priority(0)  # Late enough for the plugin's on_page_content to run first.
def on_page_content(html, page, **kwargs):
    if page.file.src_uri.endswith(".ipynb"):
        return _CSS_FOOTER.sub("", html)
