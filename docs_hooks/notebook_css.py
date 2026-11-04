"""Remove the duplicate CSS footer appended by mkdocs-jupyter 0.27."""

import re

from mkdocs.plugins import event_priority

# The custom template already emits the real Pygments styles. Match only the
# trailing list repr so this becomes a no-op when the plugin fixes its footer.
_CSS_FOOTER = re.compile(
    r"\n {8}<style>\n {8}\[['\"]pre \{ line-height: 125%; \}[^\n]*\]\n"
    r" {8}</style>\n {8}\Z"
)


@event_priority(0)  # mkdocs-jupyter supplies the notebook fragment at priority 100.
def on_page_markdown(markdown, page, **kwargs):
    if page.file.src_uri.endswith(".ipynb"):
        return _CSS_FOOTER.sub("", markdown)
