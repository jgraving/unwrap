"""
Shared matplotlib plumbing for the research notebooks.

Every research notebook draws figures with the same look (Arial, white
figures, type-42 PDF outlines, consistent ticks) and saves them to a
per-study subdirectory of ``<repo>/figures/``. The per-study
``notebooks/<study>_plots.py`` modules hold the study-specific figure
functions; this module holds the machinery those functions share so the
style and save behaviour live in exactly one place.

The save helpers report paths *relative to the repository root* (e.g.
``figures/bats/kappa_hat.pdf``) so baked notebook output is identical no
matter where the checkout lives.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl


#: Defaults shared by every study figure set.
DEFAULT_RC = {
    "font.family": "Arial",
    "font.sans-serif": ["Arial"],
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def set_style(figure_dir, extra_rc=None):
    """Apply the shared figure style and echo where figures will be saved.

    ``figure_dir`` is the absolute path of the study's figures directory (a
    subdirectory of ``<repo>/figures/``). ``extra_rc`` merges study-specific
    rcParams (font size, linewidths, tick sizes) on top of :data:`DEFAULT_RC`.
    """
    rc = dict(DEFAULT_RC)
    rc.update(extra_rc or {})
    mpl.rcParams.update(rc)
    print(f"Figures will be saved to: {repo_rel(Path(figure_dir))}")
    return rc


def save_figure(fig, figure_dir, name, pdf=True, png=False, transparent=False,
                dpi=300, pad_inches=0.01):
    """Save ``fig`` into ``figure_dir`` as ``<name>.pdf`` and/or ``<name>.png``.

    Creates ``figure_dir`` if needed and prints the repo-relative paths
    written, matching the text the notebooks' figure cells always showed.
    Returns the list of paths written (absolute).
    """
    out = Path(figure_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    written = []
    if pdf:
        p = out / (name + ".pdf")
        fig.savefig(p, dpi=dpi, bbox_inches="tight", pad_inches=pad_inches,
                    transparent=transparent)
        written.append(p)
    if png:
        p = out / (name + ".png")
        fig.savefig(p, dpi=dpi)
        written.append(p)
    if written:
        shown = [str(repo_rel(p)) for p in written]
        plural = "s" if len(shown) > 1 else ""
        print("Saved figure" + plural + " to: " + "; ".join(shown))
    return written


def repo_rel(path):
    """``path`` relative to the repository root, as a plain string.

    The figures trees live at ``<root>/figures/<study>/``; for a file under
    that the result is ``figures/<study>/<file>``, and for the study
    directory itself ``figures/<study>``. Paths with no ``figures`` ancestor
    are returned unchanged as strings.
    """
    p = Path(path).resolve()
    root = None
    for i, anc in enumerate(p.parents):
        if anc.name == "figures" and i + 1 < len(p.parents):
            root = p.parents[i + 1]
            break
    if root is None:
        return str(p)
    try:
        return str(p.relative_to(root))
    except ValueError:
        return str(p)


def tick_style(ax, width=1.0, length=None):
    """Point every tick outward with a consistent width (and optional length).

    Replaces the per-notebook ``_restore_ticks`` helpers; safe to call on any
    axis to normalise matplotlib's default inward ticks.
    """
    kwargs = dict(bottom=True, top=False, left=True, right=False, width=width)
    if length is not None:
        kwargs["length"] = length
    ax.tick_params(axis="both", which="both", direction="out", **kwargs)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
