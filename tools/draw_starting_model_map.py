"""Render the linked, full-width model map shared by the paper and docs."""
from html import escape
from math import cos, exp, pi, radians, sin
from pathlib import Path


DOCS = "https://jakegraving.com/unwrap/"
INK, MUTED, LINK = "#1f2827", "#455655", "#13665e"
RULE, RUST = "#c6d3d0", "#9e6646"
EDGE = "#738782"
OUTPUT = Path(__file__).resolve().parents[1] / "docs/assets/unwrap-model-map.svg"


def text(x, y, label, *, size=13, color=INK, bold=False, anchor="start"):
    weight = ' font-weight="bold"' if bold else ""
    return (f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" fill="{color}"{weight}>'
            f'{escape(label)}</text>')


def link(key, x, y, label, path, *, size=12.5, anchor="start"):
    url = escape(DOCS + path, quote=True)
    return (f'<a id="{key}" href="{url}" xlink:href="{url}">'
            f'{text(x, y, label, size=size, color=LINK, bold=True, anchor=anchor)}</a>')


def rule(y):
    return f'<path d="M22 {y}H704" stroke="{RULE}" fill="none"/>'


def profile(components, offset, color, dash=""):
    points = []
    for index in range(121):
        angle = index * pi / 60
        radius = offset + sum(amplitude * exp(concentration * (cos(angle - radians(mean)) - 1))
                              for mean, concentration, amplitude in components)
        points.append(f"{radius * sin(angle):.2f},{-radius * cos(angle):.2f}")
    style = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<path d="M{" L".join(points)}Z" fill="none" stroke="{color}" '
            f'stroke-width="1.35" stroke-linejoin="round"{style}/>')


def schematic(x, y, model, scale=1):
    """Sketch conditional, group-specific, or mixture densities on a common circle."""
    profiles = {
        "glm": (
            profile([(-47, 4.5, 14)], 7, LINK),
            profile([(42, 3, 11)], 7, RUST, "2.4 1.5"),
        ),
        "glmm": (
            profile([(-34, 4, 12)], 7, "#7c9d97", "2 1.5"),
            profile([(34, 4, 12)], 7, "#7c9d97", "2 1.5"),
            profile([(0, 2.7, 15)], 7, LINK),
        ),
        "mixture": (
            '<circle r="7" fill="none" stroke="#9e6646" stroke-width="1.2" stroke-dasharray="2 1.5"/>',
            profile([(-11, 7, 12)], 9, LINK),
        ),
        "bimodal": (
            profile([(-64, 6, 14), (64, 5, 12)], 7, LINK),
        ),
    }
    return (f'<g transform="translate({x} {y}) scale({scale})" aria-hidden="true">'
            f'<circle r="23" fill="none" stroke="{RULE}"/>'
            f'<path d="M0 -26v3m0 46v3m-26 -26h3m46 0h3" stroke="{RULE}" fill="none"/>'
            f'{"".join(profiles[model])}</g>')


def arrow(path):
    return (f'<path d="{path}" stroke="{EDGE}" stroke-width="1" fill="none" '
            'marker-end="url(#arrow)"/>')


def example(x, y, label, notebook):
    return [
        text(x - 97, y, label, size=12.5, bold=True),
        link(f"{notebook}-python", x + 38, y, "Python", f"notebooks/{notebook}/", size=11.8),
        link(f"{notebook}-r", x + 91, y, "R", f"notebooks_r/{notebook}/", size=11.8),
    ]


def render():
    topics = (
        (130, "glm", "Headings under different conditions",
         ("Does a condition shift direction", "or change directional consistency?"),
         "Circular GLM", "math/glms/",
         (("Regression tutorial", "circular_regression"),)),
        (363, "glmm", "Repeated headings",
         ("Do animals or groups differ in their", "preferred headings or cue responses?"),
         "Circular GLMM", "math/glmms/",
         (("Beetle orientation", "beetles"),)),
        (596, "mixture", "Several possible heading patterns",
         ("Could headings reflect directed and", "undirected travel, or competing bearings?"),
         "Circular mixture", "math/mixtures/",
         (("Bat homing", "goldshtein_bats"), ("Locust steering", "locust_direction_mixture_model"))),
    )
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        'width="178mm" height="89.49mm" viewBox="0 0 726 365" '
        'aria-labelledby="figure-title figure-desc">',
        '<title id="figure-title">Getting started with unwrap</title>',
        '<desc id="figure-desc">Panel a links three schematic heading patterns to model guides '
        'and worked examples. Animal or group effects can extend both a single circular regression '
        'and a mixture. Panel b links to foundations and a starter prompt for coding agents.</desc>',
        '<defs><marker id="arrow" viewBox="0 0 7 7" refX="6.5" refY="3.5" markerWidth="6" '
        'markerHeight="6" orient="auto" markerUnits="userSpaceOnUse">'
        f'<path d="M1 1 6.5 3.5 1 6" stroke="{EDGE}" stroke-width="1" fill="none"/>'
        '</marker></defs>',
        '<rect width="726" height="365" fill="white"/>',
        '<g font-family="Arial, Liberation Sans, sans-serif">',
        text(22, 23, "a", size=18, bold=True),
        text(41, 23, "Choose a circular model", size=14, bold=True),
        rule(147),
    ]
    for x, family, data, question, name, guide, examples in topics:
        if family == "mixture":
            parts.extend([schematic(x - 31, 80, family, .89), schematic(x + 31, 80, "bimodal", .89)])
        else:
            parts.append(schematic(x, 80, family))
        parts.extend([
            text(x, 46, data, size=12.8, bold=True, anchor="middle"),
            text(x, 120, question[0], size=12.2, color=MUTED, anchor="middle"),
            text(x, 135, question[1], size=12.2, color=MUTED, anchor="middle"),
            arrow(f"M{x} 150v17"),
            f'<rect x="{x-64}" y="169" width="128" height="25" rx="3" '
            f'fill="#f5faf8" stroke="{LINK}" stroke-width=".85"/>',
            link(f"model-{family}", x, 187, name, guide, size=13.6, anchor="middle"),
            arrow(f"M{x} 195v14"),
        ])
        for line, (title, notebook) in enumerate(examples):
            parts.extend(example(x, 226 + line * 17, title, notebook))
    parts.extend([
        arrow("M195 183H297"),
        text(246, 176, "Add group effects", size=11.55, color=MUTED, anchor="middle"),
        text(596, 261, "Mixtures can also include group effects.", size=11.6, color=MUTED, anchor="middle"),
        rule(271),
        text(22, 292, "b", size=18, bold=True),
        text(41, 292, "Choose an entry point", size=14, bold=True),
        text(22, 314, "New to statistical modeling", size=12.25, bold=True),
        link("entry-biology", 22, 332, "Start with a biological question", "#start-with-a-biological-question", size=11.7),
        text(219, 314, "New to circular data", size=12.25, bold=True),
        link("entry-angles", 219, 332, "Angles and concentration", "math/von-mises/#angular-conventions", size=11.7),
        text(397, 314, "New to Bayesian GLMs", size=12.25, bold=True),
        link("entry-glms", 397, 332, "GLMs", "models/#from-a-biological-question-to-a-glm", size=11.7),
        text(432, 332, "·", size=11.7, color=MUTED),
        link("entry-bayes", 443, 332, "Bayesian basics", "#read-a-fitted-model", size=11.7),
        text(569, 314, "Using a coding agent", size=12.25, bold=True),
        link("entry-agent", 569, 332, "Copyable prompt", "coding-agents/#copyable-prompt", size=11.7),
        rule(342),
        link("docs", 22, 352, "Package documentation", "", size=11.6),
        link("setup-python", 577, 352, "Python setup", "getting-started/", size=11.6),
        text(653, 352, "·", size=11.6, color=MUTED),
        link("setup-r", 665, 352, "R setup", "getting-started-r/", size=11.6),
        '</g>',
        '</svg>',
    ])
    return "\n".join(parts) + "\n"


if __name__ == "__main__":
    OUTPUT.write_text(render(), encoding="utf-8")
