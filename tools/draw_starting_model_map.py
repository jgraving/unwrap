"""Render the linked, full-width model map shared by the paper and docs."""
from html import escape
from math import cos, exp, pi, radians, sin
from pathlib import Path


DOCS = "https://jakegraving.com/unwrap/"
INK, MUTED, LINK = "#1f2827", "#455655", "#13665e"
RULE, RUST = "#c6d3d0", "#9e6646"
OUTPUT = Path(__file__).resolve().parents[1] / "docs/assets/unwrap-model-map.svg"


def text(x, y, label, *, size=13, color=INK, bold=False):
    weight = ' font-weight="bold"' if bold else ""
    return (f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}"{weight}>'
            f'{escape(label)}</text>')


def link(key, x, y, label, path, *, size=12.5):
    url = escape(DOCS + path, quote=True)
    return (f'<a id="{key}" href="{url}" xlink:href="{url}">'
            f'{text(x, y, label, size=size, color=LINK, bold=True)}</a>')


def rule(y):
    return f'<path d="M22 {y}H704" stroke="{RULE}" fill="none"/>'


def profile(mean, concentration, amplitude, offset, color, dash=""):
    points = []
    for index in range(121):
        angle = index * pi / 60
        radius = offset + amplitude * exp(concentration * (cos(angle - radians(mean)) - 1))
        points.append(f"{radius * sin(angle):.2f},{-radius * cos(angle):.2f}")
    style = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<path d="M{" L".join(points)}Z" fill="none" stroke="{color}" '
            f'stroke-width="1.35" stroke-linejoin="round"{style}/>')


def schematic(y, model):
    """Sketch conditional, group-specific, or mixture densities on a common circle."""
    profiles = {
        "glm": (
            profile(-47, 4.5, 14, 7, LINK),
            profile(42, 3, 11, 7, RUST, "2.4 1.5"),
        ),
        "glmm": (
            profile(-34, 4, 12, 7, "#7c9d97", "2 1.5"),
            profile(34, 4, 12, 7, "#7c9d97", "2 1.5"),
            profile(0, 2.7, 15, 7, LINK),
        ),
        "mixture": (
            '<circle r="7" fill="none" stroke="#9e6646" stroke-width="1.2" stroke-dasharray="2 1.5"/>',
            profile(-11, 7, 12, 9, LINK),
        ),
    }
    return (f'<g transform="translate(49 {y})" aria-hidden="true">'
            f'<circle r="23" fill="none" stroke="{RULE}"/>'
            f'<path d="M0 -26v3m0 46v3m-26 -26h3m46 0h3" stroke="{RULE}" fill="none"/>'
            f'{"".join(profiles[model])}</g>')


def example(y, label, notebook):
    return [
        text(507, y, label, size=12.5, bold=True),
        link(f"{notebook}-python", 637, y, "Python", f"notebooks/{notebook}/", size=11.8),
        link(f"{notebook}-r", 693, y, "R", f"notebooks_r/{notebook}/", size=11.8),
    ]


def render():
    rows = (
        ("glm", "Headings under different conditions",
         ("Do conditions change where animals point", "or how consistently they point?"),
         ("Circular GLM", "math/glms/", ("Mean direction and", "concentration regression")),
         (("Regression tutorial", "circular_regression"),)),
        ("glmm", "Repeated headings within groups",
         ("Do animals or groups differ in their", "preferred headings or cue responses?"),
         ("Circular GLMM", "math/glmms/", ("Partially pooled animal", "or group effects")),
         (("Beetle orientation", "beetles"),)),
        ("mixture", "Several possible heading patterns",
         ("Could headings reflect directed and", "undirected travel, or competing bearings?"),
         ("Circular mixture", "math/mixtures/", ("Multiple components", "with estimated weights")),
         (("Bat homing", "goldshtein_bats"), ("Locust steering", "locust_direction_mixture_model"))),
    )
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        'width="178mm" height="88.15mm" viewBox="0 0 726 360" '
        'aria-labelledby="figure-title figure-desc">',
        '<title id="figure-title">Getting started with unwrap</title>',
        '<desc id="figure-desc">Panel a pairs angular data and research questions with starting '
        'models and worked examples. Panel b links to background for readers new to circular data, '
        'Bayesian regression, or statistical modeling. All densities are schematic.</desc>',
        '<rect width="726" height="360" fill="white"/>',
        '<g font-family="Arial, Liberation Sans, sans-serif">',
        text(22, 23, "a", size=18, bold=True),
        text(41, 23, "Choose a model for your data and question", size=14, bold=True),
        text(83, 46, "Data and question", size=12.2, color=MUTED, bold=True),
        text(342, 46, "Starting model", size=12.2, color=MUTED, bold=True),
        text(507, 46, "Worked examples", size=12.2, color=MUTED, bold=True),
        rule(53),
    ]
    for index, (family, data, question, model, examples) in enumerate(rows):
        top = 53 + index * 62
        name, guide, explanation = model
        parts.extend([
            schematic(top + 31, family),
            text(83, top + 18, data, size=12.5, bold=True),
            text(83, top + 36, question[0], size=12, color=MUTED),
            text(83, top + 51, question[1], size=12, color=MUTED),
            link(f"model-{family}", 342, top + 18, name, guide, size=13.6),
            text(342, top + 36, explanation[0], size=12, color=MUTED),
            text(342, top + 51, explanation[1], size=12, color=MUTED),
        ])
        for line, (title, notebook) in enumerate(examples):
            y = top + (33 if len(examples) == 1 else 22 + line * 22)
            parts.extend(example(y, title, notebook))
        parts.append(rule(top + 62))
    parts.extend([
        text(83, 256, "Mixture models can also include animal or group effects.", size=11.6, color=MUTED),
        text(22, 283, "b", size=18, bold=True),
        text(41, 283, "Choose an entry point", size=14, bold=True),
        text(22, 305, "I know regression and GLMs", size=12.3, bold=True),
        link("entry-angles", 22, 322, "Angles and concentration", "math/von-mises/#angular-conventions", size=12),
        text(253, 305, "I know circular statistics", size=12.3, bold=True),
        link("entry-glms", 253, 322, "Circular GLMs", "models/#from-a-biological-question-to-a-glm", size=12),
        text(341, 322, "·", size=12, color=MUTED),
        link("entry-bayes", 352, 322, "Bayesian basics", "#read-a-fitted-model", size=12),
        text(507, 305, "I’m new to statistical modeling", size=12.3, bold=True),
        link("entry-biology", 507, 322, "Start with a biological question", "#start-with-a-biological-question", size=12),
        rule(332),
        link("docs", 22, 350, "jakegraving.com/unwrap", "", size=11.7),
        link("setup-python", 577, 350, "Python setup", "getting-started/", size=11.7),
        text(653, 350, "·", size=12, color=MUTED),
        link("setup-r", 665, 350, "R setup", "getting-started-r/", size=11.7),
        '</g>',
        '</svg>',
    ])
    return "\n".join(parts) + "\n"


if __name__ == "__main__":
    OUTPUT.write_text(render(), encoding="utf-8")
