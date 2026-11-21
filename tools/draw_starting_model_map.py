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


def card(x, y, width, height, *, emphasis=False):
    """Draw a card with softly squared corners and a thin outline."""
    right, bottom = x + width, y + height
    radius = 10 if height >= 45 else 8
    fill, stroke = ("#f3f9f6", "#83a5a0") if emphasis else ("white", RULE)
    d = (
        f"M{x + radius} {y}H{right - radius}"
        f"Q{right} {y} {right} {y + radius}V{bottom - radius}"
        f"Q{right} {bottom} {right - radius} {bottom}H{x + radius}"
        f"Q{x} {bottom} {x} {bottom - radius}V{y + radius}"
        f"Q{x} {y} {x + radius} {y}Z"
    )
    return f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="1"/>'


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


def connector(path, *, arrow=True):
    marker = ' marker-end="url(#arrow)"' if arrow else ""
    return (f'<path d="{path}" stroke="{EDGE}" stroke-width="1" stroke-linecap="round" '
            f'stroke-linejoin="round" fill="none"{marker}/>')


def density_schematic(x):
    """Illustrate a continuous, saturating relationship without a data scale."""
    left, right, bottom, top = x - 43, x + 43, 234, 195
    points = [
        f"{left + 4 + (right - left - 8) * t / 60:.2f},{230 - 29 * (1 - exp(-3.6 * t / 60)):.2f}"
        for t in range(61)
    ]
    return [
        f'<path d="M{left} {top}V{bottom}H{right}" stroke="{RULE}" stroke-width="1.2" fill="none"/>',
        f'<path d="M{" L".join(points)}" stroke="{LINK}" stroke-width="1.8" stroke-linecap="round" fill="none"/>',
    ]


def notebook_links(x, y, notebook):
    return [
        link(f"{notebook}-python", x, y, "Python", f"notebooks/{notebook}/", size=11.8),
        text(x + 43, y, "·", size=11.8, color=MUTED),
        link(f"{notebook}-r", x + 52, y, "R", f"notebooks_r/{notebook}/", size=11.8),
    ]


def render():
    topics = (
        (130, "glm", "Treatment groups",
         ("Does direction or consistency", "differ between conditions?"),
         (("Beetle light pollution", "beetles"), ("Zebra escape", "zebra_escape"))),
        (363, "continuous", "Continuous predictors",
         ("How does direction or consistency", "change along a gradient?"),
         (("Locust heading density", "locust_density_order"),)),
        (596, "mixture", "Mixtures",
         ("Directed and undirected travel,", "or two preferred directions?"),
         (("Bat homing", "goldshtein_bats"), ("Locust decision-making", "locust_direction_mixture_model"))),
    )
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        'width="178mm" height="97.34mm" viewBox="0 0 726 397" '
        'aria-labelledby="figure-title figure-desc">',
        '<title id="figure-title">Getting started with unwrap</title>',
        '<desc id="figure-desc">The package guide branches to foundations, setup instructions, '
        'and a coding-agent prompt. Treatment groups and continuous predictors lead to circular GLMMs. '
        'A schematic curve shows increasing predictive density at 0° with swarm order. Mixtures lead to the '
        'circular mixture GLMM. The final three boxes link to worked examples.</desc>',
        '<defs><marker id="arrow" viewBox="0 0 7 7" refX="6.5" refY="3.5" markerWidth="6" '
        'markerHeight="6" orient="auto" markerUnits="userSpaceOnUse">'
        f'<path d="M1 1 6.5 3.5 1 6" stroke="{EDGE}" stroke-width="1" fill="none"/>'
        '</marker></defs>',
        '<rect width="726" height="397" fill="white"/>',
        '<g font-family="Arial, Liberation Sans, sans-serif">',
        connector("M363 34V41M130 41H596", arrow=False),
        connector("M130 41V48"),
        connector("M363 41V48"),
        connector("M596 41V48"),
        connector("M130 104V115H596V104M363 104V115", arrow=False),
        connector("M363 115V126"),
        connector("M363 150V157M130 157H596", arrow=False),
        connector("M130 157V166"),
        connector("M363 157V166"),
        connector("M596 157V166"),
        card(267, 9, 192, 25, emphasis=True),
        link("docs", 363, 26, "Start with the package guide", "", size=13.2, anchor="middle"),
        card(22, 49, 216, 55),
        card(255, 49, 216, 55),
        card(488, 49, 216, 55),
        text(130, 65, "Learn the foundations", size=12.5, bold=True, anchor="middle"),
        link("entry-biology", 130, 81, "Biological questions", "#start-with-a-biological-question", size=11.65, anchor="middle"),
        link("entry-angles", 41, 96, "Angles", "angle-basics/", size=11.65),
        text(85, 96, "·", size=11.65, color=MUTED),
        link("entry-glms", 95, 96, "GLMs", "models/#from-a-biological-question-to-a-glm", size=11.65),
        text(131, 96, "·", size=11.65, color=MUTED),
        link("entry-bayes", 141, 96, "Bayesian basics", "bayes-basics/", size=11.65),
        text(363, 65, "Set up to run a notebook", size=12.5, bold=True, anchor="middle"),
        link("setup-python", 303, 81, "Python setup", "getting-started/", size=11.65),
        text(380, 81, "·", size=11.65, color=MUTED),
        link("setup-r", 390, 81, "R setup", "getting-started-r/", size=11.65),
        link("intro-tutorial", 363, 96, "Try a simulated regression", "examples/#the-tutorial", size=11.65, anchor="middle"),
        text(596, 65, "Use a coding agent", size=12.5, bold=True, anchor="middle"),
        link("entry-agent", 596, 81, "Copyable prompt", "coding-agents/#copyable-prompt", size=11.65, anchor="middle"),
        text(596, 96, "Bring your question and data", size=11.65, color=MUTED, anchor="middle"),
        card(248, 127, 230, 23, emphasis=True),
        text(363, 143, "Match your data and question", size=12.8, bold=True, anchor="middle"),
        connector("M130 268V281H363V268", arrow=False),
        connector("M246.5 281V295"),
        connector("M596 268V295"),
        connector("M246.5 322V331M130 331H363", arrow=False),
        connector("M130 331V340"),
        connector("M363 331V340"),
        connector("M596 322V340"),
    ]
    for x, family, data, question, examples in topics:
        parts.append(card(x - 108, 167, 216, 101))
        if family == "mixture":
            parts.extend([schematic(x - 31, 215, family, .79), schematic(x + 31, 215, "bimodal", .79)])
        elif family == "continuous":
            parts.extend(density_schematic(x))
        else:
            parts.append(schematic(x, 215, family, .84))
        parts.extend([
            text(x, 185, data, size=13.2, bold=True, anchor="middle"),
            text(x, 246, question[0], size=11.8, color=MUTED, anchor="middle"),
            text(x, 261, question[1], size=11.8, color=MUTED, anchor="middle"),
            card(x - 108, 341, 216, 46),
        ])
        if len(examples) == 1:
            (title, notebook), = examples
            parts.append(text(x, 358, title, size=11.9, bold=True, anchor="middle"))
            parts.extend(notebook_links(x - 30, 376, notebook))
        else:
            for line, (title, notebook) in enumerate(examples):
                y = 359 + line * 17
                parts.append(text(x - 98, y, title, size=11.65, bold=True))
                parts.extend(notebook_links(x + 42, y, notebook))
    parts.extend([
        card(155.5, 296, 182, 26, emphasis=True),
        link("model-glmm", 246.5, 314, "Circular GLMM", "math/glmms/", size=13.5, anchor="middle"),
        card(488, 296, 216, 26, emphasis=True),
        link("model-mixture", 596, 314, "Circular mixture GLMM", "math/mixtures/#add-animal-differences-to-a-mixture", size=13.5, anchor="middle"),
        '</g>',
        '</svg>',
    ])
    return "\n".join(parts) + "\n"


if __name__ == "__main__":
    OUTPUT.write_text(render(), encoding="utf-8")
