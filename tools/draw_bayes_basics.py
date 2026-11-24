"""Render the circular Bayesian basics diagrams using only the standard library."""

from math import atan2, cos, degrees, exp, hypot, radians, sin, tau
from pathlib import Path
import xml.etree.ElementTree as ET


ASSETS = Path(__file__).resolve().parents[1] / "docs" / "assets"
OBSERVATIONS = (35, 55, 75, 95)
PRIOR_CENTER = 0
PRIOR_CONCENTRATION = 2
HEADING_CONCENTRATION = 2
GRID = tuple(tau * i / 360 for i in range(360))
STEP = tau / len(GRID)
PRIOR = "#126A62"
LIKELIHOOD = "#94563C"
POSTERIOR = "#413970"
INK = "#192725"
MUTED = "#43534F"
GUIDE = "#A7B8B4"


def normalized(values):
    """Normalize samples of a periodic curve by trapezoidal quadrature."""
    values = tuple(values)
    integral = STEP * sum(values)
    return tuple(value / integral for value in values)


def von_mises(center, concentration):
    return normalized(exp(concentration * (cos(angle - center) - 1)) for angle in GRID)


def curves():
    prior_x = PRIOR_CONCENTRATION * cos(radians(PRIOR_CENTER))
    prior_y = PRIOR_CONCENTRATION * sin(radians(PRIOR_CENTER))
    data_x = HEADING_CONCENTRATION * sum(cos(radians(y)) for y in OBSERVATIONS)
    data_y = HEADING_CONCENTRATION * sum(sin(radians(y)) for y in OBSERVATIONS)
    data_center = atan2(data_y, data_x)
    posterior_center = atan2(prior_y + data_y, prior_x + data_x)

    prior = von_mises(radians(PRIOR_CENTER), PRIOR_CONCENTRATION)
    likelihood = von_mises(data_center, hypot(data_x, data_y))
    posterior = von_mises(
        posterior_center, hypot(prior_x + data_x, prior_y + data_y)
    )
    heading_density = von_mises(0, HEADING_CONCENTRATION)
    n = len(GRID)
    predictive = tuple(
        STEP * sum(weight * heading_density[(i - j) % n] for j, weight in enumerate(posterior))
        for i in range(n)
    )
    return prior, likelihood, posterior, predictive, degrees(data_center), degrees(posterior_center)


def label(svg, x, y, value, size, color=INK, weight="normal"):
    ET.SubElement(svg, "text", {
        "x": str(x), "y": str(y), "fill": color, "font-size": str(size),
        "font-weight": weight, "text-anchor": "middle",
    }).text = value


def polar(svg, x, title, note, values, color, scale, observations=()):
    y, reference = 120, 51

    def point(angle, radius):
        return x + radius * sin(angle), y - radius * cos(angle)

    label(svg, x, 20, title, 17, weight="bold")
    label(svg, x, 52, "0°", 11, MUTED)
    label(svg, x, 191, "180°", 11, MUTED)
    ET.SubElement(svg, "circle", {
        "cx": str(x), "cy": str(y), "r": str(reference),
        "fill": "none", "stroke": GUIDE, "stroke-width": "1",
    })
    ticks = []
    for position in (0, 90, 180, 270):
        angle = radians(position)
        inner = point(angle, reference)
        outer = point(angle, reference + 4)
        ticks.append(f"M{inner[0]:.2f},{inner[1]:.2f}L{outer[0]:.2f},{outer[1]:.2f}")
    ET.SubElement(svg, "path", {"d": "".join(ticks), "stroke": GUIDE, "stroke-width": "1"})

    path = "L".join(
        f"{px:.2f},{py:.2f}"
        for px, py in (point(angle, reference + scale * value)
                       for angle, value in zip(GRID, values, strict=True))
    )
    ET.SubElement(svg, "path", {
        "d": f"M{path}Z", "fill": color, "fill-opacity": ".09",
        "stroke": color, "stroke-width": "2", "stroke-linejoin": "round",
    })
    # Cover the interior of the density trace so its fill only shows density above zero.
    ET.SubElement(svg, "circle", {
        "cx": str(x), "cy": str(y), "r": str(reference - .5),
        "fill": "white",
    })
    for obs in observations:
        px, py = point(radians(obs), reference - 9)
        ET.SubElement(svg, "circle", {
            "cx": f"{px:.2f}", "cy": f"{py:.2f}", "r": "3",
            "fill": color,
        })
    label(svg, x, 211, note, 13, MUTED)


def make_svg(filename, title, description, panels, scale):
    width = 222 * len(panels)
    svg = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg", "viewBox": f"0 0 {width} 227",
        "width": str(width), "height": "227",
        "role": "img", "aria-labelledby": "title description",
    })
    ET.SubElement(svg, "title", {"id": "title"}).text = title
    ET.SubElement(svg, "desc", {"id": "description"}).text = description
    ET.SubElement(svg, "rect", {"width": str(width), "height": "227", "fill": "white"})
    drawing = ET.SubElement(svg, "g", {"font-family": "Arial, sans-serif"})
    for i, (panel_title, note, values, color, *observations) in enumerate(panels):
        polar(drawing, i * 222 + 111, panel_title, note, values, color, scale, *observations)
    ET.indent(svg)
    ET.ElementTree(svg).write(ASSETS / filename, encoding="unicode")


def main():
    prior, likelihood, posterior, predictive, data_center, posterior_center = curves()
    scale = 32 / max(map(max, (prior, likelihood, posterior, predictive)))
    fitted_note = f"Mean direction near {posterior_center:.0f}°"
    make_svg(
        "bayes-heading-update.svg",
        "A Bayesian update for the preferred heading",
        f"The prior favors 0 degrees. Four observations favor {data_center:.0f} degrees. "
        f"The posterior is centered near {posterior_center:.0f} degrees and is tighter than the prior.",
        (
            ("Prior", "Prior expectation at 0°", prior, PRIOR),
            ("Likelihood", f"Data favor {data_center:.0f}°", likelihood, LIKELIHOOD, OBSERVATIONS),
            ("Posterior", fitted_note, posterior, POSTERIOR),
        ),
        scale,
    )
    make_svg(
        "bayes-heading-prediction.svg",
        "Preferred heading compared with a new observed heading",
        "The posterior for the mean direction and the posterior predictive density for a new heading "
        "are centered at the same angle. The predictive curve is broader.",
        (
            ("Mean direction", "Uncertainty about the center", posterior, POSTERIOR),
            ("New heading", "Variation in a new observation", predictive, PRIOR),
        ),
        scale,
    )


if __name__ == "__main__":
    main()
