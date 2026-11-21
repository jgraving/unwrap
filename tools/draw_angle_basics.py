"""Draw the vector illustration for the introductory circular-data guide."""

from math import cos, radians, sin
from pathlib import Path


OUTPUT = Path(__file__).resolve().parents[1] / "docs/assets/angle-basics.svg"
TEAL = "#08726d"
RUST = "#a54b2d"
INK = "#193131"
RULE = "#c6d5d3"


def point(cx, cy, radius, angle):
    """Project a clockwise heading with zero at twelve o'clock."""
    theta = radians(angle)
    return cx + radius * sin(theta), cy - radius * cos(theta)


def dial(cx, cy, radius, headings):
    """Mark observed headings on a circle and the mean direction as an arrow."""
    points = [point(cx, cy, radius, angle) for angle in headings]
    dots = "\n".join(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="{TEAL}"/>'
                     for x, y in points)
    arrow = cy - radius + 12
    return f'''<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" stroke="{RULE}"/>
  <path d="M {cx} {cy} V {arrow} m -4 5 l 4 -5 4 5" stroke="{RUST}" stroke-width="2"
        stroke-linecap="round" stroke-linejoin="round" fill="none"/>
  {dots}'''


def main():
    first = dial(146, 126, 60, (350, 10))
    broad = dial(400, 126, 60, (-52, -25, 0, 25, 52))
    narrow = dial(609, 126, 60, (-20, -10, 0, 10, 20))
    left, top = point(146, 126, 68, 350)
    right, _ = point(146, 126, 68, 10)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 720 221" width="720" height="221"
     role="img" aria-labelledby="title description">
  <title id="title">The circular boundary and directional concentration</title>
  <desc id="description">Headings of 350 and 10 degrees are 20 degrees apart and have a circular
    mean of zero. Two other sets of five headings have the same mean at zero, with lower and higher
    concentration. Dots show observations and arrows show circular means.</desc>
  <rect width="720" height="221" fill="white"/>
  <path d="M 280 20 V 202" fill="none" stroke="{RULE}"/>
  <g font-family="Arial, sans-serif" fill="{INK}" text-anchor="middle">
    <text x="146" y="26" font-size="16" font-weight="bold">Across the circular boundary</text>
    <text x="505" y="26" font-size="16" font-weight="bold">Same center, different concentration</text>
    {first}
    <path d="M {left:.2f} {top:.2f} A 68 68 0 0 1 {right:.2f} {top:.2f}" fill="none"
          stroke="{TEAL}" stroke-width="2"/>
    <text x="117" y="59" font-size="14">350°</text>
    <text x="176" y="59" font-size="14">10°</text>
    <text x="146" y="111" fill="{RUST}" font-size="14" text-anchor="start" dx="7">0°</text>
    <text x="146" y="210" font-size="14">20° apart</text>
    {broad}
    {narrow}
    <text x="400" y="210" font-size="14">Lower concentration</text>
    <text x="609" y="210" font-size="14">Higher concentration</text>
  </g>
</svg>
'''
    OUTPUT.write_text(svg, encoding="utf-8")


if __name__ == "__main__":
    main()
