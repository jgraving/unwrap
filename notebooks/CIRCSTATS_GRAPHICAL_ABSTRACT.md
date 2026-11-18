# Reproduce the circstats graphical abstract

`notebooks/circstats_graphical_abstract.py` rebuilds the graphical
abstract published with the circstats (circular statistics) elife paper.
The committed `notebooks/assets/circstats_graphical_abstract/
reference.pdf` (599 × 342 points, 100 % vector) is the
Inkscape-polished target, `figures/Graphical Abstract.pdf` from
`jgraving/circstats`; every coordinate, colour, stroke width, dash
pattern and string literal in the generator comes from a measurement of
that file, pinned in `artwork.json` and SHA-256 verified by
`tests/test_circstats_graphical_abstract.py`.

The figure is a single wide strip with four panels:

- **a** — the experimental design and data collection: the two foraging
  trajectories (control, blue; treatment, green) and the 3 × 4 grid of
  hand-drawn animal circles, each ring annotated with the group's
  observed angle points.
- **b** — model specification and parameter inference: the syntax-
  highlighted `model.py` code window, the η-μ and η-κ MCMC trace panels,
  the circular modulo and softplus link curves, and the posterior
  densities with the η vector glyphs.
- **c** — model diagnostics, evaluation, and selection: the posterior
  predictive density cluster over the observed angle dots, the
  population- and individual-level KDE summary circles with their arrow
  fans, and the three-model comparison (point, interval, band).
- **d** — model summaries, scientific inference, and interpretation:
  the predictor curves with the uniform reference, the ellipsoidal
  angle loop, and the two hypothesis interval estimates against the
  uniform prior.

Between the panels, four flow arrows and a legend (data summary /
predictive density) explain the workflow direction.

## How the figure is reconstructed

PyMuPDF's `page.get_drawings()` and `page.get_text("dict")` were used to
measure the reference PDF (page coordinates in points, y down from the
top-left corner). The 526 drawings split into two classes, recorded in
the asset:

- **472 pinned drawings** — hand-made or stochastic artwork that cannot
  be re-simulated: the 12 hand-drawn animal circle units, the foraging
  trajectory (68 strokes), the code editor window (8 shapes), the MCMC
  trace wiggles (15), the posterior-cluster curves (112), the jagged
  posterior density fills (7), the observed-angle dots (8), the panel-d
  predictor curves (7), the ellipsoidal loop (12) and the KDE blob/fan
  interiors (61). These are replayed verbatim.
- **54 skeleton drawings** — the deterministic structure: axis spines
  and ticks, the modulo and softplus curves, the dashed uniform
  reference lines, the model/hypothesis interval bars, the flow arrows,
  the legend frame and swatches, the four KDE circle outlines and the η
  glyph marks. These are drawn in code from the same measured geometry.

All 236 text spans (Arial labels, the Menlo syntax-highlighted code
window, the rotated axis captions) are typeset live at their measured
baselines — each span carries its embedded font's ascender, so the
baseline sits at `y + ascender·size` independent of the substitute face
— which keeps the exports editable (`svg.fonttype = "none"`).  The two
STIX-math spans are whitespace placeholders (the η labels are the vector
glyphs above) and are skipped.

Nothing is randomized: every value is read from the hash-pinned asset,
so two renders in the same environment produce byte-identical exports.
CI renders twice with Liberation Sans (metric-compatible with the
reference Arial) and `sha256sum --check`s the SVG and PDF
(`.github/workflows/circstats-graphical-abstract.yml`).

## Generate

From the repository root, with Python 3.12+:

```bash
python -m pip install -e '.[figure]'
python notebooks/circstats_graphical_abstract.py
```

No downloads: the only inputs are the two committed asset files.
`--font-family 'Liberation Sans'` selects the CI face; the provenance
JSON records what actually resolved. The outputs default to
`assets/circstats-graphical-abstract/` (committed for review):

| Output | Content |
| --- | --- |
| `circstats_graphical_abstract.svg` | 599 × 342 pt vector figure, 149 editable text spans, no rasters |
| `circstats_graphical_abstract.pdf` | 599 × 342 pt; CreationDate/ModDate cleared and /ID pinned to the content hash |
| `circstats_graphical_abstract.png` | 1797 × 1026 px preview rendered from the PDF |
| `circstats_graphical_abstract.html` | Self-contained page with the caption |
| `circstats_graphical_abstract.json` | Provenance: input hashes, resolved fonts, counts, env, export hashes |

## Fidelity notes

- The 472 pinned drawings are byte-exact replays of the published
  vector paths, so the hand-made and stochastic content — including
  every manual Inkscape adjustment — reproduces exactly.
- Typography: the reference's Arial, Arial Bold and Helvetica render
  with Liberation Sans in CI (metric-compatible; the provenance JSON
  records the actual face) and Menlo falls back to DejaVu Sans Mono, so
  long titles may carry a fraction-of-a-point per-glyph advance width
  difference there. With the real fonts installed, all lettering
  baselines match the reference span for span.
- The PDF export is stabilised: wall-clock dates are cleared and the
  trailer `/ID` is set to the document's SHA-256, keeping committed
  refreshes diff-clean.
