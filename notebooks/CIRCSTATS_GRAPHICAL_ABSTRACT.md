# Reproduce the circstats graphical abstract

The renderer uses the 599 × 342 pt graphical abstract from
[jgraving/circstats](https://github.com/jgraving/circstats/blob/main/figures/Graphical%20Abstract.pdf).
It keeps the paper's four-panel layout, vector artwork, code window, and
workflow arrows. Panel B uses a polar predictive plot and a compact
diagram showing how model draws are transformed.

## What the renderer computes

Panel B uses 16 illustrative angles digitized from the blue dots in the
paper. The script applies independent zero-mean normal priors to the two
real-valued predictors and fits an intercept-only von Mises example with
seeded random-walk Metropolis. It uses the package's circular modulo
and softplus functions to move from the predictor scale to the
parameter scale:

- **Direction μ** wraps into `[−π, π)` with circular modulo.
- **Concentration κ** becomes positive with softplus.

Highlighted points in each trace are plotted again on the circle or
positive axis after applying the link. The polar plot shows mean prior
and posterior predictive densities, pointwise 10th–90th percentile
bands across draws, and the observed angles around the rim. This is an
illustrative model generated from the measured dots. Its results describe
that example.

## Reference artwork and transparency

The experimental illustrations and trajectories in A, the model
diagnostics in C, the summaries in D, and the code window in B reproduce
the paper's measured vector artwork from `artwork.json`. These elements
are schematic reference drawings. The renderer reads the reference
PDF's nested graphics state with PyMuPDF to preserve the opacity and clip
of each retained path. The retained text is typeset from the measured
reference spans. The 148 original plot paths in panel B are omitted from
the export.

## Generate the exports

From the repository root, with Python 3.12 or newer:

```bash
python -m pip install -e '.[figure]'
python notebooks/circstats_graphical_abstract.py --font-family 'Liberation Sans'
```

Rendering runs offline with the committed reference PDF and vector asset.
The five exports go to `assets/circstats-graphical-abstract/`: editable
SVG text and vector paths, vector PDF, 1797 × 1026 px PNG, self-contained
HTML, and a JSON record of source and export hashes, fonts, and sampling
settings. Use `--output-dir` to choose another location. The optional
font override applies to proportional text. The editor window uses a
monospace font.

The CI workflow renders twice and compares all five outputs byte for
byte. The test suite checks the source asset hashes, the wrap and
softplus bounds, the predictive densities, and retained transparency in
the vector PDF.
