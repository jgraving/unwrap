# Reproduce the circstats graphical abstract

The renderer uses the 599 × 342 pt graphical abstract from
[jgraving/circstats](https://github.com/jgraving/circstats/blob/main/figures/Graphical%20Abstract.pdf).
It keeps the paper's four-panel layout and equally weighted workflow
arrows. Panel B enlarges the original Mac-style editor window and centers
`model.py` in its title bar. It also uses a polar predictive plot and a compact
diagram showing how model draws are transformed.

## What the renderer computes

Panel B starts from the 16 blue angles digitized from the paper and uses
70% of each angle's offset from their circular mean. This gives the
illustrative data a clearer directional cluster. The script fits an
intercept-only von Mises example with seeded random-walk Metropolis.
It uses independent normal priors with means (0, 1.5) and standard
deviations (3.5, 2.5) on the direction and concentration predictors.
The higher concentration-predictor mean raises the mean κ of the seeded
prior draws from 1.25 to 2.16. The individual curves peak at different
angles. The package's circular modulo and softplus functions transform
the predictors to the parameter scale:

- **Direction μ** wraps into `[−π, π)` with circular modulo.
- **Concentration κ** becomes positive with softplus.

Three independent chains start one turn apart. The trace plots show 400
consecutive draws per chain in translucent, chain-specific colors.
For display, each direction trace uses one equivalent turn, with
within-turn variation expanded consistently across the three bands.
Eighteen marked draws per chain carry through each link curve to the
circle or positive axis, alongside the highlighted examples. The
polar plot gives each of 36 prior and 44 posterior predictive draws a
faint line. A thin blue curve averages all posterior predictive draws,
and a dashed circle marks the uniform reference. Zero density is at
radius 1, followed by a common linear density scale. The observed-angle
display uses four translucent, jittered copies of each fitted angle, with
a small gap between the dots and the inner ring.

The illustrative fit uses seed 159, 600 warmup iterations, and 1,000
retained draws per chain. Each iteration includes one direction-only
proposal with standard deviation 0.42 and two concentration-only
proposals with standard deviation 1.6. Tested on the original prior,
the extra direction-only proposal raised tail ESS for cos μ from 443 to
898. After raising the concentration-predictor prior mean, rank R-hat
and bulk/tail ESS passed checks of 1.01 and 400 on all identifiable
quantities. The posterior median κ was 3.68. HMC divergence counts are
inapplicable to random-walk Metropolis.

| Quantity | Rank R-hat | Bulk ESS | Tail ESS |
| --- | ---: | ---: | ---: |
| sin μ | 1.0015 | 1,171 | 1,491 |
| cos μ | 1.0017 | 1,017 | 1,226 |
| κ | 1.0004 | 1,157 | 908 |

## Reference artwork and transparency

The experimental illustrations and trajectories in A, the model
diagnostics in C, the summaries in D, and the code window in B reproduce
the paper's measured vector artwork from `artwork.json`. These elements
are schematic reference drawings. The renderer reads the reference
PDF's nested graphics state with PyMuPDF to preserve the opacity and clip
of each retained path. The retained text is typeset from the measured
reference baselines using shared sizes for panel titles, labels, and
legends. The model comparison and hypothesis plots use equally tall
y-axes and centered 7 pt axis titles.
Workflow arrowheads use the paper PDF's round joins. The polar
graphs in B and D share a guide style, with spokes between the two gray
rings. The 148 original plot paths in panel B are omitted from the export.

## Generate the exports

From the repository root, with Python 3.12 or newer:

```bash
python -m pip install -e '.[figure]'
python notebooks/circstats_graphical_abstract.py --font-family Arial --font-dir /path/to/arial-fonts
```

Rendering runs offline with the committed reference PDF and vector asset.
The five exports go to `assets/circstats-graphical-abstract/`: editable
SVG text and vector paths, vector PDF, 1797 × 1026 px PNG, self-contained
HTML, and a JSON record of source and export hashes, fonts, and sampling
settings. Use `--output-dir` to choose another location. The checked-in
exports use Arial regular and bold. You can omit `--font-dir` if Arial is
installed. CI renders with Liberation Sans. The editor window uses a
monospace font.

The CI workflow renders twice and compares all five outputs byte for
byte. The test suite checks the source asset hashes, the wrap and
softplus bounds, the predictive densities, and retained transparency in
the vector PDF.
