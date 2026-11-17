# Reproduce the locust density–order figure

The generator starts from the original manuscript's 750 × 525-point layout
and applies a restrained presentation pass. Image crops, proportions,
panel d's calibrated 10 cm scale bar, and panel e's complete vector schematic
retain their original geometry. Panels f–h use the same numerical inputs
and positions, with the original title wording and full predictor names.

The current styling uses 16-point panel letters, bold 11-point f–h titles,
9-point statistical labels, thin left/bottom axes, and vertical legends
without frames. Prediction lines are 1.5 points with 0.16-alpha bands.
Borders around d/e are lighter. Order legend entries use 1.0 through 0.0.
Density units are single math expressions, `(m⁻²)` and `(rad⁻¹)`, keeping
the exponent inside the parentheses.

Arial and the original red/purple, brown/green, and blue/orange palettes are
retained. Statistical curves, bands, axes, schematics, and typography remain
vector. Photos/renderings and the source scale-label backing/mask retain
their original raster content. No generated imagery or screenshots of
statistical panels are used. The faithful reproduction baseline is retained
at commit `fb1ad2d60e58facca79c16d597d0d9fbc2b1a9f8`; the current rendering
applies the approved polish to that baseline.

## Generate

From the repository root, with Python 3.12+:

```bash
python -m pip install -e '.[figure]'
python notebooks/locust_density_figure.py --download --font-dir /path/to/arial-fonts
```

Arial regular and bold must be installed or provided with `--font-dir`.
The script refuses a silent font substitution. An explicit alternative is
possible with `--font-family 'Liberation Sans'`; provenance records it.
That override affects redrawn labels; source artwork retains its original
outlined typography.
Font files are not redistributed. SVG/PDF lettering is outlined to preserve
the selected face on other computers. CairoSVG requires the Cairo system
library.

The first run fetches only `locust_density_hier.nc` and its metadata (about
29 MB) from `jgraving/unwrap-data` revision
`31eb4bdb4a0fe9a0c0bd20df7f101dbd884a8dcc`. SHA-256 checksums are verified
before reading. Subsequent runs work offline. Existing files are never
overwritten by a download. No fitting, notebook execution, or full tracking
download is required. The cache has four chains × 1,000 retained draws and
zero divergent transitions; the command rejects divergent caches.
On the eight fixed effects plus the concentration-offset scale, the pinned
cache has maximum rank-normalized R-hat 1.008955, minimum bulk ESS 729.66,
and minimum tail ESS 1150.57 (ArviZ 1.3.1). These are checks of the existing
fit, not a newly sampled fit.

To use an existing data checkout:

```bash
python notebooks/locust_density_figure.py --data-root /path/to/unwrap-data --font-dir /path/to/arial-fonts
```

The default output is the ignored `figures/locust-density/` directory.
The review snapshot in this PR is regenerated with:

```bash
python notebooks/locust_density_figure.py --download --font-dir /path/to/arial-fonts --output-dir assets/locust-density-publication
```

| Output | Content |
| --- | --- |
| `locust_density_publication.svg` | Vector figure with embedded original JPEG streams |
| `locust_density_publication.pdf` | Manuscript export with vector statistics and original raster photos |
| `locust_density_publication.png` | 2,200-pixel review preview |
| `locust_density_publication.html` | Self-contained review page and full caption |
| `locust_density_publication.json` | Input, source, font and output hashes; versions; numeric summaries |

SVG output is deterministic for identical code, fonts, inputs, and environment.
PDF creation metadata may differ. The JSON records the git revision and actual
source hashes, so uncommitted changes remain identifiable. Preserve and pin
the recorded package versions to reproduce a rendering environment exactly.

## Source images

Seven **native JPEG image streams** were extracted from page 1 of
`figures/locust_density_order_figure.pdf` in `jgraving/circstats`, revision
`9bef710eb44fd8d20a61a7119a12e4933d02ce6d`. Their dimensions, PDF object
numbers, and checksums are in
[`assets/locust_density_order/images.json`](assets/locust_density_order/images.json).
The PDF SHA-256 is
`48609fdc9185e20227ae0b1c20dd377301d5aa962c7d57078ef1bf5873d0bf21`.
The JPEG files are checked in, so normal generation does not need access to
the manuscript repository. `original_panels.svg` additionally preserves
the original a–e artwork and all PDF clipping/image transforms, rather than
inventing new image crops. The source image files remain intact.

The preparation command removes the old f–h statistical drawings while
keeping their panel letters/titles, and removes the old unit line. It
exports the retained source artwork as SVG paths and
images. A pixel comparison at two pixels per point verifies **zero changes**
throughout a–e, excluding only the unit line being replaced. The scale bar
keeps its original endpoints, 139.509506 to 163.978531 PDF points, alongside
the unmodified photo transform and 10 cm label. This avoids recalibrating a
scale bar after changing an image crop.

The baseline asset stays unchanged by styling. The renderer parses a copy,
shrinks only the identified panel-letter glyphs, replaces f–h title glyphs
with bold typesetting of the same wording, and lightens only the d/e frame
strokes. Targeted structure checks guard against editing unexpected paths.
Every source image attribute and ancestor transform/clipping attribute,
and every scale-bar path attribute, is checked for equality before and after
this presentation pass. Source-preparation pixel comparisons refer to the
preserved baseline, before these deliberate typographic/frame edits.

To repeat the extraction from a local copy of the **pinned** source PDF and
notebook (PyMuPDF is needed only for this optional preparation step):

```bash
python -m pip install 'pymupdf>=1.24,<2'
git show 075cb7274f8d80f284346dbe30bdc3e0f20dfd89:notebooks/locust_density_order.ipynb > /tmp/locust-density-source.ipynb
python tools/prepare_locust_figure_assets.py --pdf /path/to/locust_density_order_figure.pdf --notebook /tmp/locust-density-source.ipynb
```

The command validates both source checksums and extracts image streams
without resampling. It also parses the saved screening result described
below. Photo credit follows the manuscript: Centre for the Advanced Study
of Collective Behaviour, used with permission. No broader licence is
inferred from extraction; keep the credit with the figure.

## What the statistical panels show

**f — Model comparison.** Five specifications are compared on **440 recording
mean directions**, one per retained recording. The numerical input is the
saved output of cell 23 in the Python notebook at unwrap revision
`075cb7274f8d80f284346dbe30bdc3e0f20dfd89`, automatically parsed into
`assets/locust_density_order/selection.json`. It includes ELPD, differences
from the best model, standard errors, maximum Pareto k, and source provenance.
Bars are ±1 SE. Its six-decimal print precision is retained. Density and
Order are written out in full, including both interaction labels.

This panel reproduces a **saved screening table**, rather than recomputing
PSIS-LOO from cached Python draws: unwrap-data contains the separate initial
screening fit but does not contain the five Python sweep fits. The exact
original recomputation path remains `locust_density_posthoc.model_sweep` on
`posthoc.screening_table(merged_df)`, using draws=1000 and seed=123, after
the notebook's original preprocessing. This requires the full raw tracking
parts and repeats the five original fits. No new screening run is claimed.

**g — Fixed effects.** The saved hierarchical fit models **18,641 retained
headings grouped into 440 recordings**. Points and bars show population
fixed-effect posterior means and 95% highest-density intervals. The mean
for the direction intercept is circular. Coefficients describe direction
in radians before wrapping, or concentration before softplus. Predictor
coordinates are `o = 2*order - 1`, `d = log2(density)/3 - 1`, and `d*o`.
The graph check described below confirms the stored coefficient coordinates.

**h — Density in the march direction.** Reference predictions set the
recording concentration offset to zero. They do not average fitted offsets
or integrate over a new recording's offset. The curves show
`exp(mean(log f(0)))`, with 90% HDIs formed on log density then exponentiated,
matching the notebook's intended geometric summary. The label explicitly
uses density per radian. This quantity is neither a probability nor an
ELPD score; the dashed uniform reference is `1/(2*pi)`. Curves show uncertainty
in the reference prediction, rather than the spread of individual headings.

## Correcting prediction coordinates

The current notebook explicitly warns that its manual helpers subtract
training means from **already uncentred** stored Bambi coefficients, then
form an interaction from centred marginal predictors. That changes the
fitted surface. The new generator corrects this reconstruction: it uses the
stored intercept with the original scaled main effects and their product.
The scientific fit, priors, and sampling settings are untouched.

This follows the fitted Bambi 0.21 implementation: its backend subtracts
`mean(design columns) @ slopes` from the sampled intercept when exporting
the posterior. To evaluate its centred PyMC graph again,
`Model._re_center_intercept` adds that factor back. Consequently, predictions
from the exported coefficients use the **original** design columns.
`tests/test_locust_density_figure.py` verifies this against both components
of a built Bambi/PyMC graph, without sampling, using an unbalanced design
where the interaction mean differs from the product of marginal means.
It evaluates the actual circular-modulo and softplus links, bypassing the
known Bambi von Mises prediction special case.

This correction affects the new publication curves. The existing notebook
helpers, notebook prose, and saved notebook outputs are preserved in this PR.
Their broader prediction/PPC correction and rebake remain separate work.
Do not interpret a visual difference from the manuscript's old curves as a
new fit or a change in the model specification.

## Validation

```bash
python -m pytest tests/test_locust_density_figure.py -q
```

Checks cover the Bambi/PyMC graph comparison, stable density against SciPy,
HDIs against ArviZ, periodicity, dimension order, the analytic uniform case,
input integrity, original image hashes, deterministic SVG generation, and
SVG/PDF export, unique composed SVG IDs, and the checksum of the preserved
a–e artwork. The presentation check verifies photo geometry and calibrated
scale-bar preservation. Source preparation also verifies pixel preservation
of the original baseline panels. Review the PNG and an independently rendered PDF at
intended manuscript size. This PDF contains the source photographs and the
original scale-label backing/mask; statistical content remains vector.
