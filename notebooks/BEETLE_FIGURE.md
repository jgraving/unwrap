# Reproduce the beetle publication figure

The six-panel figure combines the original arena illustration with newly
rendered statistical panels. It retains Arial and the existing grey/gold
palette (`#9e9e9e` / `#f4a300`), with aligned headings, restrained axes,
consistent spacing, and a shared condition legend. The statistical
estimands and 89% highest-density intervals match `beetles_plots.py`.

A rendered review snapshot (SVG, PDF, PNG, HTML, and provenance) is checked in
under [`assets/beetles-publication/`](../assets/beetles-publication/). Regenerate
that snapshot with `--output-dir assets/beetles-publication`.

## Generate

From the repository root, in a Python 3.12+ environment:

```bash
python -m pip install -e '.[figure]'
python notebooks/beetles_figure.py --download
```

The command downloads only the raw beetle table, cached fit, and fit metadata
(about 4 MB). It pins `jgraving/unwrap-data` revision
`31eb4bdb4a0fe9a0c0bd20df7f101dbd884a8dcc` and verifies all three SHA-256
checksums before reading them. It does not sample or execute the notebook.
Subsequent runs work offline. Existing files are never replaced by a download.

To use an existing copy of that data revision:

```bash
python notebooks/beetles_figure.py --data-root /path/to/unwrap-data
```

Arial regular and bold must be installed, or supplied locally:

```bash
python notebooks/beetles_figure.py --download --font-dir /path/to/arial-fonts
```

Fonts are not distributed with the repository. Missing Arial is an error,
so exports cannot silently change the paper's typeface. An intentional
alternative is possible with `--font-family 'Liberation Sans'`; that choice
is recorded in the provenance file. Text is converted to vector outlines
in SVG and PDF, preserving the chosen typeface on other computers.
CairoSVG also requires the Cairo system library (see its installation
instructions if your platform does not already provide it).

## Outputs

By default, outputs go to the gitignored `figures/beetles/` directory.
Choose another directory with `--output-dir`.

| File | Purpose |
| --- | --- |
| `beetles_publication.svg` | Complete vector figure, including arena artwork |
| `beetles_publication.pdf` | Vector manuscript export, same layout as SVG |
| `beetles_publication.png` | 1800-pixel review preview |
| `beetles_publication.html` | Self-contained figure and caption for browser review |
| `beetles_publication.json` | Input/source/font/output hashes, package versions, numerical summaries |

The SVG is deterministic for the same inputs, code, fonts, and rendering
environment. The PDF may contain renderer creation metadata, so scientific
and visual reproducibility does not imply identical PDF bytes. The JSON
records the git revision **and hashes of the actual source files**, including
uncommitted changes, plus the installed dependency versions. For exact
environment replication, preserve the JSON and pin those versions.

## What each panel means

- **a:** The original schematic of repeated arena exits. Ten beetles have
  ten trials in each condition; four light-on bearings are missing.
- **b:** Von Mises kernel densities of the posterior population mean
  direction, using the notebook's bandwidth concentration of 80 and a common
  radial density scale. Dots are **observed circular mean bearings per beetle**,
  as in the notebook, not fitted individual means. The curve width represents
  uncertainty about the population direction, not the dispersion of exits.
- **c:** Population-level within-individual concentration, expressed as
  `I1(kappa) / I0(kappa)` after applying the softplus link to the fixed effects.
- **d:** The resultant length across the ten fitted individual mean directions,
  calculated separately for every posterior draw and condition.
- **e:** Von Mises probability within 45° of the southwest floodlight on either
  side. It uses the population direction and concentration, with individual
  offsets set to zero; it is not a predictive distribution averaged over
  individuals. The uniform reference is 0.25.
- **f:** On-minus-off sector probability, paired within posterior draw. Zero
  indicates no change. Its axis expands below zero if the displayed mean or
  interval extends below zero.

Panels c–f show posterior means and 89% highest-density intervals, using the
same helper as the notebook. There are no new priors, sampling settings,
bootstrap estimates, or refits. The current pinned posterior has four chains
with 2,000 retained draws each and zero divergent transitions; the generator
refuses a cache with divergent transitions. A new fit/data revision requires
an explicit revision-and-checksum update and renewed scientific validation.

The HTML caption is suitable as a starting point for the manuscript caption.
It explicitly distinguishes observed means in b from fitted means in d,
and defines the population-level probability in e–f.

## Validation

```bash
python -m pytest tests/test_beetles_figure.py
```

The tests cover circular boundary crossing, matched posterior contrasts,
dimension ordering, an analytic uniform distribution, data integrity, and
vector export. Review the generated PNG/PDF at manuscript size after visual
changes. The original notebook and its outputs remain intact.
