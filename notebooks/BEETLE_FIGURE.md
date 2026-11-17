# Reproduce the beetle publication figure

The six-panel figure combines the original arena and floodlight illustrations with newly
rendered statistical panels. Panel a composes matched light-off and light-on
arenas from the same SVG paths, at the same scale, with the floodlight beside
the light-on arena. The grey/gold exit paths illustrate the measurement
protocol: identical schematic geometry is reused across conditions. Each
vector copy receives unique IDs. All placement and recolouring are in
`beetles_figure.py`; no manually edited or generated illustration is needed.

The figure retains Arial and the existing grey/gold palette
(`#9e9e9e` / `#f4a300`), with concise panel headings centred over the plotting
areas and measurement labels on the c–f y-axes.
The polar plot has light crosshairs and N/E/S/W labels. A condition legend
sits above the plot, while a separate key beside it distinguishes the
population-level mean curves from the individual-level mean points.
White space separates the two rows.
The statistical panels have ticked axes with only zero and uniform reference
lines. Explanatory details live in the caption. The statistical
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

- **a:** Matched light-off and light-on schematics show the two conditions
  experienced by the same ten beetles. The original arena and beetle artwork
  is reused in both; exit paths are illustrative, not observed trajectories.
  Each beetle has ten trials per condition; four light-on bearings are missing.
- **b:** Von Mises kernel densities of the posterior population mean
  direction, using the notebook's bandwidth concentration of 80 and a common
  radial density scale. Dots are **observed circular mean bearings per beetle**,
  as in the notebook, not fitted individual means. The curve width represents
  uncertainty about the population direction, not the dispersion of exits.
- **c:** Population-level within-individual concentration, expressed as
  `I1(kappa) / I0(kappa)` after applying the softplus link to the fixed effects.
- **d:** The resultant length across the ten fitted individual mean directions,
  calculated separately for every posterior draw and condition.
- **e:** Von Mises probability within the 157.5°–247.5° analysis sector,
  centred on the nominal south-southwest bearing (202.5°). It uses the
  population direction and concentration, with individual
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

## Check against the original experiment

The original study is Foster et al. (2021),
[*Light pollution forces a change in dung beetle orientation behavior*](https://doi.org/10.1016/j.cub.2021.06.038).
Its Figure 1 caption describes an **upward sky view** with north at the top,
east on the left, south below, and west on the right. Figure 3 uses that same
convention. This is a reflected view of the sky, not the clockwise plan view
used by this figure and the beetle analysis notebook.

The authors' [bearing-analysis script](https://github.com/JJFosterLab/light-pollution/blob/2dadc3ab3034251e190b04fd9b3cce4cc618669a/Behaviour/LPHeadings-SecurityLightBergsig.R#L42-L54)
explicitly puts zero at the top and makes angles increase clockwise. Their
reported easterly heading of 74° also establishes the compass interpretation
of these values. The original `SecurityLightBergsig.txt` is byte-identical to
the pinned `unwrap-data` input (SHA-256 recorded above). Thus, the current
plot's north-up, east-right orientation is consistent with the bearing data;
it is not intended to reproduce the reflected sky photographs. The N/E/S/W
labels and light crosshairs follow this convention, which is also stated
in the caption.

**Lamp position versus analysis sector:** the original paper describes the
floodlight as south-southwest and reports a collective **beetle heading** of
210°. That statistic is not an independent measurement of the lamp azimuth.
The previous 225° (−135°) centre represented southwest. We replace it with
**202.5°**, the nominal compass bearing of south-southwest. This is an
explicit approximation from the directional description, **not a recovered
surveyed azimuth**. Neither the article nor the authors' bearing-analysis
script specifies an exact lamp azimuth. The archived Figshare image data
were also inspected, but the files inspected did not establish a compass
orientation for measuring the lamp independently. No pixel-derived bearing
is claimed.

[`beetles_study.json`](beetles_study.json) is the single setting shared by
the publication generator, Python notebook/helper, and R notebook/Quarto
twin. The icon and arc move together, and both probabilities are integrated
over that same sector. A surveyed bearing, if obtained, should replace this
approximation there and trigger regeneration. The caption and provenance
identify its basis. The statistical intervals condition on this centre;
they do not propagate uncertainty about the lamp position.

### Sector centre sensitivity

Using the **same pinned Python posterior**, with the 45° half-width fixed:

| Centre | Off probability | On probability | On − off | 89% HDI of difference |
| --- | ---: | ---: | ---: | ---: |
| 191.25° | 0.328 | 0.758 | 0.430 | −0.183 to 0.917 |
| **202.5° (nominal SSW)** | **0.229** | **0.778** | **0.550** | **0.110 to 0.932** |
| 213.75° | 0.149 | 0.756 | 0.607 | 0.231 to 0.930 |
| 225° (previous SW setting) | 0.091 | 0.694 | 0.603 | 0.218 to 0.929 |

191.25° and 213.75° bracket the conventional 22.5° compass sector for SSW;
they are a sensitivity range, **not a confidence interval for the lamp**.
The positive mean contrast persists, but its HDI crosses zero at the lower
endpoint. Therefore, the evidence for a strictly positive sector contrast
is sensitive to the assumed centre. Panels c–d and the posterior direction
curves are unchanged; only the marked sector and e–f estimates change.

Reproduce this table from the repository root after generating the figure:

```python
import sys
sys.path.insert(0, "notebooks")
import numpy as np
import xarray as xr
from beetles_figure import DATA_REVISION, posterior_arrays, _mean_hdi

path = f".cache/beetle-figure/{DATA_REVISION}/fit_cache/beetles_fit.nc"
with xr.open_dataset(path, group="posterior", engine="h5netcdf") as ds:
    post = ds.load()
for degrees in (191.25, 202.5, 213.75, 225):
    arrays = posterior_arrays(post, sector_center=np.deg2rad(degrees))
    print(degrees, arrays["mass"].mean(axis=1),
          _mean_hdi(arrays["difference"], hdi_prob=0.89))
```

## Notebook rebake handoff

Updated sources are `notebooks/beetles.ipynb`, `notebooks_r/beetles.ipynb`,
and the generated `notebooks_r/beetles.qmd`. Only the changed Python plot
cells (zero-based 16–17) and R summary cell (19) have their stale outputs
and execution counts cleared. All other cells, outputs, model priors,
sampling specifications, and notebook metadata are preserved.

For the local rebake:

1. Pull this PR branch and use the existing study environments and data
   checkout. The sector change is post-processing; it does not require a
   refit. To reuse matching caches, set `UNWRAP_PY_FIT_CACHE=1` and
   `UNWRAP_PY_FIT_CACHE_DIR=/path/to/unwrap-data/fit_cache` for Python, and
   `UNWRAP_R_FIT_CACHE=1` and `UNWRAP_R_FIT_CACHE_DIR` to that same directory
   for R. Confirm cache hits for `beetles_fit` and `beetles_fit_v2`; a cache
   miss follows the notebooks' existing refit behavior. Keep the full-budget
   settings, and report whether outputs came from a reused fit or a refit.
2. Execute and save both beetle notebooks. Preserve all model and sampling
   settings and check their diagnostic gates. If a gate fails, record it
   before changing any settings, following `AGENTS.md`.
3. Check that both use 202.5° ±45° and a uniform mass of 0.25. Compare the R
   summaries to its own fit; the numbers in the sensitivity table above
   come from the pinned Python fit. The plots now allow negative difference
   intervals instead of clipping them at zero.
4. Run `python tools/ipynb_to_qmd.py`, then
   `python tools/ipynb_to_qmd.py --check`. The `.ipynb` is the R source of
   truth; do not edit the generated `.qmd` separately.
5. Rebuild the publication exports with the command above. A different fit
   intentionally requires updating its pinned revision and checksums.

The revised publication figure was rebuilt from the saved Python fit.
Full notebook execution and R rendering are left for the local rebake;
they are not represented as completed checks in this PR.

## Validation

```bash
python -m pytest tests/test_beetles_figure.py
```

The tests cover circular boundary crossing, independent numerical
integration at the revised centre (including the R notebook's quadrature),
matched posterior contrasts, dimension ordering, an analytic uniform
distribution, negative difference intervals, data integrity, and vector
export. Review the generated PNG/PDF at manuscript size after visual changes.
