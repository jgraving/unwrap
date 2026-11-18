# Reproduce the Goldschtein bat figure

`notebooks/goldshtein_bats_publication_figure.py` rebuilds the circstats
paper's bat composite figure (`goldshtein_bat_figure.pdf`, panels a–h,
583 × 626 points) the way the zebra figure does: **programmatically,
from the notebooks' own plotting code, the pinned model fits and the
pinned data**, so changes to the notebook flow into the paper figure.
The committed `notebooks/assets/goldshtein_bats_publication/reference.pdf`
is the Inkscape-polished author artwork the layout is measured against;
only panel b and the panel letters have no notebook counterpart and are
pinned verbatim from it.

Each panel is rendered by the very `goldshtein_bats_plots.py` function
the `goldshtein_bats.ipynb` cell calls, with that cell's arguments, at
its natural matplotlib size:

- **a** — `plot_destination_panels_compact` (cell 23)
- **c** — `plot_total_elpd` (cells 57/58)
- **d** — `plot_delta_total_and_by_treatment` (cells 60/61)
- **e** — `plot_mean_hdi_horizontal` (cells 70/73)
- **f** — `plot_mixture_polar_pdfs` (cell 68)
- **g** — `plot_population_curve_smooth_compact` (cells 100/104)
- **h** — `plot_tracks_pgoal_small_bigfont` (cells 94/98)

…then embedded in the 583 × 626 pt page at the Inkscape-measured embed
scale, anchored on the panel's own rendered axes box (below). Panel b
(the tagged-bat photograph, the treatment-sense icon grid, the
"Treatment groups" label, the "Placeholder" watermark) plus the panel
letters a.–h are pinned chrome.

## Generate

From the repository root, with the notebook/figure environment:

```bash
python -m pip install -e '.[figure,test]'
python notebooks/goldshtein_bats_publication_figure.py \
    [--data-root DIR] [--download] [--output-dir DIR] \
    [--font-dir DIR] [--font-family FAMILY]
```

- The seven data/fit inputs are pinned to `jgraving/unwrap-data`
  revision `86c41b5f43eb61191405b9bc85baf730a492a06f` with SHA-256
  checksums (`INPUT_HASHES`): the two tracking workbooks, the Agamon–
  Hula lake outline, and the two cached fits (`bats_vm_v1`,
  `bats_vmu_v1`) with their metadata. `--download` fetches any missing
  file from that revision; existing files are verified in place and
  never overwritten.
- The committed artwork assets (five files) are SHA-256-verified
  (`check_assets`) and never modified.
- Fonts: the house style is Arial; `--font-dir` registers licensed font
  files, `--font-family` selects the family (a missing face raises
  instead of silently substituting). CI renders with
  `--font-family 'Liberation Sans'`; the committed preview is rendered
  with Arial.

Nothing is sampled or refit. The slim fit snapshots strip exactly the
per-observation log-likelihood (vMU also the per-observation
deterministics), which the generator recomputes from the stored draws
through the notebook's own restoration paths: the baseline's
log-likelihood from the rebuilt bambi graph (cell 41,
`graph_free_rvs` + `pm.compute_log_likelihood` against `model.backend.
model`, working around bambinos/bambi#1084); the mixture's
deterministics via `posthoc.restore_per_obs_deterministics` (the cell
53 data-prep `tag_idx`/`treat_idx`) and its log-likelihood from the
rebuilt raw-PyMC graph (cell 53, `sample=False`).

The outputs default to `assets/goldshtein-bats-publication/`
(committed for review):

| Output | Content |
| --- | --- |
| `goldshtein_bats_publication.svg` | Vector composite (583 × 626 pt, fonts as text) with the original basemap/photo JPEGs re-emitted as data URIs |
| `goldshtein_bats_publication.pdf` | 583 × 626 pt composite with the original basemap/photo JPEG streams restored onto their image objects, the photo's uniform 0.4684 SMask re-attached, dates cleared and `/ID` content-pinned |
| `goldshtein_bats_publication.png` | 1749 × 1878 px (3×) preview rendered from the PDF |
| `goldshtein_bats_publication.html` | Self-contained page with the published caption |
| `goldshtein_bats_publication.json` | Provenance: data revision and input hashes, asset hashes, fit metadata and convergence diagnostics, per-panel scales/placements, fonts, package versions and every export SHA-256 |

## Measured layout

The Inkscape composite embeds each notebook panel figure as a uniformly
scaled vector page. The per-panel scale is the reference line width of
the panel's own chrome divided by the notebook code's line width of the
same stroke (exact `artwork.json` values in `SCALE_*`):

| panel | stroke (code lw → reference lw) | embed scale |
| --- | --- | --- |
| a | border 1.0 → 1.0014435 | 1.0014435 |
| c, d, e | spines 0.8 → 0.9209633 | 1.1512041 |
| f | outer rings 0.8 → 1.0818636 | 1.3523296 |
| g | spines 0.8 → 1.1081346 | 1.3851683 |
| h | borders 1.5 → 1.2558924 | 0.8372616 |

The scales are cross-checked against the reference typography: the
embedded fonts equal the current-code fonts times the measured scale
(exactly for the panel-a title 9 pt × 1.0014435 ≈ 9 pt, the c ylabel
8 × 1.1512041 = 9.21 pt, the g tick 6 × 1.3851683 = 8.31 pt, the h
colorbar 9 × 0.8372616 = 7.54 pt), which is why embedding the current
output at these scales reproduces the reference typography while only
fit-derived numbers and current-cell labels deviate.

Anchors are the measured axes boxes (page points, y down), from the
reference image placements and spine/border boxes; the generator places
each panel by translating its rendered first axes box onto the measured
box:

| panel | measured box(es) (x0, y0, x1, y1) | source |
| --- | --- | --- |
| a | (17.15, 16.48, 105.12, 238.91) / (106.81, 16.48, 194.78, 238.91) | `image_45` placements |
| c | (404.10, 25.87, 565.07, 90.61) | spine boxes |
| d | (404.10, 135.22, 565.07, 215.50) | spine boxes |
| e | (34.26, 268.96, 195.23, 338.88) | spine boxes |
| f | four polar boxes, e.g. (20.735, 385.740, 118.795, 483.800) | outer ring centres, ring radius 2.0 of `RMAX` 2.03 |
| g | (281.68, 259.09, 575.00, 340.46) | spine boxes |
| h | four axes boxes, e.g. (249.72, 378.30, 330.70, 583.06) | `image_23` placements |
| h colorbar bar | (364.46, 586.24, 465.74, 591.26) | `image_36` placement |
| b photo | (204.74, 6.62, 305.06, 92.57), alpha 0.4684 | `image_279` placement |

Panel h's colorbar needs one extra trick: the reference bar is 101.28
pt wide (figure fraction 0.30) while the current cell's default is
0.60. The colorbar axis box is a figure fraction, so `draw_panel_h`
solves it in two deterministic passes: pass one at the cell's own
fractions fixes the panel's translate (anchored on the first axes
box); pass two passes the fraction that lands the gradient bar exactly
on the measured bar.

The two flight-map panels take their CartoDB basemaps from the
published figure's own tile mosaics (`image_45` / `image_23`) through
the additive `basemap_img` parameter of the map functions (an
`imshow` over the unified data extent, skipping the network fetch);
the identical tracking data gives identical extents, so the stretch
matches the reference placement.

## Pinned chrome

From `artwork.json` (the full verbatim extraction: 459 drawings, 8 image
objects / 12 placements, 64 spans; schema checked on load):

- the panel-b icon grid: the 69 drawings whose bounding-box centre sits
  in the band below the photo and above panel g, under the photo's x
  extent (count asserted);
- the chrome text: the panel letters a.–h (20 pt bold), "Treatment
  groups" (10 pt) and the "Placeholder" watermark (17.4 pt, alpha
  0.4684, rotated along its measured line direction);
- the photograph JPEG, drawn into the exports with a uniform 0.4684
  SMask (the reference stores no SMask for it; its "1.4cm" ruler is
  baked into the JPEG).

## Expected deviations from the reference

The reference was assembled from an older notebook state; by design the
rebuilt figure follows the **current** cell code and the **pinned** fits:

- panels c/d ELPD values and deltas (reference ≈ −6800 / −6890, Δ≈90;
  pinned fits −6710.4 ± 60.3 / −6829.2 ± 59.3, Δ 118.8), panel e's
  density HDIs, panel f's loops and panel g's curves move with the
  pinned fits;
- panel e's xlabel reads "Predictive density (mixture)" (the current
  cell); panel h's colorbar reads "Posterior directedness" and its
  tracks are coloured by `evi_post_mean` (the current cell calls);
- panel d shows the current per-treatment labelled rows (the older
  reference shows bracketed "Within-treatment");
- panel h's N label and scale bar sit in the rightmost subpanel (the
  current code's placement; the reference has them in subpanel 1);
- panels a/h basemaps and tracks and the panel-h colorbar are
  code-drawn vector/gradient rasters where the Inkscape reference
  rasterised them (its four track overlays and their SMask channels
  and the viridis gradient raster are no longer committed assets).

Layout fidelity (panel positions, scales, fonts) is the match
criterion; the 3× raster diff against `reference.pdf` highlights the
fit-dependent deltas above.

## Determinism

- SVG: `svg.hashsalt = "unwrap-goldshtein-bats-publication-v1"`,
  `savefig(format="svg", metadata={"Date": None})`;
- PDF: creation/modification dates cleared, the trailer `/ID` normalised
  to a constant before hashing and set to the SHA-256 of the document
  (so the random ID MuPDF stamps on write leaves no trace);
- `plot_mixture_polar_pdfs` draws at `rng(123)` (≤ 4000 samples); both
  fits come from the cached snapshots (no sampling);
- environment: `MPLBACKEND=Agg`, `OMP_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`.

Two renders from the same checkout produce byte-identical SVG and PDF
(verified in the CI double-render job, under
`--font-family 'Liberation Sans'`; the committed snapshots in
`assets/goldshtein-bats-publication/` are rendered with Arial).

## Tests

`tests/test_goldshtein_bats_publication_figure.py` is hermetic (no
download, no licensed fonts, no sampling): asset-pin and artwork-
structure checks (459 drawings / 8 images / 12 placements / 64 spans,
the panel-b curation and chrome-span sets), the layout constants
against the asset's measured image placements, `prepare_inputs`
tamper-proofing in `tmp_path`, the placement/affine/bar-solver math,
the SVG id-prefixing and basemap-swap on synthetic parts, the
`_stabilize_pdf` unit test, and the chrome layer rendering from the
pinned artwork on the default font stack.

CI (`.github/workflows/goldshtein-bats-figure.yml`) runs the hermetic
tests on relevant pull requests. Run the workflow manually to render the
full pipeline twice (Liberation Sans; first render `--download`) and
compare both output directories.
