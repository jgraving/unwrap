# Reproduce the zebra escape-response figure

`notebooks/zebra_publication_figure.py` rebuilds the circstats paper's
zebra composite figure from stored data only. The committed
`notebooks/assets/zebra_publication/reference.pdf` (900 × 1080 points) is
a reduced copy of the author-supplied composite. The layout is
measured from it, and the exports redraw panels c–f from the cached
posterior while retaining the composite's page geometry, typography and
palette.

- **a** — the escape map from PR #154 (`draw_escape_map` in
  `notebooks/zebra_plots.py`): the aerial mosaic, the 18 tracked
  trajectories of `observation088` (focal in black, group in white) and
  the focal animal's body-orientation arrows every 5 m of straight-line
  travel (orange pre-approach, purple post). The two circular callout
  insets show the rotated photographs and arrows from the published artwork.
  The cone tangents connect the data-derived anchors to the inset rings.
- **b** — the species photograph (author-provided asset).
- **c** — group-level mean direction: posterior KDE of the equally
  weighted mean of the 18 fitted individual mean directions, radial scale
  0–1 (relative density).
- **d** — group-level alignment: resultant length of the 18 fitted
  individual mean directions, posterior mean with 95% HDI.
- **e** — individual-level mean directions: per-track posterior KDEs plus
  each animal's observed circular-mean heading.
- **f** — individual-level concentration: per-track posterior mean
  resultant length (A(kappa)) with 95% HDI, paired within animal
  (house jitter, seed 123, sd 0.15).

c–f come from `zebra_posthoc.summary_stats` — the same path the
`zebra_escape.ipynb` notebook uses, so the composite and the notebook
figure cannot disagree on the numbers.

## Generate

From the repository root, with Python 3.12+:

```bash
python -m pip install -e '.[figure]'
python notebooks/zebra_publication_figure.py --download --font-dir /path/to/arial-fonts
```

Arial regular and bold must be installed or provided with `--font-dir`;
the script refuses a silent font substitution. An alternative family is
possible with `--font-family 'Liberation Sans'` (what CI renders); the
provenance JSON records it. Font files are not redistributed, and all
lettering is outlined to paths in the SVG/PDF. CairoSVG requires the
Cairo system library.

The first run fetches 48 pinned files (about 488 MB) from
`jgraving/unwrap-data` revision
`86c41b5f43eb61191405b9bc85baf730a492a06f` (the zebra-map-data
merge; the zebra fit blobs in it are byte-identical to the beetle/locust
figures' pin `31eb4bdb…`): the fit cache (`zebra_escape_fit.nc` +
metadata), the block table, the mosaic manifest and eight crop tiles, the
18 position parquets and the 18 raw-bearing parquets of `observation088`.
SHA-256 checksums are verified before reading; existing files are never
overwritten. No fitting, sampling, or notebook execution happens.

To use an existing data checkout instead:

```bash
python notebooks/zebra_publication_figure.py --data-root /path/to/unwrap-data --font-dir /path/to/arial-fonts
```

The outputs default to `assets/zebra-publication/` (committed for review):

| Output | Content |
| --- | --- |
| `zebra_publication.svg` | Vector figure, lettering outlined, four raster embeds |
| `zebra_publication.pdf` | 900 × 1080 pt composite with the same four JPEGs, transforms and clips as the SVG |
| `zebra_publication.png` | 1800 × 2160 px preview rendered from the PDF |
| `zebra_publication.html` | Self-contained page with the caption |
| `zebra_publication.json` | Provenance: data/code revisions, hashes, fonts, packages, fit diagnostics, layout constants, and every summary number below |

## Fit and data

- Fit: hierarchical von Mises GLM on the 5,417 block-averaged headings
  (3,222 pre / 2,195 post) of the 18 tracked `observation088` zebra;
  4 chains × 1,000 draws, tune 2,000, `target_accept` 0.8, seed 42
  (recorded in the fit metadata), **zero divergent transitions**.
- Pinned-cache diagnostics recomputed with ArviZ: maximum rank-
  normalized R-hat 1.004215, minimum bulk ESS 1351.49, minimum tail ESS
  1886.37 (ArviZ 1.3.0; recorded per variable in the provenance JSON).
  The notebook's `circular_summary` gate
  (`ess_bulk < 400 | r_hat > 1.01`) returns zero flagged rows; the
  generator rejects the fit otherwise.
- Map window: UTM 37N `x ∈ [267525, 267680]`, `y ∈ [37300, 37770]`,
  mosaic read window `y ∈ [37290, 37770]` (the 10 m strip at the bottom
  is clipped by the panel), frames capped at 90,698; mosaic assembled at a maximum
  side of 2600 px, then reduced to 1940 px before the map window is
  clipped. The visible mosaic is 627 × 1900 px at about 129 ppi.
  The source mosaic resolution is 0.02707 m/px. The calibrated 20 m scale bar and
  the north arrow sit at their measured positions (45.1 pt = 20 m at the map's
  0.4430 m/pt).
- The trajectories are vector. In the exports their vertices are
  decimated at 1 m of straight-line travel (`MAP_VERTEX_SPACING_M`,
  passed as `draw_escape_map`'s `decimate_m`; one meter is 2.26 pt at
  the panel's 0.4430 m/pt, and the published composite flattens its own
  trajectories into the ~640-px-wide raster, so no finer detail is
  carried). The 18 tracks cover ~48.4 km of travel, and 1 m of vertices
  keeps the map layers at ~48k vector segments -- the size at which the
  SVG backend (one element per collection segment) stays in the house
  artwork range. The heading arrows (every 5 m) are computed from the
  full-resolution trajectory exactly as in the map section of `zebra_escape.ipynb`.
- The pre/post arrow split on the map uses the per-frame `pre_scare`
  flag of `load_zebra_bearings_raw()`, exactly as the map section of `zebra_escape.ipynb`
  does: it flips at frame **53,902** for all 18 `observation088` tracks.
  `scare_frames.csv` lists `observation088` later (frame 72,703); the
  discrepancy is recorded in the provenance and left as published.

## Inset anchor frames

The published figure's callout cones are a presentation choice: their
apexes are not recoverable from the data by any sampling or turn-rate
rule. Comparing apexes fitted from the published composite against the
focal track (`observation088-5`, keypoint 0) gives:

| Callout | Anchor frame | Distance to fitted apex | 0.5 s window circular-mean bearing |
| --- | --- | --- | --- |
| pre (orange) | 17,623 | ~7.8 m | 93.2° (instantaneous 102.3°) |
| post (purple) | 80,780 | ~5 m (path bend) | 214.1° (instantaneous 215.3°) |

The paper's hand-drawn cone lines deviate 2.5-4° from the exact
tangents from the anchor positions to the inset disc; the generator
draws the exact tangents. The inset arrows follow the visible animals in
the rotated photographs.
The provenance records the measured compass bearings separately.

The 09.32.58 drone crop is the orange pre-stimulus inset. Its 208.78 pt
square rotates 90° counterclockwise about (115.04, 331.45) pt, behind a
65.25 pt radius mask. The 09.35.58 crop is the purple post-stimulus inset.
Its 214.90 pt square rotates 45° counterclockwise about (115.04, 743.02)
pt, behind a 67.16 pt radius mask. Both JPEGs retain the author's original
160 × 160 pixels. The rings, arrows and cone lines remain vector.

## Baseline numbers (from the pinned fit)

Directions are reported in degrees on the figure's polar-axis scale
(zero at the top, clockwise) -- the scale the panel c–e loops and
markers are drawn on. The house notebook's figure uses the same axis
(`set_theta_direction(-1)` over the model's angles) and the paper's
45°-step tick labels run clockwise from north, so the on-axis readings
below are what the published panels display. An on-axis reading is not
the geographic compass bearing of the fitted direction (that would be
(90° − model angle) mod 360°); the trajectory arrows in panel a use
true compass bearings from the per-frame data.

- Group mean direction (on-axis): **104.5° pre → 232.0° post**
  (posterior mean of the equally weighted per-track means).
- Group alignment (resultant length): **0.934 [0.909, 0.955] pre**,
  **0.912 [0.876, 0.944] post**.
- Observed headings in the fitted table: **3,222 pre / 2,195 post** of
  5,417.
- Paired per-animal turn (circular mean, on-axis, of the per-track
  post-minus-pre observed-mean directions): **130.9°**.
- Inset anchor bearings (true compass, from the bearing data):
  **93.2° pre (frame 17,623), 214.1° post (frame 80,780)**.
- Per-track values (empirical mean directions, concentration means and
  HDIs for all 18 animals) are in the provenance JSON.

## Tests

`tests/test_zebra_publication_figure.py` is hermetic (no data download,
no licensed fonts): synthetic-fit statistics cross-checked against
direct circular computations and 2π-wrap invariance, a double-render
SHA-256 determinism + no-`<text>` check of the vector export with DejaVu
Sans, the `prepare_inputs` tamper/refusal behavior, the committed-asset
manifest, a pixel-mismatch rejection test, and an end-to-end comparison
of the final photo and inset JPEG streams, signed placement matrices and
rendered image patches against the author's reference.

`tests/test_zebra_map.py` covers the `draw_escape_map` extraction itself
(the refactor keeps `plot_escape_map` a thin wrapper and
`_finish_escape_map` gained `scale_bar`/`compass`/`legend_loc`/
`tight_layout`/`legend_anchor`/`fancybox`/`framealpha` switches with the
standalone defaults preserved; all 16 pre-existing tests pass).

CI (`.github/workflows/zebra-figure.yml`) runs the hermetic tests and a
deterministic-render job that downloads the pinned data, renders the
figure twice with Liberation Sans, and compares both output directories; the CI
render is not byte-compared against the committed Arial exports.
The PDF export removes the Cairo creation timestamp. The provenance JSON
records the SHA-256 of each graphical export.
