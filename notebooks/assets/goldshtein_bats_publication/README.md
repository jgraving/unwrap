# Goldschtein bat figure inputs

Pinned inputs for
`notebooks/goldshtein_bats_publication_figure.py`. The published composite
(`goldshtein_bat_figure.pdf`, polished in Inkscape, committed below as
`reference.pdf`) is rebuilt *programmatically*: every notebook panel is
rendered by the notebook's own plotting code from the pinned data and the
pinned model fits, and embedded in the 583 × 626 pt page at the
Inkscape-measured scale and anchor (see
[the reproduction notes](../../GOLDSHTEIN_BATS_FIGURE.md)). The committed
artwork supplies only what has no notebook counterpart:

* `artwork.json` — the full verbatim extraction of the reference (schema
  1: 459 vector drawings, 8 image objects with 12 placements, 64 text
  spans), from which the generator takes the panel-b artwork (the
  treatment-sense icon grid: the 69 drawings whose bounding-box centre
  sits in the icon-grid band, checked on load) and the chrome text (the
  panel letters a.–h, "Treatment groups" and the "Placeholder" photo
  watermark). All other reference vector content is superseded by the
  code-drawn panels.
* `image_23.jpeg` / `image_45.jpeg` — the published figure's own CartoDB
  tile mosaics, placed under the panel-h (×4) and panel-a (×2) flight
  maps. Feeding them back through the map functions' additive
  `basemap_img` parameter reproduces the reference basemaps offline
  (identical tracking data → identical extents → identical stretch).
* `image_279.jpeg` — the tagged-bat photograph under panel b (page
  alpha 0.4684, a uniform SMask the generator synthesises; the
  reference stores no SMask for it).

The reference's other rasters (the four panel-h track overlays
`image_26/29/32/35` with their SMask channels `smask_284–287`, and the
viridis colorbar `image_36`) are superseded by the code-drawn vector
tracks and colorbar and are no longer committed.

The generator verifies every file against its SHA-256 before drawing and
never overwrites or re-downloads them.

## Files

| File | Role | SHA-256 |
| --- | --- | --- |
| `artwork.json` | Schema-1 asset: 459 vector drawings, 8 image objects (12 placements, each with `rect_ydown`, `transform`, `paths_before`, `fill_alpha`), 64 text spans (each with baseline `origin`, line `dir`, `flags` and, for the 60 drawn spans, a measured `fill_alpha`), and the extraction provenance note. | `9da0eb2639fea22d6af64820e1072524de6b159705b42cfe8d12c05e3d851730` |
| `reference.pdf` | The 583 × 626 pt authored figure (Inkscape-polished `goldshtein_bat_figure.pdf`). | `849d8f791e9c53ffbcf559f4fb5a839fb770d44709b63991245a29cf35ae264f` |
| `image_23.jpeg` | 161 × 409 panel-h tile mosaic, placed 4× (h subpanels 1–4, image_23 placements). | `f5f042faee6ad39359f3cf757db286072139733e5241f2b458050fdfb7f01d1c` |
| `image_279.jpeg` | 201 × 189 panel-b tagged-bat photograph (fill alpha 0.4684). The "1.4cm" ruler is baked into this JPEG. | `6ffabea4e91fc7e53c7a97a97cb651642c8f8c392b1e4add3ae2fc7fad269734` |
| `image_45.jpeg` | 175 × 444 panel-a tile mosaic, placed 2× (Home / Agamon panels). | `95c074d54e570b973f1b0c969d6e54d0f9ae56636a2f568664723919bf49975b` |

The `image_*.jpeg` names are the PDF xrefs of the original image objects
(23, 279, 45).

## Pinned data and fits

Everything else the figure needs is pinned to
`jgraving/unwrap-data` revision `86c41b5f43eb61191405b9bc85baf730a492a06f`
with SHA-256 checksums in the generator's `INPUT_HASHES` (the two
tracking workbooks, the Agamon–Hula lake outline, and the two cached
fits `bats_vm_v1` / `bats_vmu_v1` with their metadata), and fetched on
demand with `--download`. The slim fit snapshots strip exactly the
per-observation log-likelihood (vMU also the per-observation
deterministics), which the generator recomputes from the stored draws
through the notebook cells' own restoration paths (the rebuilt bambi
graph for the baseline, `restore_per_obs_deterministics` plus the rebuilt
raw-PyMC graph for the mixture).
