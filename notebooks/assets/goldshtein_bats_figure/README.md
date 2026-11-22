# Goldshtein bat figure materials

Shared materials for the bat homing figure pipeline:
`goldshtein_bats_figure_core.py` (pinned inputs, artwork loading, font
resolution) and `goldshtein_bats_layout_review.py` (the publication
composite renderer, written up in
[`assets/goldshtein-bats-layout-review/README.md`](../../assets/goldshtein-bats-layout-review/README.md)).

The composite draws everything it shows — including the eight flight maps,
which sit on live CartoDB "Positron" tiles — from code and pinned data.
What lives here is only what has no code counterpart:

* `artwork.json` — verbatim vector extraction of the authored original
  (schema 1: 459 vector drawings, 8 image objects with 12 placements,
  64 text spans). The renderer takes the panel-b treatment-sense icon
  grid (the 69 drawings whose bounding-box centre sits in the icon-grid
  band, checked on load) from it; all other reference content is
  superseded by code-drawn panels.
* `image_279.jpeg` — the tagged-bat photograph panel b sits on (page
  alpha 0.4684, a uniform mask the renderer synthesises; the "1.4cm"
  ruler is baked into the JPEG).
* `basemap_osm.geojson.gz` — the pinned vector basemap: water, roads and
  buildings clipped to the map extents from the Geofabrik
  israel-and-palestine free-shp extract (retrieved 2026-11-22), written
  by `notebooks/build_bats_vector_basemap.py`. The renderer draws this
  file directly — no map service or download at render time. Features
  are © OpenStreetMap contributors (ODbL), so the figure caption should
  credit "Basemap © OpenStreetMap contributors".
* `vanimo.txt` — Fabio Crameri's Vanimo colour table (256 stops) used
  for the directedness scale, so renders do not depend on the Scientific
  Colour Maps package being installed.
* `reference.pdf` — the 583 × 626 pt authored original
  (`goldshtein_bat_figure.pdf`, polished in Inkscape). Not used at render
  time; it is the layout reference the composite's page geometry was
  measured from.

The original figure's own basemap tile mosaics (xref JPEGs 23 and 45)
were committed here until the composite started drawing live CartoDB
tiles itself (see the renderer's `draw_track_map`); they are no longer
needed. The `image_*.jpeg` names are the PDF xrefs of the original image
objects, and `artwork.json` still lists all of them as extraction
provenance.

The renderer verifies every file above against its SHA-256
(`ASSET_SHA256` in `goldshtein_bats_figure_core.py`) before drawing and
never overwrites or re-downloads them.

## Pinned data and fits

Everything else the figure needs is pinned to
`jgraving/unwrap-data` revision `86c41b5f43eb61191405b9bc85baf730a492a06f`
with SHA-256 checksums in the core module's `INPUT_HASHES` (the two
tracking workbooks, the Agamon–Hula lake outline, and the two cached
fits `bats_vm_v1` / `bats_vmu_v1` with their metadata), and fetched on
demand with `--download`. The slim fit snapshots strip exactly the
per-observation log-likelihood (vMU also the per-observation
deterministics), which `build_context` recomputes from the stored draws
through the notebook cells' own restoration paths (the rebuilt bambi
graph for the baseline, `restore_per_obs_deterministics` plus the rebuilt
raw-PyMC graph for the mixture).
