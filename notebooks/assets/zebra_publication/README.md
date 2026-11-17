# Zebra publication figure inputs

Raster and reference inputs for
`notebooks/zebra_publication_figure.py`, all hash-pinned in `images.json`
(the generator refuses to run if any file is missing or altered).

| File | Role | Notes |
|---|---|---|
| `reference.pdf` | The published composite figure (900 x 1080 pt), the layout the generator reproduces | Provided by the authors; Keynote/Quartz production. Layout constants in the generator are measured from this file. |
| `drone_inset_pre.jpg` | Drone crop of the focal zebra before the approach (callout inset a, orange) | 160 x 160 px, author-provided crop; 160 px is marginal at 300 dpi, so the inset is intentionally small in the layout. |
| `drone_inset_post.jpg` | Drone crop of the focal zebra after the approach (callout inset a, purple) | 160 x 160 px, author-provided crop. |
| `zebra_photo.jpg` | Species photograph for panel b | 2000 x 1333 px, author-provided. |
| `images.json` | SHA-256 manifest for the files above plus `reference.pdf` | Checked by `zebra_publication_figure.check_assets()`. |

The insets are raw crops: the circular ring, the callout cones, and the
body-heading arrows are drawn by code on top of them. The heading arrows
are the circular-mean compass bearing over the +/- 0.5 s window centred on
the recovered anchor frame on the focal track (the paper's hand-drawn
arrows differ by up to ~20 degrees; see `notebooks/ZEBRA_FIGURE.md`).

Everything else the figure draws (mosaic, trajectories, bearings, fit)
comes from `jgraving/unwrap-data` at the pinned revision; see
`notebooks/ZEBRA_FIGURE.md` for the data and fit provenance.
