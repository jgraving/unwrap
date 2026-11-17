# Locust density–order source assets

Seven photo/rendering rasters for the author's manuscript figure. Each file
matches the image checked in at
`jgraving/circstats/figures/locust_density_order_figure.pdf` (revision
`9bef710eb44fd8d20a61a7119a12e4933d02ce6d`) in content and placement, but
carries the higher-resolution pixels of the author's Keynote export of the
same figure (`locust_density_order_figure.pdf`, SHA-256 recorded in
`images.json`). These are source pixels, not rendered PDF page crops or
AI-generated replacements. `tools/upgrade_locust_figure_images.py` performs
the swap: each raster is resized to at most 400 effective pixels per display
point (never above native), encoded at JPEG quality 88, and re-anchored in
the SVG with a per-image scale transform so every display box, clip path and
ancestor transform stays exactly as published.

| File | PDF object | Export pixels | Stored pixels | Use |
| --- | ---: | ---: | --- | --- |
| `real_locust.jpeg` | 12 | 3840 × 2160 | 2179 × 1226 | Real focal locust |
| `virtual_locust_front.jpeg` | 14 | 1920 × 1080 | 866 × 487 | Virtual agent, anterior view |
| `virtual_locust_lateral.jpeg` | 13 | 1920 × 1080 | 1165 × 653 | Virtual agent, lateral view |
| `virtual_locust_rear.jpeg` | 15 | 1920 × 1080 | 863 × 485 | Virtual agent, posterior view |
| `vr_arena.jpeg` | 11 | 3840 × 1655 | 2233 × 962 | Panoramic VR apparatus |
| `real_swarm.jpeg` | 20 | 720 × 1080 | 720 × 1080 | Real swarm (native cap) |
| `virtual_swarm.jpeg` | 21 | 3840 × 3840 | 1235 × 1235 | Rendered virtual swarm |

The lateral and rear exports carry a few boundary rows the manuscript crops
omitted; both are cropped to the published framing before resizing.
`images.json` records the low- and high-resolution source PDF SHA-256s, the
resize policy, the manuscript payload hashes, and each stored image hash.
`reference.pdf` remains the checksum-pinned original for typography and crop
comparison; `tools/refresh_locust_figure_artwork.py` re-extracts the artwork
from it and deterministically re-applies the high-resolution payloads.
`selection.json` is automatically parsed from the versioned Python notebook's
saved five-model PSIS-LOO table; it is not inferred from plotted coordinates.

Photo credit: Centre for the Advanced Study of Collective Behaviour. The
manuscript caption states used with permission; no additional redistribution
licence is inferred. Preserve the credit with the figure.

`original_panels.svg` preserves panels a–e at the original 750 × 525-point
page coordinates, including clipping paths, image transforms, the 10 cm
scale bar, lettering, and the complete vector schematic. Original f–h panel
letters and titles are also retained. Only the statistical drawings and the
old unit line were removed. Source preparation
checks that retained a–e pixels are unchanged, excluding the unit line.
The original scale-label backing and its transparency mask are retained.
`images.json` records the SVG hash and preservation result. The generator
rebuilds statistical panels f–h in code and inserts a correctly typeset unit.
The presentation pass edits a parsed copy of this asset to reduce panel
letters, replace f–h titles with bold versions of the same wording, and
lighten d/e borders. Image geometry, clipping, and the complete scale-bar
path are preserved; the source asset itself is not rewritten.

To repeat extraction, see [the figure instructions](../../LOCUST_DENSITY_FIGURE.md).
