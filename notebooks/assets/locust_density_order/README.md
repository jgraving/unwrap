# Locust density–order source assets

Seven native JPEG image streams from the author's manuscript figure,
`jgraving/circstats/figures/locust_density_order_figure.pdf` at revision
`9bef710eb44fd8d20a61a7119a12e4933d02ce6d`. These are source pixels,
not rendered PDF page crops or AI-generated replacements.

| File | PDF object | Source pixels | Use |
| --- | ---: | --- | --- |
| `real_locust.jpeg` | 12 | 784 × 441 | Real focal locust |
| `virtual_locust_front.jpeg` | 14 | 311 × 175 | Virtual agent, anterior view |
| `virtual_locust_lateral.jpeg` | 13 | 419 × 235 | Virtual agent, lateral view |
| `virtual_locust_rear.jpeg` | 15 | 310 × 174 | Virtual agent, posterior view |
| `vr_arena.jpeg` | 11 | 803 × 346 | Panoramic VR apparatus |
| `real_swarm.jpeg` | 20 | 490 × 735 | Real swarm |
| `virtual_swarm.jpeg` | 21 | 445 × 445 | Rendered virtual swarm |

`images.json` records the source PDF SHA-256 and each extracted image hash.
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
