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

The generator defines all display crops and positions in code. To repeat
extraction, see [the figure instructions](../../LOCUST_DENSITY_FIGURE.md).
