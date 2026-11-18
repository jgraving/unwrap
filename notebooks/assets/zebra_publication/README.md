# Zebra publication figure inputs

The author-supplied Keynote PDF defines the page geometry, image placements
and circular clips. The reduced reference keeps its text, vectors and clips.
Only the panel-b photo and aerial mosaic image streams were replaced.

| Working file | Source and preparation |
| --- | --- |
| `reference.pdf` | 900 × 1080 pt author composite. The mosaic retains its 640 × 1940 pixels at JPEG quality 92 with 4:4:4 chroma. Panel b uses the working photo below. |
| `zebra_photo.jpg` | `2009_Jul30_IMG_5242.JPG`, resized from 2048 × 1365 to 1600 × 1066 with Lanczos, JPEG quality 92, 4:4:4 chroma and its ICC profile. At the published panel size this gives 238 ppi. |
| `drone_inset_pre.jpg` | Original 160 × 160 JPEG from 09.32.58. The orange inset rotates it 90° counterclockwise. |
| `drone_inset_post.jpg` | Original 160 × 160 JPEG from 09.35.58. The purple inset rotates it 45° counterclockwise. |
| `images.json` | SHA-256 hashes for the author files and the committed working files. The generator verifies the working files before drawing. |

The exported map uses the georeferenced tiles, trajectories and bearings
from `jgraving/unwrap-data`. The statistics use the saved fit at the same
pinned revision. See [the reproduction notes](../../ZEBRA_FIGURE.md).
