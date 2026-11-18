# circstats_graphical_abstract assets

Pinned inputs for `notebooks/circstats_graphical_abstract.py`.

| File | SHA-256 | Description |
| --- | --- | --- |
| `artwork.json` | `d2d6991cb522ac180dea76bf130511049a3562cded5b2b501cb822685a9fd612` | Every vector drawing of the reference PDF (526 paths, page points, y down from the top-left) plus the 236 measured text spans, the 12 locust-unit grids, the pinned/skeleton region split and the four diagonal circle outlines. |
| `reference.pdf` | `5874a6da22930ff3b0e53f37c3ffd0074c2e88ddb25cdd4b909f7187a9e4e55c` | The ground-truth target: `figures/Graphical Abstract.pdf` from `jgraving/circstats` after its manual Inkscape polish (599 x 342 pt, 100 % vector). |

## Provenance of `artwork.json`

`artwork.json` was extracted from `reference.pdf` with PyMuPDF
(`page.get_drawings()` + `page.get_text("dict")`), so every coordinate,
colour, stroke width and dash pattern is the exact value of the published
figure.  The 476-pinned / 50-skeleton split records which paths are
hand-made or stochastic artwork (replayed verbatim) and which are the
deterministic skeleton (drawn from the same measured constants).

Re-extract with (from the circstats checkout):

```
python - <<'PY'
import json, pymupdf
p = pymupdf.open("figures/Graphical Abstract.pdf")[0]
print(len(p.get_drawings()), "drawings")
PY
```

The committed hash is verified by
`tests/test_circstats_graphical_abstract.py`.
