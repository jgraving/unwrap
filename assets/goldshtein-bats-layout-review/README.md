# Bat homing composite

[PDF](bat_composite.pdf) · [PNG](bat_composite.png) · [Provenance](bat_layout_review.json)

## Figure caption

**Hierarchical von Mises–uniform mixture modeling of bat homing after sensory
manipulation.** Free-flying Kuhl’s pipistrelles (*Pipistrellus kuhlii*) were
translocated within their home range and tracked during homing. Headings of
60 m steps were measured relative to the instantaneous bearing to the goal,
with 0° toward the goal and ±180° directly away.

**(a)** Homing trajectories after translocation. Control bats carried a tracking
tag only (n = 29). Vision was occluded alone (n = 14), in combination with a
magnetic manipulation (n = 15), or in combination with both magnetic
manipulation and anosmia (n = 9). The cyan polygon is the lake (Agamon), one of
the two homing destinations; the star marks the roost (the other destination),
crosses mark release sites, and dashed segments indicate tracking gaps.
The tagged-bat photograph is reproduced from the original manuscript figure
with permission. Basemap © OpenStreetMap contributors (ODbL).

**(b)** The von Mises–uniform mixture (vMU) represents headings as a mixture
of a directed von Mises component and an undirected uniform component.
Predictive performance is compared with a single von Mises model (vM) using
stepwise PSIS–LOO expected log predictive density. Points show vMU minus vM
overall and within each treatment. Whiskers show one standard error.

**(c)** Polar curves show the posterior mean heading density under the full
mixture for each treatment, with pointwise 95% highest-density intervals.
At the goal (0°), population contrasts subtract the control density within
each posterior draw, setting individual offsets to zero. Each bat’s density
at the goal is estimated under its observed treatment and compared with the
uniform density, 1/(2π). Points and bars show posterior means and 95%
highest-density intervals.

**(d)** Directedness compares the directed and undirected heading likelihoods
under equal prior odds, with 0.5 indicating equal evidence. Step scores were
smoothed over up to ten consecutive usable steps, with greater weight on recent
steps. Maps and heatmaps show posterior mean scores. Each heatmap row is one bat,
aligned by path distance remaining. White cells indicate unavailable steps.
Treatment curves average bats within each posterior draw and show pointwise 90%
highest-density intervals where at least two bats contribute. The polar ring
shows how the score varies with heading under the control population fit. A
sustained goalward bout begins at the first of five consecutive usable steps
with bearing error less than 90° and smoothed directedness above 0.5. Solid dots
and bars mark its median onset and central 90% interval among draws with a
qualifying bout. Open circles mark bouts beginning at the first usable step.

## Rendering

From the repository root, install the figure dependencies and render from the
pinned fit and data snapshots:

```bash
uv pip install -e '.[figure,notebook]'
python notebooks/goldshtein_bats_layout_review.py --download \
    --font-family Arial --font-dir /path/to/arial-fonts
```

Pass `--data-root` to use verified local inputs without downloading. The
renderer writes the PDF, PNG, and provenance JSON to this directory by default.
Pass `--output-dir` to choose a separate destination and `--photo` to supply
the tagged-bat photograph.

## Basemaps

The eight flight maps (panels a and d) sit on a **vector basemap**
drawn from a pinned OpenStreetMap extract (water, roads, buildings —
clipped to the map extents from the Geofabrik israel-and-palestine
shapefile by `notebooks/build_bats_vector_basemap.py`, committed as
`notebooks/assets/goldshtein_bats_figure/basemap_osm.geojson.gz`).
Because the basemap is real vectors rather than map tiles, it stays
sharp at any zoom in the PDF and renders fully offline.  Basemap
features are © OpenStreetMap contributors (ODbL) — the render JSON
records the extract hash and this attribution; a figure credit
("Basemap © OpenStreetMap contributors") belongs in the caption.

The pinned lake outline is drawn on top of the OSM water in a cyan that the
basemap's pale water never matches, so the destination reads at a glance, and
panel a's key carries a matching colour box ("Lake"); its interior
rings render as holes, so islands inside the outline stay land.  Both colours
live once in `goldshtein_bats_plots.py` (`LAKE_FILL`, `LAKE_EDGE`), which the
notebook's own map cells use too.  If the vector extract is ever missing,
the renderer falls back to live CartoDB "Positron" raster tiles at a
pinned z16 (`CARTO_API_KEY` in the environment for clean tiles — never
stored in the repository — otherwise the anonymous pool answers with
watermarked "API key required" tiles), and finally to the plain map
face with a printed notice.  The data, fits and geometry behind every
panel remain the pinned inputs recorded in `bat_layout_review.json`.

The directedness scale uses Fabio Crameri’s Vanimo in reverse order, with
yellow-green for undirected evidence and purple-pink for directed evidence.
[Scientific Colour Maps](https://www.fabiocrameri.ch/colourpalettes/)
