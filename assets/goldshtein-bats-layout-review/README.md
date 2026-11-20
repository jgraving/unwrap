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
manipulation and anosmia (n = 9). Crosses mark release sites, the star marks
the roost, and dashed segments indicate tracking gaps. Tagged-bat photograph
by [Photographer]. © [Year] [Copyright holder].

**(b)** The von Mises–uniform mixture (vMU) represents headings as a mixture
of a directed von Mises component and an undirected uniform component.
Predictive performance is compared with a single von Mises model (vM) using
stepwise PSIS–LOO expected log predictive density. Points show vMU minus vM
overall and within each treatment. Whiskers show one standard error.

**(c)** Population-level posterior predictive heading densities under the full
mixture. Curves and bands show posterior means and pointwise 95% highest-density
intervals. The lower plot shows treatment-minus-control differences in density
at 0° (rad⁻¹), calculated within each posterior draw. Points and bars show
posterior means and 95% highest-density intervals. The dashed line marks zero
difference.

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

The directedness scale uses Fabio Crameri’s Vanimo in reverse order, with
yellow-green for undirected evidence and purple-pink for directed evidence.
[Scientific Colour Maps](https://www.fabiocrameri.ch/colourpalettes/)
