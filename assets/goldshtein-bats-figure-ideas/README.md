# Candidate bat figure panels

[Goal-cone and flight panels PDF](bat_figure_ideas.pdf) · [PNG preview](bat_figure_ideas.png)

**A. Heading within the goal cone.** The upper plot shows the probability of
a heading falling within ±α of the local goal bearing. The inset marks α
from the goal bearing to one edge of the cone. The lower plot shows paired
treatment-minus-control differences on the probability scale. Curves are
posterior means and bands are pointwise 95% highest-density intervals.
At α = 45°, the fitted probabilities are 0.83 for control and 0.72–0.74
for the eye-cover groups. Predictions use the full vMU mixture with
individual offsets set to zero.

**B. Homing, flight by flight.** The 67 rows show the same posterior
directedness measure and trailing smoother as the main composite. Rows are
sorted by modeled path length within treatment. Zero is the last modeled step
of each flight. The horizontal scale uses 60 m per modeled step. Gray cells
mark distance beyond each flight's modeled path. Smoothing resets at tracking
gaps.

## Composite with flight heatmaps

[Composite PDF](bat_composite_heatmaps.pdf) · [PNG preview](bat_composite_heatmaps.png)

Panel b shows the geography of the 67 routes, including release sites, the
roost, the lake and tracking gaps. Panel d aligns the corresponding treatment
heatmaps below the maps. Each row is one bat, and all four tiles show the
last 6.6 km of modeled path on the same distance and directedness scales.
Gray cells mark distance beyond an individual's modeled path. The
group-average curves and their uncertainty are candidates for a supplement.

Render from the same verified tracks and model fits as the main composite.

```bash
uv pip install -e '.[figure,notebook]'
python notebooks/goldshtein_bats_figure_ideas.py --download \
    --font-family Arial --font-dir /path/to/arial-fonts
```
