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

**B. Directedness by individual.** The 67 rows show the same posterior
directedness measure and trailing smoother as the main composite. Rows are
sorted by modeled path length within treatment. Zero is the last modeled step
of each flight. The horizontal scale uses 60 m per modeled step. White cells
mark distance beyond each flight's modeled path, and white separators mark
individual rows. The purple and olive-green ramps use even Oklab steps from black
at 0.5, with matched lightness and chroma across the midpoint. At 0.5 the
directed and undirected components have equal support. Values toward 0 favor
the undirected component and values toward 1 favor the directed component.
Smoothing resets at tracking gaps.

## Composite with flight heatmaps

[Composite PDF](bat_composite_heatmaps.pdf) · [PNG preview](bat_composite_heatmaps.png)

Panel b shows the geography of the 67 routes, including release sites, the
roost, the lake and tracking gaps. Shades of each treatment color distinguish
individual routes. Panel d aligns the corresponding treatment
heatmaps below the maps. Each row is one bat, and all four tiles show the
last 6.6 km of modeled path on the same distance and directedness scales.
White cells mark distance beyond an individual's modeled path. The
group-average curves and their uncertainty are candidates for a supplement.

Panel e shows full-mixture polar densities and two sets of paired contrasts
at the 0° goal bearing. The upper contrast combines directed and undirected
densities using the fitted mixture weight. The lower contrast uses the
conditional density of the directed component. Both use the same posterior
draws, with individual offsets set to zero. The blue Control point marks
zero in each plot.

| Group vs. control | Full mixture | Directed component | Full minus directed |
| --- | ---: | ---: | ---: |
| Eye cover | −0.132 [−0.267, −0.008] | −0.123 [−0.250, 0.019] | −0.009 [−0.077, 0.045] |
| + Magnet | −0.157 [−0.285, −0.029] | −0.150 [−0.286, −0.022] | −0.007 [−0.065, 0.048] |
| + Anosmia | −0.169 [−0.308, −0.014] | −0.192 [−0.344, −0.034] | 0.024 [−0.029, 0.070] |

Entries are posterior means and 95% highest-density intervals in inverse
radians. Differences are calculated within each draw before summarizing.

Within each draw, multiplying the directed-density contrast by the two
groups' mean mixture weight gives its contribution to the full contrast.
Multiplying the weight contrast by the two groups' mean directed density
above the undirected baseline gives the remaining contribution. The
posterior means of the weight contribution are −0.022, −0.022, and 0.010
inverse radians in table order. Each 95% interval spans zero.

Render from the same verified tracks and model fits as the main composite.

```bash
uv pip install -e '.[figure,notebook]'
python notebooks/goldshtein_bats_figure_ideas.py --download \
    --font-family Arial --font-dir /path/to/arial-fonts
```
