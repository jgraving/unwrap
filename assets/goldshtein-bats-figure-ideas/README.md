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

In panel a, faded icons identify the senses manipulated in each group.

Panel b shows the geography of the 67 routes, including release sites, the
roost, the lake and tracking gaps. Shades of each treatment color distinguish
individual routes. Panel d aligns the corresponding treatment
heatmaps below the maps. Each row is one bat, and all four tiles show the
last 6.6 km of modeled path on the same distance and directedness scales.
White cells mark distance beyond an individual's modeled path. The
group-average curves and their uncertainty are candidates for a supplement.
The blue Control point in panel e marks the zero difference from control.

Render from the same verified tracks and model fits as the main composite.

```bash
uv pip install -e '.[figure,notebook]'
python notebooks/goldshtein_bats_figure_ideas.py --download \
    --font-family Arial --font-dir /path/to/arial-fonts
```

## Panel b alternative

[Panel PDF](bat_model_panel.pdf) · [PNG preview](bat_model_panel.png)

The map puts all 67 routes on the same geographic scale, colored by
treatment. Crosses mark release sites and the star marks the roost.
The plots compare the fitted von Mises and von Mises–Uniform heading
distributions with the observed distribution of absolute heading error.
Each posterior draw averages 15° bin probabilities over the 4,942 modeled
steps, using the fitted effects for each bat. Bands are pointwise 95%
highest-density intervals for the fitted bin probabilities. The lower
plot subtracts the observed proportion in each bin.

Within 30° of the goal bearing, the observed share is 0.413. The fitted
shares are 0.381 for von Mises and 0.412 for von Mises–Uniform.

```bash
python notebooks/goldshtein_bats_model_panel.py --download \
    --font-family Arial --font-dir /path/to/arial-fonts
```

## Compact-map composite and model diagnostics

[Composite PDF](bat_composite_methods.pdf) · [PNG preview](bat_composite_methods.png) ·
[Mixture parameter plot](bat_mixture_parameters.pdf)

The four maps use the same geographic bounds and show each treatment on a
smaller scale. Panel c pairs absolute PSIS-LOO ELPD with the mean paired
change in ELPD per step (von Mises–Uniform minus von Mises) in six 30° bins
of observed absolute heading error. Positive values favor the mixture. The
largest total gain comes from steps within 30° of the goal bearing. The
largest gain per step comes from the 150°–180° bin. Panel d gives the
individual heatmaps the full width of the composite.

The parameter plot shows one pair of posterior means per bat. The
horizontal axis is the fitted directed-component weight. The vertical
axis is the fitted probability of a heading falling within ±30° of the
directed component’s own mean, conditional on that component.

```bash
python notebooks/goldshtein_bats_methods_review.py --download \
    --font-family Arial --font-dir /path/to/arial-fonts
```
