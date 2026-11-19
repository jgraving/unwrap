# Candidate bat figure panels

**A. Heading within the goal cone.** The upper plot shows the predicted share
of headings within a cone centered on the local goal bearing, across cone
half-angles from 0° to 180°. The lower plot shows treatment-minus-control
differences within each posterior draw. Curves are posterior means and bands
are pointwise 95% highest-density intervals. The vertical guide marks ±45°,
where the predicted shares are 83% for control and 72–74% for the eye-cover
groups. Predictions use the full vMU mixture with individual offsets set to
zero.

**B. Homing, flight by flight.** The 67 rows show the same posterior
directedness measure and trailing smoother as the main composite. Rows are
sorted by modeled path length within treatment. Zero is the last modeled step
of each flight. The horizontal scale uses 60 m per modeled step. Gray cells
mark distance beyond each flight's modeled path. Smoothing resets at tracking
gaps.

Render from the same verified tracks and model fits as the main composite.

```bash
uv pip install -e '.[figure,notebook]'
python notebooks/goldshtein_bats_figure_ideas.py --download \
    --font-family Arial --font-dir /path/to/arial-fonts
```
