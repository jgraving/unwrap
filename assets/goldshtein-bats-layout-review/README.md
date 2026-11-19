# Bat homing composite

The main-text composite pairs the tagged-bat photo placeholder and treatment
icons (a) with all 67 flight trajectories colored by posterior directedness
(b). The bottom row shows absolute step-level PSIS-LOO ELPD with one SE per
model (c), posterior directedness over remaining path distance (d), and
changes in goal alignment from control (e).

Panel e integrates the full von Mises–Uniform (vMU) predictive heading density
within ±45° of the goal, using fitted group parameters with individual offsets
set to zero. It subtracts the control probability from each treatment
probability within the same posterior draw. Points and whiskers show posterior
means and 95% highest-density intervals. The dashed line at zero marks the
control reference. The JSON also records the pooled contrast across the three
groups with eye cover, weighted by their modeled bat counts (14, 15, and 9).

From the repository root, install the figure dependencies and render from the
pinned fit and data snapshots:

```bash
uv pip install -e '.[figure]'
python notebooks/goldshtein_bats_layout_review.py --download \
    --font-family Arial --font-dir /path/to/arial-fonts
```

Pass `--data-root` to use verified local inputs without downloading. The
renderer writes one composite as PDF and PNG and a JSON record to this directory by
default. Pass `--output-dir` for a separate destination, or `--photo` to use
the final tagged-bat photograph when it is available.
