# Bat homing composite

The main-text composite pairs the tagged-bat photo placeholder and treatment
icons (a) with all 67 flight trajectories colored by posterior directedness
(b). The bottom row shows absolute step-level PSIS-LOO ELPD with one SE per
model (c), posterior directedness over remaining path distance (d), and
changes in predictive density at the goal bearing from control (e).

Panel e evaluates the full von Mises–Uniform (vMU) predictive heading density
at 0° relative to the local goal bearing, using fitted group parameters with
individual offsets set to zero. The diagram marks the height of an illustrative
density at 0°. The plot subtracts the control density from each treatment
density within the same posterior draw. Points and whiskers show posterior
means and 95% highest-density intervals. The dashed line at zero marks the
control reference.

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
