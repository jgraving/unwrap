# Bat composite layout review

The landscape composite pairs the tagged-bat photo and treatment icons with
all 67 flight trajectories. The bottom row shows absolute step-level PSIS-LOO
ELPD with one SE per model, posterior directedness over remaining path distance,
and full vMU density at 0° relative to control.

The goal-window figure shows the full vMU predictive probability within ±45°
of the goal for each of the four arms. Points and whiskers show posterior
means and 95% HDIs. The separate eye-cover sensitivity plot compares the
eye-cover-only arm with control across angular windows. The figures evaluate
the fitted group parameters with individual offsets set to zero. The JSON
records draw-wise contrasts and the pooled contrast across the three groups
with eye cover, weighted by their modeled bat counts (14, 15, and 9).

From the repository root, install the figure dependencies and render from the
pinned fit and data snapshots:

```bash
uv pip install -e '.[figure]'
python notebooks/goldshtein_bats_layout_review.py --download \
    --font-family Arial --font-dir /path/to/arial-fonts
```

Pass `--data-root` to use verified local inputs without downloading. The
renderer writes PDF and PNG previews and a JSON record to this directory by
default. Pass `--output-dir` for a separate destination, or `--photo` to use
the final tagged-bat photograph when it is available.
