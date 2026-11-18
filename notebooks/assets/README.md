# Beetle artwork

`beetle_arena.svg` contains the original vector arena and beetles used in
panel a of the beetle figure. It retains the original colours and paths;
the panel text is supplied by `../beetles_figure.py`.

Source: `figures/BeetlesModel.svg` in the author's `jgraving/circstats`
manuscript repository, revision `9bef710eb44fd8d20a61a7119a12e4933d02ce6d`.
Only paths `path2` through `path42` were extracted, with a view box of
`20 8 262 122`. No manuscript text, fitted plot paths, or other figure
panels are included. The illustration is schematic, not measured tracks.

Keeping this small vector asset alongside the generator makes the complete
figure reproducible without access to the manuscript repository.

`beetle_floodlight.svg` restores the floodlight graphic from the same source
revision. It contains paths `path251` through `path258`, extracted from group
`g124` without its group translation or text. Its 33 × 33 view box is centred
on the graphic's circular lamp body (353.67169, 132.2807). The generator
positions that centre on the ±45° sector's central spoke after resolving
the polar axes' final aspect ratio, and supplies the label.
It is a schematic icon, not a measurement of the lamp's angular extent.

# Locust composite figure artwork

`locust_circ_stats/panel_a.json` holds the geometry of panel a of the published
locust circular-statistics figure, read from the figure PDF: the twelve
virtual-swarm image placements (CTM plus clip window each), the thirty-four
vector paths in the panel's region, and the fourteen text spans.
`locust_circ_stats/virtual_swarm.png` is the sprite those placements stamp,
cropped from the same page. `../locust_circ_stats_figure.py` replays them onto
a 900 x 600 pt stage.

The path list spans y = 0-321 pt: panel a's bottom swarm has two rows, and a
first extraction pass that stopped at y = 300 pt left the second row's two grey
flight arrows out.

The published PDF is not committed; only this panel's artwork and geometry are.
No fitted plot data from the published panels is taken from the PDF - panels b,
c and d are computed from the `unwrap-data` fits, and the layout constants in
the script record the numbers measured off the page so they can be checked.

The review snapshots in `../../assets/locust-circ-stats/` are the output of
`python ../locust_circ_stats_figure.py` (PDF, PNG, SVG, an HTML review page and
a provenance JSON), committed for the same reason as
`../assets/beetles-publication/`: so a reviewer can look at the rendered figure
without running the fit. `comparison_vs_published.png` and
`panel_comparison.png` place that output next to the published page for the same
check.
