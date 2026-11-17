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
