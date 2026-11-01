"""Data loading for notebooks/zebra_escape.ipynb.

The study fits a hierarchical von Mises GLM to block-averaged body
orientations of a focal group of Grevy's zebra recorded by an overhead drone
at the Mpala Conservancy. The preprocessed table (one row per 10 s window,
per track and condition) is published in the ``jgraving/unwrap-data``
repository; :func:`load_zebra_model_table` pulls it end to end through
:func:`unwrap.datasets.load_zebra_bearings` (the data repository is
shallow-cloned to ``~/.cache/unwrap-data`` on first use; afterwards the load
resolves to local files with no network access).

The block averaging was done once when the table was published and is
documented here for provenance: body orientations are circularly averaged per
track and condition over non-overlapping 300-frame windows (10 s at 30 Hz),

    theta = atan2(mean_t(sin b_t), mean_t(cos b_t)),

giving approximately independent observations; each window is labelled
pre-approach (baseline) or post-approach (escape response); and tracks
appearing in only one condition are dropped. The raw source is still
available frame by frame through ``unwrap.datasets.load_zebra_bearings_raw``
(4 zebras, 59 tracks, 3,897,536 frames) if the averaging ever needs to be
reproduced.
"""

import pandas as pd

from unwrap.datasets import load_zebra_bearings


def load_zebra_model_table():
    """Load the block-averaged bearings table (end to end).

    Returns
    -------
    pandas.DataFrame
        Columns ``Track``, ``Observation``, ``Condition`` (ordered
        categorical pre-approach / post-approach, exactly as the notebook's
        models expect), ``HeadingRadians`` (in $(-\\pi, \\pi]$) and
        ``PostStimulus`` (0/1 flag).
    """
    df = load_zebra_bearings().copy()
    df.reset_index(drop=True, inplace=True)
    return df
