"""
Preprocessing for the dung-beetle security-light notebook.

Turns the raw security-light table (``load_beetle_security_light``) into the
model-ready frame the notebook's Bambi cells consumed: a binary ``Light``
indicator (1 = 'ON'), ``Heading`` in radians wrapped to ``(-pi, pi]`` via
``shifted_modulo``, and a categorical ``Individual`` column. Only the relevant
columns are kept and rows with missing headings dropped. Deterministic.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from unwrap.datasets import load_beetle_security_light
from unwrap.link.numpy import shifted_modulo


def build_beetle_table(df):
    """Apply the notebook's data cleanup to ``df`` and return the model frame."""
    # Binary indicator for the Light (1 for 'ON', 0 for 'OFF')
    df = df.copy()
    df["Light"] = df["Experiment"].apply(lambda x: 1 if "ON" in x else 0)

    # Heading: degrees -> radians
    df["Heading"] = pd.to_numeric(df["Heading"], errors="coerce")
    df["HeadingRadians"] = np.deg2rad(df["Heading"])

    # Wrap to (-pi, pi]
    df["HeadingRadians"] = df["HeadingRadians"].apply(
        lambda x: shifted_modulo(x) if pd.notnull(x) else x
    )

    # Keep the relevant columns, drop missing, reset index
    df = df[["Light", "Beetle", "Trial", "HeadingRadians"]].dropna().copy()
    df.reset_index(drop=True, inplace=True)

    # Rename Beetle -> Individual and make it categorical
    df = df.rename(columns={"Beetle": "Individual"})
    df["Individual"] = df["Individual"].astype("category")
    return df


def load_beetle_model_table():
    """Convenience: load the raw table and run the transformation."""
    return build_beetle_table(load_beetle_security_light())
