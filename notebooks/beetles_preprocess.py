"""Data loading and preprocessing for notebooks/beetles.ipynb.

The study fits von Mises hierarchy models to beetle heading angles from the
2019 security-light experiment. The raw per-trial table lives in the
``jgraving/unwrap-data`` repository; :func:`load_beetle_model_table` pulls
it via :func:`unwrap.datasets.load_beetle_security_light` (first use
shallow-clones the data repo to ``~/.cache/unwrap-data``) and applies the
notebook's original transforms: a binary light indicator, a
degrees-to-radians conversion of ``Heading``, wrapping to $(-\\pi, \\pi]$
via :func:`unwrap.link.numpy.shifted_modulo`, and the column rename /
categorisation the models expect.
"""

import numpy as np
import pandas as pd

from unwrap.datasets import load_beetle_security_light
from unwrap.link.numpy import shifted_modulo


def load_beetle_model_table():
    """Load the security-light trials and build the modelling table.

    Returns
    -------
    pandas.DataFrame
        Columns ``Light`` (0/1), ``Individual`` (category of ``Beetle``
        labels), ``Trial``, ``HeadingRadians`` (in $(-\\pi, \\pi]$); rows
        with non-numeric or missing headings are dropped. Bit-identical to
        the notebook's former inline data cell.
    """
    df = load_beetle_security_light()

    # Create a binary indicator for the Light (1 for 'ON', 0 for 'OFF')
    df['Light'] = df['Experiment'].apply(lambda x: 1 if 'ON' in x else 0)

    # Convert the 'Heading' column from degrees to radians
    df['Heading'] = pd.to_numeric(df['Heading'], errors='coerce')  # Convert non-numeric values (e.g., 'NA') to NaN
    df['HeadingRadians'] = np.deg2rad(df['Heading'])

    # Normalize the heading radians to be between -pi and pi
    df['HeadingRadians'] = df['HeadingRadians'].apply(lambda x: shifted_modulo(x) if pd.notnull(x) else x)

    # Keep only the relevant columns and drop rows with missing values
    df = df[['Light', 'Beetle', 'Trial', 'HeadingRadians']].dropna().copy()
    df.reset_index(drop=True, inplace=True)

    # Rename columns: change 'Beetle' to 'Individual'
    df = df.rename(columns={'Beetle': 'Individual'})
    df['Individual'] = df['Individual'].astype('category')

    return df
