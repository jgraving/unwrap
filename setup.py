
from setuptools import setup, find_packages

# Read the requirements from the requirements.txt file
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='unwrap',
    version='0.0.1',
    author='Your Name',
    author_email='your.email@example.com',
    description='A brief description of your package',
    packages=find_packages(),
    install_requires=requirements,
    extras_require={
        "test": ["pytest>=8,<10"],
        "notebook": [
            "jupyterlab>=4,<5",
            "nbclient>=0.10,<1",
            "ipykernel>=6,<8",
            "seaborn>=0.13",
            "pyarrow>=14",
            "joblib>=1.3",
            "numba>=0.60",
            "tqdm>=4.66",
            "xarray>=2024.1.0",
            # bats (Goldshtein) dataset loaders: legacy .xls + geojson
            "python-calamine>=0.2",
            "geopandas>=1.0",
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)
