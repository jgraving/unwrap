
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
        # parquet engine needed by the dataset loader tests; >=16 is the
        # first line compatible with numpy 2, which the minimum-deps CI
        # job resolves
        "test": ["pytest>=8,<10", "pyarrow>=16"],
        # publication-figure generators (zebra/locust/beetles): SVG->PDF
        # raster export, PDF image-object restore, and netcdf fit caches
        "figure": [
            "cairosvg>=2.7,<3",
            "pyarrow>=16",
            "seaborn>=0.13",
            "xarray>=2024.1.0",
            "h5netcdf>=1.3",
            "h5py>=3.10",
            "pymupdf>=1.26.6,<2",
        ],
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
            # bats (Goldshtein) dataset loaders + notebook: legacy .xls, geojson, map panels
            "contextily>=1.6",
            "geopandas>=1.0",
            "shapely>=2.0",
            "pyproj>=3.6",
            "python-calamine>=0.3",
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',
)
