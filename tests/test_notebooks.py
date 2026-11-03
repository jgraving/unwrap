"""Local-only execution tests for the repository notebooks.

Run every notebook end to end with its committed cell contents and
check that each one still renders its inline figures::

    UNWRAP_TEST_NOTEBOOKS=1 python -m pytest tests/test_notebooks.py -v

The flag is required on purpose: a plain ``pytest -q`` (the unit test
suite, and what the repository's test workflow runs) must never spend
hours fitting research models. These tests are skipped under
``GITHUB_ACTIONS`` as well - notebook execution is local by design,
not a CI coverage target.

Select a subset by notebook stem (comma separated)::

    UNWRAP_TEST_NOTEBOOKS=1 UNWRAP_NOTEBOOKS=circular_regression,beetles \
        python -m pytest tests/test_notebooks.py -v

Python notebooks run on the ``python3`` kernel and the R notebooks on
the ``ir`` kernel (``Rscript -e 'IRkernel::installspec()'`` once). A
notebook whose kernel is not installed is skipped, not failed. The R
budget cells read ``UNWRAP_R_QUICK``: set it to ``1`` for the short
two-chain smoke budget, and leave it unset for full sampling.
``UNWRAP_NOTEBOOK_TIMEOUT`` overrides the per-cell timeout in seconds
(default 4 hours).
"""

import os
import sys
from pathlib import Path

import pytest

nbclient = pytest.importorskip("nbclient")
nbformat = pytest.importorskip("nbformat")

REPO_ROOT = Path(__file__).resolve().parents[1]

_LOCAL_ONLY = "notebook execution is local-only (GITHUB_ACTIONS guard)"


def _notebook_cases():
    cases = []
    for directory, kernel in (("notebooks", "python3"), ("notebooks_r", "ir")):
        base = REPO_ROOT / directory
        if not base.is_dir():
            continue
        for path in sorted(base.glob("*.ipynb")):
            cases.append(pytest.param(path, kernel, id=f"{directory}/{path.name}"))
    return cases


def _selected(path):
    wanted = os.environ.get("UNWRAP_NOTEBOOKS", "").strip()
    if not wanted:
        return True
    stems = {part.strip() for part in wanted.split(",") if part.strip()}
    return path.name in stems or path.stem in stems


def _kernel_installed(kernel_name):
    from jupyter_client import KernelSpecManager

    try:
        KernelSpecManager().get_kernel_spec(kernel_name)
        return True
    except Exception:
        return False


_TIMEOUT = int(os.environ.get("UNWRAP_NOTEBOOK_TIMEOUT", "14400"))


@pytest.mark.parametrize("path,kernel_name", _notebook_cases())
def test_notebook_executes(path, kernel_name):
    if os.environ.get("GITHUB_ACTIONS"):
        pytest.skip(_LOCAL_ONLY)
    if not os.environ.get("UNWRAP_TEST_NOTEBOOKS"):
        pytest.skip(
            "set UNWRAP_TEST_NOTEBOOKS=1 to run notebook execution "
            f"(UNWRAP_NOTEBOOKS=stem,stem selects a subset); see {path.name}"
        )
    if not _selected(path):
        pytest.skip(f"not selected by UNWRAP_NOTEBOOKS ({path.name})")
    if not _kernel_installed(kernel_name):
        pytest.skip(f"kernel '{kernel_name}' not installed")
    notebook = nbformat.read(path, as_version=4)
    client = nbclient.NotebookClient(
        notebook,
        timeout=_TIMEOUT,
        kernel_name=kernel_name,
        # Run from the notebook's own directory: the notebooks resolve
        # the repository root (and their data) relative to it.
        resources={"metadata": {"path": str(path.parent.resolve())}},
    )
    client.execute()
    nbformat.validate(notebook)
    # Every committed notebook renders at least one inline figure.
    assert any(
        "image/png" in output.get("data", {})
        for cell in notebook.cells
        for output in cell.get("outputs", [])
    ), f"no inline figures in {path}"
