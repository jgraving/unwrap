# Development

Keep each PR focused: separate documentation, numerical changes, reorganizations,
and analysis updates. Include the reason for the change and how it was checked.

## Tests

```bash
uv venv --python 3.12
uv pip install -e ".[test]"
uv run --no-project python -m pytest -q
```

The test workflow covers Python 3.12–3.14 on Linux, Python 3.12 on macOS/Windows,
lowest-compatible direct dependencies on Python 3.12, and source/wheel installs.
The test suite checks numerical properties and backend agreement as well as a
short Bambi fit. That short fit tests execution, not convergence or parameter
recovery. Exact resolved dependencies are recorded in CI logs.

## Notebook execution (local only)

No GitHub workflow executes a notebook: the fits are too long for CI and
the results need a real machine anyway. `tests/test_notebooks.py`
executes every committed notebook (Python and R) end to end and checks
that each one still renders its inline figures. It is skipped under
`GITHUB_ACTIONS` and additionally requires a local opt-in flag, so the
unit suite above never triggers it:

```bash
uv pip install -e ".[notebook]"  # nbclient/nbformat + notebook stack
# everything (hours; research notebooks download data and fit full models):
UNWRAP_TEST_NOTEBOOKS=1 uv run --no-project python -m pytest tests/test_notebooks.py -v
# a subset, by notebook stem:
UNWRAP_TEST_NOTEBOOKS=1 UNWRAP_NOTEBOOKS=circular_regression,beetles \
    uv run --no-project python -m pytest tests/test_notebooks.py -v
```

The R notebooks need the `ir` kernelspec (see the
[R setup](getting-started-r.md)); set `UNWRAP_R_QUICK=1` for the short
two-chain budget while iterating. Before a change that touches a
notebook is merged, that notebook is rerun end to end from the checkout
this way, and the outputs stored in the repository come from that run.

## Build and preview documentation

Docs dependencies are separate from the modeling environment:

```bash
uv venv .venv-docs --python 3.12
uv pip install --python .venv-docs -r requirements-docs.txt
uv run --no-project --python .venv-docs python sync_docs_notebooks.py
uv run --no-project --python .venv-docs mkdocs build --strict
uv run --no-project --python .venv-docs mkdocs serve
```

The study notebooks render on the site from their committed outputs only -
the docs build never executes a notebook, and the files the build reads are
plain copies made by `sync_docs_notebooks.py` (gitignored; the notebooks in
`notebooks/` are the source of truth).

The documentation workflow builds with strict link validation and uploads the
static site as a GitHub Actions artifact for every PR. Pushes to `main` also
publish the site to <https://jgraving.github.io/unwrap/> through GitHub Pages.
Deployment only runs after a successful build, using the `github-pages`
environment. PRs do not deploy or receive deployment permissions.

Repository **Settings → Pages → Build and deployment → Source** must be set to
**GitHub Actions**. To republish without a code change, run the **Documentation**
workflow manually from the Actions tab with `main` selected.
