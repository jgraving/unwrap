# Development

Keep each PR focused: separate documentation, numerical changes, reorganizations,
and analysis updates. Include the reason for the change and how it was checked.

## Tests

```bash
uv venv --python 3.12
uv pip install -e ".[test]"
uv run --no-project python -m pytest -q
```

The test workflow covers Python 3.10–3.14 on Linux, Python 3.12 on macOS/Windows,
lowest-compatible direct dependencies on Python 3.10, and source/wheel installs.
The test suite checks numerical properties and backend agreement as well as a
short Bambi fit. That short fit tests execution, not convergence or parameter
recovery. Exact resolved dependencies are recorded in CI logs.

The separate notebook workflow executes the introductory notebook and uploads
its executed copy for review. Full research analyses and R/MATLAB are not yet CI
coverage targets.

## Build and preview documentation

Docs dependencies are separate from the modelling environment:

```bash
uv venv .venv-docs --python 3.12
uv pip install --python .venv-docs -r requirements-docs.txt
uv run --no-project --python .venv-docs mkdocs build --strict
uv run --no-project --python .venv-docs mkdocs serve
```

The documentation workflow builds with strict link validation and uploads the
static site as a GitHub Actions artifact for every PR. Pushes to `main` also
publish the site to <https://jgraving.github.io/unwrap/> through GitHub Pages.
Deployment only runs after a successful build, using the `github-pages`
environment. PRs do not deploy or receive deployment permissions.

Repository **Settings → Pages → Build and deployment → Source** must be set to
**GitHub Actions**. To republish without a code change, run the **Documentation**
workflow manually from the Actions tab with `main` selected.
