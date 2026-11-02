"""Copy the study notebooks into the docs tree for the documentation build.

The mkdocs site renders every notebook in ``notebooks/`` inline
(``mkdocs-jupyter``), so the build needs the files under ``docs/``.
They are copied - not moved or edited - to keep the notebooks under
``notebooks/`` the single source of truth. The copies are gitignored
and are rendered from their committed outputs only (the docs build
never executes a notebook).

Run before ``mkdocs build`` / ``mkdocs serve`` (the documentation
workflow does this for you).
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    dest = ROOT / "docs" / "notebooks"
    dest.mkdir(parents=True, exist_ok=True)
    count = 0
    for src in sorted((ROOT / "notebooks").glob("*.ipynb")):
        target = dest / src.name
        same = (
            target.exists()
            and not target.is_symlink()
            and target.read_bytes() == src.read_bytes()
        )
        if not same:
            if target.exists() or target.is_symlink():
                target.unlink()
            target.write_bytes(src.read_bytes())
        count += 1
    print(f"docs/notebooks: {count} notebooks in sync")


if __name__ == "__main__":
    main()
