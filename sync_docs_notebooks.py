"""Copy the study notebooks into the docs tree for the documentation build.

The mkdocs site renders every notebook under ``notebooks/`` (Python) and
``notebooks_r/`` (R) inline (``mkdocs-jupyter``), so the build needs the
files under ``docs/``. They are copied - not moved or edited - to keep the
notebooks under their source folders the single source of truth. The
copies are gitignored and are rendered from their committed outputs only
(the docs build never executes a notebook).

Run before ``mkdocs build`` / ``mkdocs serve`` (the documentation
workflow does this for you).
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent

# notebook source folder (relative to the repository root) -> destination
# inside docs/, mirroring the folder name so rendered URLs stay stable
SOURCES = {
    "notebooks": "docs/notebooks",
    "notebooks_r": "docs/notebooks_r",
}


def sync_folder(src_dir: Path, dest_dir: Path) -> int:
    dest_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for src in sorted(src_dir.glob("*.ipynb")):
        target = dest_dir / src.name
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
    return count


def main():
    total = 0
    for src_name, dest_name in SOURCES.items():
        count = sync_folder(ROOT / src_name, ROOT / dest_name)
        total += count
        print(f"{dest_name}: {count} notebooks in sync")
    if total == 0:
        raise SystemExit("no notebooks found to sync")


if __name__ == "__main__":
    main()
