"""Re-run a study notebook end to end and bake outputs back in.

Usage: rerun_nb.py <path/to/notebook.ipynb>
Kernel cwd is the notebook's own directory (metadata.path) and the kernel
inherits PYTHONPATH=<worktree> so the in-tree `unwrap` (with the new
modules) wins over whatever the venv has installed.
"""
import os
import sys
import time

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
_WORKTREE = "/tmp/opencode/wt-unified"
os.environ["PYTHONPATH"] = (
    _WORKTREE + os.pathsep + os.environ["PYTHONPATH"]
    if os.environ.get("PYTHONPATH") else _WORKTREE
)

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager

PY = "/home/jake/.venvs/unwrap-env/bin/python"
NB_REPO = os.path.abspath(sys.argv[1])
NB_DIR = os.path.dirname(NB_REPO)


class VenvKM(KernelManager):
    """Run the kernel with the unwrap-venv python explicitly."""

    @property
    def argv(self):
        return [PY, "-m", "ipykernel_launcher", "-f", self.connection_file]


nb = nbformat.read(NB_REPO, as_version=4)
# Force the venv's python3 kernelspec: study notebooks may carry stale
# environment-specific kernel names (e.g. "locust") that are not installed.
# The VenvKM below supplies the launch command, so only the spec lookup
# needs to succeed.
nb.metadata["kernelspec"] = {
    "display_name": "Python 3 (ipykernel)", "language": "python",
    "name": "python3",
}
t0 = time.time()
NotebookClient(
    nb,
    timeout=16200,
    kernel_manager_class=VenvKM,
    resources={"metadata": {"path": NB_DIR}},
).execute()

imgs, saves, errs, fig_cells = 0, 0, 0, []
for i, c in enumerate(nb.cells):
    for o in c.get("outputs", []):
        if o.get("output_type") in ("display_data", "execute_result"):
            if "image/png" in o.get("data", {}):
                imgs += 1
                fig_cells.append(i)
        if o.get("output_type") == "stream":
            txt = o.get("text", "")
            txt = "".join(txt) if isinstance(txt, list) else txt
            if "Saved figure" in txt:
                saves += 1
        if o.get("output_type") == "error":
            errs += 1

nbformat.validate(nb)
nbformat.write(nb, NB_REPO)
print(f"DONE in {time.time()-t0:.0f}s: images={imgs} saved-figures={saves} errors={errs}")
print("figure cells:", fig_cells)
print("size:", os.path.getsize(NB_REPO), "bytes")
sys.exit(1 if errs else 0)
