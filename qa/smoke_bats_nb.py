"""Smoke-test the slim goldshtein_bats notebook end to end at reduced size.

Clones the slim notebook from the worktree, shrinks the data (2 tags per
treatment) and every MCMC / draw-heavy call to a few hundred draws, then
executes the whole thing in the unwrap venv. Purpose: shake errors out of
the notebook -> scripts wiring (all 27 figure functions must run and save)
before the hours-long full re-bake.
"""
import json
import os
import sys
import time

SLIM = "/tmp/opencode/wt-unified/notebooks/goldshtein_bats.ipynb"
OUT = "/tmp/opencode/smoke_bats/smoke.ipynb"

os.makedirs(os.path.dirname(OUT), exist_ok=True)

nb = json.load(open(SLIM))
cells = nb["cells"]


def src(i):
    return "".join(cells[i]["source"])


def patch(i, old, new, count=1):
    s = src(i)
    assert s.count(old) == count, f"c{i:02d}: expected {count}x {old!r}, found {s.count(old)}"
    cells[i]["source"] = s.replace(old, new).splitlines(keepends=True)


# ---- 2 tags per treatment (all four treatments must survive) ----------
tag_filter = '''
# [smoke] only 2 tags per treatment
_smoke = (
    meta_df.groupby("treatment", observed=False)["tag"]
           .apply(lambda s: s.dropna().unique()[:2])
           .explode()
)
_smoke_tags = list(_smoke)
raw_track_df = raw_track_df[raw_track_df["tag"].isin(_smoke_tags)]
meta_df = meta_df[meta_df["tag"].isin(_smoke_tags)]
print("[smoke] tags:", len(_smoke_tags), "rows:", len(raw_track_df))
'''
s = src(4)
assert s.strip() == "raw_track_df = load_bats_tracks()\nraw_track_df"
cells[4]["source"] = (s + tag_filter).splitlines(keepends=True)

patch(40, "tune=2000, draws=2000,", "tune=200, draws=200,")
patch(52, "tune=2000,", "tune=200,")
patch(52, "draws=2000,", "draws=200,")
patch(69, "max_draws=2000,", "max_draws=200,")
patch(82, "max_draws=2000,", "max_draws=200,", count=2)
patch(87, "N_BOOT_OBS = 1000", "N_BOOT_OBS = 50")
patch(87, "N_PPC_DRAWS = 4000", "N_PPC_DRAWS = 200")
patch(90, "n_cols=8)", "n_cols=4)")
patch(93, "max_draws=2000,", "max_draws=200,")
patch(98, "max_draws=2000,", "max_draws=200,")

json.dump(nb, open(OUT, "w"), indent=1, ensure_ascii=False)

# ---- execute in the venv -------------------------------------------------
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ["PYTHONPATH"] = "/tmp/opencode/wt-unified"

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager

PY = "/home/jake/.venvs/unwrap-env/bin/python"


class VenvKM(KernelManager):
    @property
    def argv(self):
        return [PY, "-m", "ipykernel_launcher", "-f", self.connection_file]


nbn = nbformat.read(OUT, as_version=4)
nbn.metadata["kernelspec"] = {
    "display_name": "Python 3 (ipykernel)", "language": "python", "name": "python3",
}
t0 = time.time()
try:
    NotebookClient(
        nbn,
        timeout=16200,
        kernel_manager_class=VenvKM,
        resources={"metadata": {"path": os.path.dirname(OUT)}},
    ).execute()
except Exception as e:
    print("EXECUTE RAISED:", type(e).__name__)

errs = []
saved = 0
for i, c in enumerate(nbn.cells):
    for o in c.get("outputs", []):
        if o.get("output_type") == "error":
            tb = "".join(o.get("traceback", []))
            print(f"ERROR c{i:02d}: {o.get('ename')}: {o.get('evalue')}")
            print(tb[-2500:])
            errs.append(i)
        if o.get("output_type") == "stream":
            t = o.get("text", "")
            t = "".join(t) if isinstance(t, list) else t
            if "Saved figure to:" in t:
                saved += t.count("Saved figure to:")

print(f"SMOKE DONE in {time.time()-t0:.0f}s: errors={errs or 'none'} saved-lines={saved}")
