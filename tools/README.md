# Editing notebooks one cell at a time

`notebook_cells.py` creates a working copy of each cell's source without printing or extracting saved outputs into the source files. It uses Python's standard library. The complete notebook is parsed in the local process, but only the requested cell needs to be read by an editor or assistant.

From the repository root:

```sh
python tools/notebook_cells.py split notebooks/zebra_escape.ipynb .notebook-work/zebra
cat .notebook-work/zebra/index.txt
cat .notebook-work/zebra/cells/0000.py
```

The index lists filenames and sizes without dumping cell content. Read and edit one source file at a time. For a long cell, read a bounded line range, for example `sed -n '1,80p' .notebook-work/zebra/cells/0000.py`. These are editing fragments, not standalone runnable scripts. Keep the original notebook as the execution entry point.

Before applying changes:

```sh
python tools/notebook_cells.py check .notebook-work/zebra
python tools/notebook_cells.py apply .notebook-work/zebra
```

Code changes are rejected by default. After reviewing an intentional code edit, pass `--allow-code` to both commands. Outputs and execution counts are cleared only for edited code cells, since their saved results may be stale. Markdown changes preserve all outputs. Cell order, IDs, metadata, attachments, and all other notebook fields are retained. No-op rebuilds preserve the exact original bytes.

The workspace contains an exact backup and its SHA-256 digest. Missing or extra source files, an altered backup, or a notebook changed since extraction cause the operation to fail. Application uses an atomic file replacement so an interrupted write cannot leave half a notebook. Cell insertion, deletion, splitting, and reordering are deliberately outside this helper's scope.

Use one worker per notebook. After applying, inspect the diff and commit a small, complete unit of work. Use a new workspace for the next edit; existing workspaces are never overwritten. The temporary workspaces are ignored by Git and are not durable checkpoints. Push completed commits to preserve progress across sessions.

For assistant workflows, never reconstruct a notebook from displayed tool output: that output may be truncated. Read the file locally, edit the extracted source, and let this helper preserve the untouched content. Keep saved images and other outputs out of conversation context. A successful preservation check does not establish scientific correctness or replace runtime validation for code changes.

Run the preservation tests with:

```sh
python -m unittest discover -s tests -p test_notebook_cells.py
```
