# TASK — Notebooks → Scripts (FOCUS ANCHOR — re-read before acting)

## WHAT THE USER WANTS (verbatim intent)
Focused PRs, ONE PER NOTEBOOK. For each notebook:
1. Move ALL bloat OUT of the notebook:
   - pre/post-processing  -> notebooks/<name>_preprocess.py (separate)
   - ALL plotting code    -> notebooks/<name>_plots.py  (separate)
2. Notebook keeps only: data-load cells, MODEL FITTING + post-hoc
   analysis (inference stays inline), and thin cells calling the scripts.
3. Figures must be DISPLAYED in the notebook at LOW RESOLUTION inline
   (unwrap.plots: figure.dpi=100 baked via plt.show()) AND the same cell
   saves a HIGHER-RES publication file (save_figure -> savefig.dpi=300
   vector PDF, figures/<study>/<name>.pdf).
4. Notebook small enough to load on GitHub (few hundred KB) WITH figures.

## DECIDED (with user, do not relitigate)
- Base = **annotated explain-branch notebooks** (the explanatory-markdown
  versions), NOT plain main notebooks. "newer doesn't matter as much as
  correct matters" (user, 2026-10-31); per-branch notebook links confirmed
  by the user.
- **Phase 0** = merge the annotated explain branches into main, one at a
  time, then close #48-#51 with explanatory comments (wrong base).
- One focused PR per notebook: `<study>_preprocess.py` +
  `<study>_plots.py`; the notebook = explanatory markdown + data calls +
  MODEL FIT + evaluation/LOO, inline; preprocessing transforms, ALL
  plotting, and heavy post-hoc math (extraction loops, bootstrap, KDE) go
  to scripts and are called as one-liners.
- Plot functions: `save` kwarg ON by default + `save_path` kwarg; notebook
  calls pass save_path explicitly. Two resolutions: 100 dpi inline
  (renders on GitHub) / 300 dpi saved file. Saved figure files are
  .gitignored so the repo doesn't bloat.
- **E2E** for zebra, locust_density_order, locust_direction_mixture_model,
  goldshtein_bats: raw download via `unwrap.datasets.ensure_clone()` →
  preprocess → model → evaluate → post-hoc, all within the notebook.
- CI (.github/workflows/notebooks.yml) stays limited to executing
  circular_regression; E2E notebooks re-baked locally.
- Keep open: #40 (user's zebra predict fix; zebra order = annotated PR
  first, then #40), #43 (draft, unrelated).
- Infra #47 MERGED (main = 3417cb7): unwrap/{plots,circular,datasets,
  preprocess}.py + tests + .gitignore (figures/**, notebooks/figures/,
  unwrap-data/).


## STATE - updated 2026-10-31 (post-pivot plan locked; duplicate tracks resolved by evidence; NEXT = user confirm, then Phase 0) - NEXT CHAT STARTS HERE
- BOARD: main = 3417cb7 (#47 infra merged). Open: #48 beetles @0b03a22
  (5.24 MB), #49 periodic @242db24 (53 KB), #50 policy @6d7d191 (35 KB),
  #51 circular @6798174 (694 KB) - all built on the WRONG base (un-annotated
  main notebooks); to be CLOSED after Phase 0. Salvage from them:
  (a) #50 commit a85d8a5 = additive save_figure svg=False kwarg + test;
  (b) #51 cwd-independent notebook import (CI runs the kernel from repo
  root while local bakes use notebooks/ cwd):
      sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
          unwrap.__file__)), "notebooks"))
  (c) #51's verified CI recipe (nbclient + python3 kernelspec, repo-root
  cwd, image/png assert). #40 (user's zebra fix; head b5ea56e =
  pm.compute_deterministics + own zebra_plots.py + c16 split) stays open.
  #43 draft untouched.
- EXPLAIN BRANCHES (all fetched locally as pr/explain-*; cells/md):
  #33 pr/explain-circular-regression 20/13 (NOTES: duplicated step
  headings "## Step 1 ..." + "## 1. ..." - inspect when reworking),
  #34 pr/explain-periodic 6/3, #35 pr/explain-von-mises 3/2,
  #36 pr/explain-beetles 25/13 (DEFINITIVE beetles base; already
  thin+scripted: plot cells 26 B/166 B/239 B; own noncentered=False;
  gold/grey palette; PPC figure; NO *_preprocess.py anywhere),
  #37 pr/explain-locust- 38/12 (locust_direction annotated base).
  #41 and #42 are absorbed (see resolved pairs below).
  #38 pr/explain-zebra ALREADY MERGED -> main's zebra_escape.ipynb is the
  annotated 18-cell version.
- DUPLICATE TRACKS (resolved 2026-10-31 by diff evidence; awaiting user
  nod then merge):
  * BATS: #46 superset of #42. goldshtein_bats.ipynb sources IDENTICAL at
    both heads (96 cells / 541,206 B, zero per-cell diff); #46 = #42 head
    + 2 commits (db4db33 set_style prints figure dir relative to repo
    root; 6429614 _save() emits figure-dir-relative paths). Both already
    load via unwrap.datasets (commit 465eb49 is in #42 too). USE #46.
  * LOCUST DENSITY: #45 superset of #41. Same 9 annotated md cells (3
    updated to E2E wording: "rediscretization runs in the code cells of
    this notebook"); 39 of 42 code cells verbatim; the 3 dropped cells
    were the bundled-parquet load, replaced by E2E dfXY load +
    rediscretization; bundled notebooks/datasets/*.parquet removed;
    plots script carries the path-output fixes. USE #45.
  * #44 (feat/locust-data-infra, closed): fully superseded by merged #47 -
    main's datasets.py contains ALL #44 functions (plus
    bats_lake_outline_path; 432 vs 417 L), preprocess.py byte-identical,
    .gitignore already has unwrap-data/. DO NOT merge.
- All six explain *_plots.py files: single DPI=300, own set_style(),
  FIG_DIR=Path(__file__).parent/"figures"/<study> - none use
  unwrap.plots; all need the two-res migration + save/save_path kwargs.
- ZEBRA (annotated, on main, 18 cells): c04 reads
  datasets/zebra_bearings_blocks.parquet (E2E -> switch to
  load_zebra_bearings_raw); c09 model.fit(tune=2000, target_accept=0.8)
  is UNSEEDED - reproducibility gap, RAISE WITH USER when redoing (baked
  outputs differ every run); c13 heavy post-hoc (7.4 KB) -> moves to
  scripts per the confirmed boundary; c15 2x2 summary + c16 6-panel PPC
  (bootstrap seeds 123/1000+j/2000+j; B_GROUP=2000/B_INDIV=1500) are the
  two plot cells; #40 fixes the model.predict(kind="response_params")
  kappa bug (bambinos/bambi#1084) via pm.compute_deterministics - kept
  verbatim. Zebra order (user): annotated first (done - #38 merged),
  then #40, then my E2E redo.
- Upstream bugs kept verbatim (do NOT silently fix): beetles old c19
  az.summary reads idata_treat (idata_nolight_km seems intended); old
  c38 per-draw loop's first except-branch writes the intercept failure
  into the light slot.
- Per-prior noncentered (bambi >=0.19) deferred to the PyMC 6 migration;
  model-level noncentered=False stays where needed (beetles, zebra).
- TERMINAL TASK (final PR after the notebook series): PyMC 6 / ArviZ 1 /
  bambi 0.19+ migration. Pin bumps (pymc>=6,<7; drop arviz <1 cap;
  bambi>=0.19; pytensor>=3). One real in-package risk:
  `from pymc.distributions.transforms import LogExpM1` in
  unwrap/prior.py (path moves under pytensor 3). Per-prior
  noncentered=False on VonMises group priors (affects ALL notebook
  fits). Check bambi#1084 fixed upstream (then drop #40-style
  workarounds). ArviZ 1.0 API surface used is tiny (hdi, summary, loo,
  plot_trace, plot_posterior, compare, circular_hdi). Re-bake ALL
  notebooks afterwards (sampling defaults shift). Do NOT leave the venv
  on bambi >=0.19/pymc 6 (pin clash; an upgrade probe was rolled back).
- VENV: /home/jake/.venvs/unwrap-env - bambi 0.17.2, pymc 5.28.5,
  arviz 0.23.4, pytensor 2.38.3, pip check clean. Leave it pinned there.
- RE-BAKE RECIPE (proven, /tmp/opencode/rerun_beetles.py): nbclient
  NotebookClient + KernelManager.argv override to
  [/home/jake/.venvs/unwrap-env/bin/python, -m, ipykernel_launcher, -f,
  conn_file]; resources={"metadata": {"path": <repo>/notebooks}} (kernel
  cwd); LAUNCH env: PYTHONPATH=<worktree>, MPLBACKEND=Agg,
  OMP_NUM_THREADS=1. GOTCHAS: serialized stream outputs store `text` as a
  LIST of lines (join before substring checks); FIG_DIR must be
  ABSOLUTE: Path(__file__).resolve().parent.parent / "figures" /
  "<study>" (a relative string resolves against kernel cwd=notebooks/
  and the file lands in notebooks/figures/ while the printed repo_rel
  still looks right - it lied).
- VERIFY GATES: (a) clean end-to-end run, 0 errors; (b) every figure cell
  has image/png AND a "Saved figure... figures/<study>/..." stream line;
  (c) no plotting/preproc residue in notebook source (no `def `/
  plt.figure/plt.subplots/sns./contourf/fill_between); (d) pytest
  tests/ green; (e) ipynb <~1 MB (diagnostic-only figures: user
  discretion, documented in PR).
- Helper scripts in /tmp/opencode: build_beetles_nb.py (builder pattern:
  cell keep/remap), rerun_beetles.py, beetles_map.py, beetles_dump.py,
  verify_beetles_nb.py (kept-cells byte identity + residue regex),
  pr48_beetles_body.md (PR body pattern).
- Work repo /tmp/opencode/wt-unified (clean tree; currently on
  notebooks-to-scripts-zebra, an old-base branch to be redone). User's
  live repo /home/jake/projects/unwrap UNTOUCHED (docs/notebooks-
  explained @ 277c7ae + 6 uncommitted WIP files; user confirmed the PR
  branches, not this WIP, are authoritative).
- Handoff mirror: this file <-> origin/docs/notebooks-to-scripts-handoff
  (push after each update: switch the branch, copy the file, commit,
  push, switch back).

## STATE AT FREEZE (superseded; kept for history)

- Board: only #44 #43 #40 open. (10 superseded PRs closed.)

- Next action: build & open the INFRA PR, then notebook #1 = **beetles**.

- Work repo: use /tmp/opencode/wt-unified if it exists, else `git clone

  https://github.com/jgraving/unwrap.git` and `git switch -c work

  origin/feat/notebooks-to-scripts` (= ca04985). origin/main = 1e21fee.

  ALL committed work for this task is on origin/feat/notebooks-to-scripts

  — nothing important exists only in a working tree. NOTE:

  /home/jake/projects/unwrap is the user's live repo (docs/notebooks-

  explained + uncommitted WIP) — do not commit/switch/stash in it.

ify with `git status` + `gh pr list` first.

- Python env: **/home/jake/.venvs/unwrap-env/bin/python**

  (unwrap installed; datasets load from the git data repo; matplotlib

  must be Agg; MCMC runs with OMP_NUM_THREADS=1, 2 chains, seeds per

  notebook — always seeded).

- Reference architecture files (READ first):

  - unwrap/plots.py  (127 L): DEFAULT_DISPLAY_DPI=100, PUBLICATION_DPI=300,

    set_style(figure_dir, extra_rc), save_figure(fig, figure_dir, name,

    pdf=True, png=False), repo_rel(), tick_style().

  - notebooks/circular_regression.ipynb — the thin tutorial example

    (no real data; already clear of baked output; 2 simple plots that

    still need the scripts treatment).

- Study scripts that ALREADY EXIST (on tmp-unified; reuse/extend):

  - goldshtein_bats_plots.py (230 L), goldshtein_bats_preprocess.py (1344 L)

  - locust_density_plots.py (1348 L — uses a SINGLE DPI constant, NOT the

    two-res system; needs migration to unwrap.plots),

    locust_density_preprocess.py (245 L)

  - beetles_preprocess.py (47 L, done). NO beetles_plots.py yet.

  - NO scripts for: zebra_escape, locust_direction_mixture_model,

    periodic_not_circular, von_mises_policy, circular_regression.

- Per-notebook remaining-bloat notes:

  - beetles: plot/data cells c5 c33 c36 c38 c41 (c41 is 11.9 KB!) -> move to

    new beetles_plots.py. Model fits c7 c16 c22 + summaries stay inline.

  - zebra_escape: fits are seeded; c13 = post-hoc extraction (stay inline),

    c15 = 2x2 main summary figure, c16 = 6-panel PPC bootstrap-vs-posterior

    (needs the bambi `model` for model.predict — a fresh refit is required;

    bambi Model objects do NOT joblib-pickle; refit in one script pass).

  - locust_density_order: notebook still has big inline helper blocks

    (c16/17 cell_color etc., c23/47/48/49/50/51/52/53/54 = duplicated

    normalize/softplus/circular_modulo/posterior_* helpers + a 2000+-line

    plot-heavy tail) — scripts exist but notebook is NOT yet thin.

  - bats: scripts exist; notebook still holds huge inline plot cells

    (c14,16,18,20,22,24,26,27,45,47,49,50,54,56,57,59,61,63,64,66,68,70,71,

    73,74,75,77,79,80,82,83,84,86,87,88,90).

  - periodic_not_circular, von_mises_policy: 1 figure cell each, small.

  - locust_direction_mixture_model: no scripts; check what it can share

    with locust_density_preprocess.py.

- My last 3 commits on feat/notebooks-to-scripts (b21e228, ca04985, 85c5af2)

  CLEARED baked PNG images from bats+locust+beetles notebooks to shrink

  them. That was the WRONG intermediate fix relative to the goal (user

  wants low-res figures KEPT after code moves out to scripts). The image

  removal is harmless (figures get re-baked low-res when the notebooks are

  re-run after refactoring); do not revert, do not re-add high-res bakes.

- Stale scratch in /tmp/opencode (ignore/delete as needed): zebra_ref/,

  zebra_*.py, zebra_*.joblib, bats_before_clear.ipynb, nb_backup/,

  msg*.txt — leftovers from an abandoned "zebra module" side-track that

  was stopped because the architecture is inline-analysis + scripts.



## PER-NOTEBOOK PROCEDURE (repeat, do not skip verification)

0. Re-read THIS FILE. `git fetch; git status; gh pr list`.

1. Branch from origin/main: `notebooks-to-scripts-<study>`.

2. Map code cells: data/prepost -> _preprocess.py; ALL plotting ->

   _plots.py (study script uses unwrap.plots set_style/save_figure,

   figures written to <repo>/figures/<study>/); model fit/analysis stays.

3. Thin the notebook to: imports of the two scripts, data cells, model

   cells, one-line plot cells: `plots.plot_x(...); plt.show()` where the

   function internally also calls save_figure().

   IMPORTANT: low-res inline render must actually appear in the saved

   ipynb — run the notebook to bake outputs (Agg + %matplotlib inline

   equivalent via nbconvert execute).

4. Verify ALL of: (a) notebook re-runs end-to-end without errors in the

   venv; (b) every figure cell has an image/png output (low res) and a

   stdout "Saved figure to: figures/<study>/..."; (c) notebook file

   ~<1 MB; (d) no plotting/prepost code remains inline (grep for

   `def ` / plt.figure in notebook); (e) `pytest tests/` passes.

5. Commit 1-3 focused commits on the branch, push, open ONE PR

   (title: "refactor(<study>): move plotting/preprocessing to scripts;

   low-res inline figures + publication PDFs"), link it. Then next notebook.



## NOTEBOOK WORK ORDER (post-pivot, annotated bases)

0. PHASE 0 (once user confirms the three evidence resolutions): merge the
   annotated branches into main one at a time, checking mergeable vs
   current main after each (watch .gitignore / notebooks/README.md
   cross-branch conflicts); fold the #40 merge in before the zebra redo
   (user order: annotated first, then #40). Then close #48-#51 with
   explanatory comments. Salvage a85d8a5 (save_figure svg kwarg + test)
   and #51's cwd-independent import pattern into the first new PR.

1. **beetles** (TEMPLATE) - off #36 base: migrate beetles_plots.py to the
   two-res unwrap.plots (save/save_path kwargs); extract preprocessing +
   kappa-mu post-hoc into beetles_preprocess.py / scripts; thin to markdown
   + one-line calls; re-bake; verify gates; 1-2 commits; one PR.

2. circular_regression - off #33 base. CI-executed tutorial: the PR must
   keep .github/workflows/notebooks.yml green (repo-root kernel, image/png
   assert); inspect the duplicated step headings first.

3. periodic_not_circular - off #34 base.

4. von_mises_policy - off #35 base.

5. zebra_escape (E2E) - main's annotated 18-cell + #40 fix; raw download
   via load_zebra_bearings_raw; RAISE the unseeded c09 fit with the user.

6. locust_direction_mixture_model (E2E) - off #37 base; raw data via
   load_directional() + preprocess.rediscretize_directional (both in
   main via #47, ported from the closed #44 work).

7. locust_density_order (E2E) - off #45 base (already annotated + E2E);
   migrate locust_density_plots.py to the two-res system + save/save_path;
   re-bake; verify gates.

8. goldshtein_bats (E2E) - off #46 base (already annotated + E2E); the
   user's ask: markdown + refactored plotting + pre/post-processing into
   scripts; largest notebook (96 cells) - biggest job, last.
