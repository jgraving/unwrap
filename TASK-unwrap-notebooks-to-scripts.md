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
- PR strategy: **INFRA PR first**, then each notebook PR forks from main.
- Infra PR = unwrap/{plots.py, circular.py, datasets.py, preprocess.py}
  + their tests, from origin/main.
- After infra PR is up, close #44 (feat/locust-data-infra; bundles the same
  datasets.py/preprocess.py with unrelated notebook edits).
- Keep open: #43 (feat/circular-rhat), #40 (pr/fix-zebra-ppc-kappa).
- All other PRs already CLOSED (30,33,34,35,36,37,41,42,45,46). Branches kept.

## STATE — updated 2026-10-31 (infra #47 MERGED; beetles #48 open + green on main; next = periodic_not_circular) — NEXT CHAT STARTS HERE
- Board: #48 (beetles, open, ALL CI green), #43 draft, #40 open (user's). #44 closed
  (superseded by #47). #47 merged into main as merge commit 3417cb7 (infra modules +
  tests + .gitignore + repo_rel Windows as_posix fix @ db5cc15).
- #48 = notebooks-to-scripts-beetles @ 0b03a22, rebased DIRECTLY ON TOP OF main after
  #47 merged (force-pushed). 2 commits: scripts (3c2bb0c), thin re-baked notebook
  (0b03a22). Notebook 5.24 MB — user-approved as-is ("if they're just trace plots I
  don't care. Just push it we can fix later"). Trace-figure shrink (scalar-only
  var_names / explicit figsize / lower inline dpi) remains a noted follow-up.
  PR body: /tmp/opencode/pr48_beetles_body.md (now documents the per-prior analysis).
- PER-PRIOR NONCENTERED (investigated at user request): YES, it exists —
  `Prior(..., noncentered=None)` added in bambi 0.19.0 (verified absent in 0.18.0,
  present in 0.19/0.20/0.21; backend: `effective_nc = prior.noncentered if
  prior.noncentered is not None else self.noncentered`). BUT every 0.19+ release
  requires pymc>=6 + pytensor>=3, conflicting with the repo's deliberate pins
  (pymc>=5.18.0,<6; arviz>=0.20.0,<1) → #48 KEEPS the model-level
  `noncentered=False`; per-prior form is the natural first step of the PyMC 6
  migration. WARNING: a probe upgraded the venv to bambi 0.21/pymc 6.3/arviz 1.3 —
  RESTORED to the pinned stack afterwards: bambi 0.17.2, pymc 5.28.5, arviz 0.23.4,
  pytensor 2.38.3, pip check clean. Do NOT leave it on 0.21.
- TERMINAL TASK (agreed with user 2026-10-31, after the notebook PR series):
  PyMC 6 / ArviZ 1 / bambi 0.19+ migration. Checklist (audit-verified surface):
  1. Pin bumps: requirements.txt → pymc>=6,<7; arviz cap <1 dropped (1.x);
     bambi>=0.19 (0.21); pytensor>=3. CI matrix 3.10-3.14 already fine (pymc 6
     drops 3.9 only).
  2. unwrap package — one real risk: `from pymc.distributions.transforms import
     LogExpM1` in unwrap/prior.py (path moves under pytensor 3; try
     pytensor.distributions.transforms). Everything else is stable API
     (pm.VonMises/Mixture/Model in distributions.py, data.py; pt.* in link/pymc.py).
     CircularUniformPrior is a bmb.Prior(...) FACTORY, not a subclass → already
     compatible with the new noncentered kwarg. Verify via pytest.
  3. bambi >=0.19: swap model-level noncentered for per-prior
     `noncentered=False` on the VonMises group priors — affects ALL notebook model
     fits (beetles 3 cells + zebra + locust x2 + bats + policy + regression). Also
     check whether the bambi 0.17.2 predict() kappa bug (bambinos/bambi#1084) is
     fixed upstream → then decide to drop #40 / the zebra predict workarounds.
  4. ArviZ 1.0: API surface used is tiny (hdi x109, summary, loo, plot_trace,
     plot_posterior, compare, circular_hdi) — one empirical pass expected enough.
  5. Re-bake ALL 9 notebooks (PyMC 5→6 sampling defaults shift baked figures).
     Recipe reusable; ~half a day of run-and-fix total. Land the notebook PR
     series first; migration is one final focused PR.
- NEXT: notebook #3 = periodic_not_circular.ipynb (handoff says: small, expected
  1 figure cell — verify against the real file). Then von_mises_policy (1 figure
  cell), circular_regression (tiny; CI-executed tutorial — its PR must keep
  .github/workflows/notebooks.yml green: it asserts an image/png output exists),
  then zebra, locust x2, bats (heavier; bats/locust have pre-existing large
  scripts on tmp-unified to reuse/migrate — some single-DPI ones need the
  two-res system). Follow-on branches all off origin/main (now includes infra).
  Per notebook: branch notebooks-to-scripts-<study>; <study>_preprocess.py +
  <study>_plots.py on unwrap.plots (FIG_DIR=figures/<study>); thin notebook
  (setup cell with %matplotlib inline + imports + plots.set_style(); data cell;
  model fits verbatim minus imports + noncentered comment; thin one-line plot
  cells); re-bake; verify gates; 1-2 commits; one PR (body pattern:
  /tmp/opencode/pr48_beetles_body.md).
- Re-bake recipe (PROVEN, /tmp/opencode/rerun_beetles.py — copy & retarget per
  notebook): nbclient NotebookClient, KernelManager.argv override →
  [/home/jake/.venvs/unwrap-env/bin/python, -m, ipykernel_launcher, -f, conn_file],
  resources={"metadata": {"path": <repo>/notebooks}} (kernel cwd!), in the
  LAUNCHING process env: PYTHONPATH=<worktree> (THE venv's installed unwrap lacks
  the new modules — importing from repo root masks this; the kernel sees it),
  MPLBACKEND=Agg, OMP_NUM_THREADS=1. GOTCHA: serialized stream outputs store
  `text` as a LIST of lines — join before substring checks.
- Verify gates: clean end-to-end run (0 errors); every figure cell has image/png
  AND a "Saved figure to: figures/<study>/..." stream line; no plotting residue in
  notebook source (no `def ` / plt.figure / plt.subplots / sns. / contourf /
  fill_between); pytest tests/ green; ipynb <~1 MB (diagnostic-only figures: user
  discretion, documented in PR).
- Helper scripts in /tmp/opencode: build_beetles_nb.py (builder pattern),
  rerun_beetles.py, beetles_map.py (cell index map), beetles_dump.py (dump cell
  sources), verify_beetles_nb.py (kept-cells byte-identity + residue regex),
  debug_trace.py, pr48_beetles_body.md.
- Work repo /tmp/opencode/wt-unified (clean tree, on notebooks-to-scripts-beetles).
  User's live repo /home/jake/projects/unwrap left UNTOUCHED (docs/notebooks-
  explained @ 277c7ae, 6 uncommitted WIP files).
- Upstream bugs kept verbatim in beetles (not fixed): old c19 az.summary reads
  idata_treat where idata_nolight_km seems intended; c38 per-draw loop's first
  except-branch writes the intercept failure into the light slot.

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



## NOTEBOOK WORK ORDER

1. INFRA PR (unwrap/plots.py circular.py datasets.py preprocess.py + tests)

   — from main; close #44 once infra PR is up.

2. beetles.ipynb (smallest; preprocess script done)

3. periodic_not_circular.ipynb

4. von_mises_policy.ipynb

5. circular_regression.ipynb (tutorial; tiny; CI runs it — keep CI green)

6. zebra_escape.ipynb

7. locust_direction_mixture_model.ipynb

8. locust_density_order.ipynb (big; scripts exist, migrate locust_density_

   plots.py to unwrap.plots two-res, thin 6 MB -> small)

9. goldshtein_bats.ipynb (largest; scripts exist, biggest job)
