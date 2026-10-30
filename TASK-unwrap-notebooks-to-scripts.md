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

## STATE AT FREEZE (next chat starts here)
- Board: only #44 #43 #40 open. (10 superseded PRs closed.)
- Next action: build & open the INFRA PR, then notebook #1 = **beetles**.
- Working clone: **/home/jake/projects/unwrap** (persistent; this file is committed

  on branch docs/notebooks-to-scripts-handoff there). Repo state there:

  origin/feat/notebooks-to-scripts = ca04985 (latest infra+work),

  origin/main = 1e21fee. Stash holds the old tmp-clone notebook edits.
  Re-verify with `git status` + `gh pr list` first.
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
