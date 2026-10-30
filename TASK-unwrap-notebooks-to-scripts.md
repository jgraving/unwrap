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

## STATE — updated 2026-10-30 late (beetles PR #48 OPEN; next = periodic_not_circular) — NEXT CHAT STARTS HERE
- Board: #47 (infra, GREEN, awaiting merge/decision), #48 (NEW beetles, just opened),
  #43 draft, #40 open. #44 closed (superseded by #47).
- #47 = notebooks-to-scripts-infra @ db5cc15: 4 modules + 2 test files + .gitignore
  + a repo_rel Windows fix (as_posix). #48 = notebooks-to-scripts-beetles @ e30517e,
  parented ON the infra branch (2 focused commits on top: scripts, then thin
  re-baked notebook). When #47 merges, rebase both #48 and all follow-on branches
  onto main so each PR shows only its own delta.
- Beetles DONE (user-approved state, 5.24 MB ipynb accepted for now):
  - scripts: notebooks/beetles_preprocess.py (47 L, load_beetle_model_table) +
    notebooks/beetles_plots.py (set_style, plot_data_scatter, plot_trace,
    plot_posterior_pdfs, kappamu_individual, kappamu_per_draw, plot_summary).
  - notebook: 43 cells = title/setup md, setup code (%matplotlib inline + all
    imports + plots.set_style()), data cell, 3 model fits (verbatim minus imports;
    + noncentered=False with comment — bambi>=0.17 noncentred default breaks
    VonMises group priors), summaries, TODO quantile cell, LOO, thin plot cells.
  - Re-bake recipe that WORKS: /tmp/opencode/rerun_beetles.py — venv kernel via
    KernelManager.argv override, kernel cwd = notebooks (metadata.path),
    PYTHONPATH=<worktree> in the LAUNCHING process env (kernel inherits; venv's
    installed unwrap lacks the new modules), MPLBACKEND=Agg, OMP_NUM_THREADS=1,
    ~200-240 s.
  - KNOWN FOLLOW-UP (user: fix later): the 4 arviZ trace figures are 12x24 in
    (20+ panels; vec params expand per-individual; compact=True does NOT shrink
    them in arviz 0.23.4 — a direct az.plot_trace call gives the same figure;
    posterior dims are named "Individual__factor_dim" by bambi). Options when
    revisiting: scalar-only var_names, or explicit figsize, or lower inline dpi
    for trace cells. Trace SAVE files are 300 dpi raster PNG (per user; 1-5 MB
    each, gitignored).
- NEXT: notebook #3 = periodic_not_circular.ipynb (small; 1 figure cell), then
  von_mises_policy (1 figure cell), then circular_regression (tiny; CI runs it —
  keep notebooks.yml green). For each: branch notebooks-to-scripts-<study> off
  notebooks-to-scripts-infra (or main once #47 merged); preprocess ->
  notebooks/<study>_preprocess.py, plots -> notebooks/<study>_plots.py on
  unwrap.plots (FIG_DIR=<repo>/figures/<study>); thin notebook; re-bake with the
  rerun recipe (adjust paths); verify gates; 1-2 commits; push; one PR.
  NOTE circular_regression is the CI-executed tutorial notebook — its PR must
  keep the CI execute job passing (it asserts an image/png output exists).
- Helper scripts in /tmp/opencode: beetles_map.py, beetles_dump.py,
  verify_beetles_nb.py (byte-check kept cells), build_beetles_nb.py,
  rerun_beetles.py, debug_trace.py (fits treat model, probes arviz).
- User's live repo /home/jake/projects/unwrap left UNTOUCHED.
- Work repo /tmp/opencode/wt-unified (clean tree), branches:
  - notebooks-to-scripts-infra @ cc497e3 (pushed = PR #47)
  - **notebooks-to-scripts-beetles @ 09bedbc (LOCAL ONLY, not pushed)**, parented on the
    infra commit so it is exactly main + infra + WIP. Contains:
    1. notebooks/beetles_preprocess.py (from ca04985 — done, 47 L)
    2. notebooks/beetles.ipynb — the OLD pre-refactor notebook (ca04985) still in tree
    3. notebooks/beetles_plots.py — NEW, written this session (imports OK; not yet
       exercised with real idata). Functions: set_style(), plot_data_scatter(df),
       plot_trace(idata, name, var_names=None, circ_var_names=None),
       plot_posterior_pdfs(idata_treat, df), kappamu_individual(idata_treat, df) -> dict,
       kappamu_per_draw(post_hoc) -> dict(+3 posterior arrays),
       plot_summary(df, idata_treat, post_hoc). FIG_DIR=<repo>/figures/beetles;
       uses unwrap.plots set_style/save_figure (100 dpi inline / 300 dpi vector PDF).
- **Trace-figure size (DECIDED with user, 2026-10-30):** the 4 arviZ trace
  figures are what blows the size. Their figures are ~12in x 20-40in (arviz
  expands vector params into per-index panels: 4 var_names -> 20+ panels;
  rcParams['plot.max_subplots']=20 caps them). They were ALWAYS 100 dpi —
  the size was inches, not dpi (my earlier "200 dpi" guess was wrong;
  az.plot_trace in arviz 0.23 returns an ndarray of axes, not the fig).
  Recipe: plot_trace uses compact=True (1 row/panel; 12x8in for 20 panels;
  inline ~370KB vs ~1.1MB) and saves the saved file as 300dpi RASTER PNG
  (user said "rasterize them to save space?") instead of vector PDF
  (vector was 2-2.8 MB). The 3 non-trace figures keep vector PDF.
  Expected notebook size after bake: ~2 MB (data_scatter 416KB + summary
  113KB + pdfs 84KB + 4 compact traces ~1.1MB). If user wants the strict
  <1MB gate: options = lower inline dpi for traces (72), or fewer panels.
- REMAINING for beetles (exact plan):
  1. ~~Rewrite notebooks/beetles.ipynb~~ DONE via /tmp/opencode/build_beetles_nb.py
     (input = /tmp/opencode/beetles_old.ipynb; 43 cells: title/setup md+code,
     data cell, thin plot cells, 3 model fits kept verbatim EXCEPT imports
     (moved to setup cell) + `noncentered=False` added to each bmb.Model call
     (bambi 0.17 defaults noncentered=True which only supports Normal group
     priors; VonMises group priors -> NotImplementedError; same fix as
     pr/explain-beetles track). Kept cells byte-identical (+1 trailing \n from
     nbformat normalisation; verified via /tmp/opencode/verify_beetles_nb.py).
  2. Re-run to bake: /tmp/opencode/rerun_beetles.py — CRITICAL: kernel env needs
     PYTHONPATH=/tmp/opencode/wt-unified (venv's installed unwrap lacks
     datasets/plots; running from repo root shadows locally, but kernel cwd is
     notebooks/) AND kernel cwd = notebooks (metadata.path) for the
     beetles_preprocess/beetles_plots imports. Venv python kernel via
     KernelManager.argv override (in the script). MPLBACKEND=Agg,
     OMP_NUM_THREADS=1. Prior full bake ~150-185 s.
  3. Verify: 7 figure cells each have image/png output AND stdout "Saved figure to:
     figures/beetles/<name>.pdf"; ipynb < 1 MB; no `def `, plt.figure/plt.subplots,
     seaborn, scatter/fill_between left in notebook source; pytest tests/ green.
  4. Split the WIP commit into 2 commits on top of infra (reset --soft to
     notebooks-to-scripts-infra; commit1 = the two scripts, commit2 = the
     re-baked notebook), push, open ONE PR:
     "refactor(beetles): move plotting/preprocessing to scripts; low-res inline figures
     + publication PDFs" — body draft: /tmp/opencode/pr48_beetles_body.md —
     note-based on #47 (rebase onto main once #47 merges).
- Flag for user (NOT fixed; kept byte-faithful): (a) old c19 calls az.summary on
  idata_treat where idata_nolight_km seems intended (upstream copy-paste bug);
  (b) old c38's first except-branch writes 100 into the LIGHT slot when the INTERCEPT
  fit fails (upstream bug, practically unreachable); (c) dead code dropped in the move
  (rho_circular, df_kappamu_delta frame, unused imports; "alternative estimator"
  alternatives kept as comments; figure title strings kept verbatim incl. unicode κ).
- Helper scripts in /tmp/opencode: beetles_map.py (cell index map), beetles_dump.py
  (dump cell sources by index).
- User's live repo /home/jake/projects/unwrap left UNTOUCHED (docs/notebooks-explained
  @ 277c7ae, 6 uncommitted WIP files).
- After restart: first check `gh pr checks 47`; if CI red, fix #47 before continuing.

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
