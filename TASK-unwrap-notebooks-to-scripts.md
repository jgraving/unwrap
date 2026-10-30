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

## STATE — updated 2026-10-30 (infra PR #47 OPEN; beetles in progress) — NEXT CHAT STARTS HERE
- Board: **#47 (NEW infra PR) OPEN** (`notebooks-to-scripts-infra` @ cc497e3, from main 1e21fee):
  the 4 modules + 2 test files + .gitignore (figures/, notebooks/figures/, unwrap-data/).
  105 tests green in venv. setup.py bats-extras deliberately NOT included (goes in the
  bats PR). #44 CLOSED (comment -> #47). #43 draft + #40 still open.
  Check `gh pr checks 47` (test jobs were still pending when this was written).
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
- REMAINING for beetles (exact plan):
  1. Rewrite notebooks/beetles.ipynb with nbformat:
     - NEW title md + Setup md (notebook currently opens "# Import unwrap" with a
       pip-uninstall/reinstall cell — drop that cell entirely; setup md says
       `pip install -e '.[notebook]'` from repo root).
     - Setup code cell: `%matplotlib inline` + imports (np, pd, plt, bmb, az,
       LINKS, inverse_softplus, softplus, CircularUniformPrior,
       `import beetles_preprocess`, `import beetles_plots as plots`) + `plots.set_style()`.
     - Data cell: `df = beetles_preprocess.load_beetle_model_table(); df.head()`.
     - Model cells c07/c16/c22 verbatim MINUS their import blocks (imports live in
       setup cell); keep every comment, the fit kwargs (seeds 123, target_accept).
     - Summary cells c09/c10, c18/c19, c24/c25 verbatim (they use az from setup).
     - Keep c14 ("#TODO delete not helpful" kappa-mu effect quantile) verbatim.
     - Keep LOO cells c28-c30 verbatim (incl. `import pymc as pm` in c28).
     - All markdown cells (c12, c15, c21, c27, c31, c34-c35, c37, c39, c40 etc.) kept.
     - Thin plot cells replacing old c05/c11/c13/c20/c26/c33/c36/c38/c41:
       plots.plot_data_scatter(df); plt.show()
       plots.plot_trace(idata_treat, "trace_treat", circ_var_names=["Intercept","Light","1|Individual","Light|Individual"]); plt.show()
       plots.plot_trace(idata_treat, "trace_treat_kappamu", var_names=["1|Individual_kappa","Light|Individual_kappa","kappa_Light"]); plt.show()
       plots.plot_trace(idata_nolight_km, "trace_nolight_km", circ_var_names=["Intercept","Light","1|Individual","Light|Individual"]); plt.show()
       plots.plot_trace(idata_nomu, "trace_nomu", circ_var_names=["Light","1|Individual","Light|Individual"]); plt.show()
       plots.plot_posterior_pdfs(idata_treat, df); plt.show()
       post_hoc = plots.kappamu_individual(idata_treat, df)
       post_hoc = plots.kappamu_per_draw(post_hoc)
       plots.plot_summary(df, idata_treat, post_hoc); plt.show()
  2. Re-run to bake (MPLBACKEND=Agg, OMP_NUM_THREADS=1, venv python). Kernel CWD MUST BE
     notebooks/ (else `import beetles_preprocess` fails): nbclient NotebookClient with
     a KernelManager subclass overriding the argv property to
     [/home/jake/.venvs/unwrap-env/bin/python, "-m", "ipykernel_launcher", "-f", conn_file]
     and resources={"metadata": {"path": "<repo>/notebooks"}}. Prior full bake took
     ~150-185 s (3 fits, 2000 draws each).
  3. Verify: 7 figure cells each have image/png output AND stdout "Saved figure to:
     figures/beetles/<name>.pdf"; ipynb < 1 MB; no `def `, plt.figure/plt.subplots,
     seaborn, scatter/fill_between left in notebook source; pytest tests/ green.
  4. Commit (amend/replace the WIP commit into 1-2 focused commits), push, open ONE PR:
     "refactor(beetles): move plotting/preprocessing to scripts; low-res inline figures
     + publication PDFs" — body: based on #47 (rebase onto main once #47 merges).
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
