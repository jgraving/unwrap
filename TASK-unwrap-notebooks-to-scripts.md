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


## STATE - updated 2026-11-01 (#55 UPGRADE + #56 CIRCULAR + #57 PERIODIC + #58 VMP + #59 ZEBRA E2E STUDY PRS MERGED = main 6beabf9; NEXT = locust_direction study PR) - NEXT CHAT STARTS HERE
- BOARD: main = ac2858e. Open PRs: only #43 (user's draft, unrelated).
  Merged history (oldest first): #48-#51 (annotated-branch merges),
  #52, #53 (beetles template), #54 (template polish: ppc inside
  plots.plot_ppc, uniform preprocess/posthoc/plots aliases,
  circular_modulo, within/between kappa markdown, no attribution
  lines). #55 MERGED (stack upgrade: PyMC 6 / ArviZ 1 / bambi 0.21,
  non-centered kappa priors -> 0 divergences, seeded fits,
  circular_summary, py3.12+; see UPGRADE PROGRESS below). #56 MERGED
  (circular_regression study PR = the CI-executed notebook, see
  CIRCULAR STUDY PR note below). #57 MERGED (periodic_not_circular study
  PR, 3 commits: script / notebook re-bake / README row,
  see PERIODIC note below). #58 MERGED (von_mises_policy study
  PR, 3 commits, see VON MISES POLICY note below). #59 MERGED
  (zebra E2E study PR, see ZEBRA note below). Work repo is on
  main, clean.  on main, clean.- USER DIRECTIVES 2026-10-31 (template polish; beetles is THE template -
  "we should get this right"):
  1. plot_ppc: do NOT pass the 14 parsed dict entries - put ALL the
     post-hoc computation INSIDE the plot function; the notebook cell is
     one call + plt.show(). (First suggested ppc-dict arg, then: "why
     have it in the notebook at all".) Final: plots.plot_ppc(
     idata_treat, model_treat, df, save_path=...) computes via
     posthoc.ppc_bootstrap_stats internally; markdown above the cell
     explains what it is (bootstrap-of-data vs posterior-of-model
     comparison).
  2. Import aliases - SAME in every study notebook:
     `import <study>_preprocess as preprocess`,
     `import <study>_posthoc as posthoc`, `import <study>_plots as
     plots`. All call sites use the aliases.
  3. Terminology (get it EXACTLY right, it is the point of the beetles
     notebook): kappa = WITHIN-individual concentration (how tightly an
     individual's own trial headings cluster around its own mean
     direction; model parameter, light effect + per-individual random
     effect). kappa_mu = BETWEEN-individual concentration = "the
     concentration of the individual means across the population" (user
     words); the paper's change-in-dispersion test = does the stimulus
     change kappa_mu. Markdown c00/c09/c14/c17 must state this split
     (c00 currently muddles them; c14 maps the figure panels to it).
  4. Naming: use `circular_modulo` (the primary name in
     unwrap/link/numpy.py + pymc.py; `shifted_modulo` is only an alias)
     in the notebook + all study scripts; "we'll eventually change to
     wrap but let's switch to circular_modulo for now". Swap includes
     the plot-fn parameter name and the model cell's LINKS["..."] key
     (registry is built from globals() so both keys resolve to the same
     func). The other 6 notebooks keep `shifted_modulo` until their own
     study PR.
  5. NO attribution lines: c23 md had "A supplementary posterior
     predictive check (contributed by James) compares ..." - remove the
     "(contributed by James)" and never add such lines to PR bodies /
     commits / notebooks.
  6. DIVERGENCES (user: "the model has divergences you need to check
     for ... rhat is meaningless if there are divergences"): the model
     needs a diagnostics pass that reports the divergence count FIRST
     and only trusts Rhat when it is zero. Fixing the divergences
     requires the stack upgrade: bambi with the PER-PRIOR `noncentered`
     kwarg (bambi >=0.19) -> pymc 6 / pytensor 3 / arviz, and
     noncentered=True on the NORMAL priors (the kappa random effects,
     e.g. softplus('1|Individual_kappa'); the mu random effects are
     VonMises). USER DECIDED ORDER: land the polish PR FIRST
     ("push the polish first so you don't lose my comments"), THEN the
     upgrade PR. Upgrade PR then re-bakes the beetles notebook on the
     new stack (outputs shift) and the remaining 5 studies run on the
     new stack.
     NOTE (found during polish): the beetles c10 fit is UNSEEDED
     (`model_treat.fit(tune=2000, draws=2000, target_accept=0.9,
     idata_kwargs={"log_likelihood": True})` - no random_seed): the
     kappa_mu point estimate moved 3.39 -> 5.32 between two bakes of
     the same code. The upgrade PR should add an explicit random_seed
     (and raise it for zebra's c09 too, which is also unseeded).
- BEETLES TEMPLATE (PR #53 MERGED = main 0e6e43b; branch
  notebooks-to-scripts-beetles @ 1612f9b was the head). Gates were all
  green; merged with CI fully passing (build + execute + full matrix;
  the notebooks.yml "execute" job = circular_regression re-run, still
  green - the #36-style cwd-import anchor in setup is the pattern for
  CI-executed studies). 3 commits on top of old main (57e5b0a):
  * 937bbdc unwrap/plots.py: save_figure gains additive svg=False kwarg
    (+ test) - salvaged from closed #50 (a85d8a5); beetles saves 3 of 4
    figures as SVG.
  * 9213e97 study scripts: beetles_preprocess.py (load_beetle_model_table,
    verbatim E2E data cell); beetles_posthoc.py (kappamu_individual,
    kappamu_per_draw, ppc_bootstrap_stats - verbatim incl. upstream
    quirks: intercept kappa_mu override via inverse_bessel_ratio, per-
    draw intercept-failure-into-light-slot fallback; dead to_360/
    circ_xerr_from_interval helpers + never-displayed HDI tables
    dropped); beetles_plots.py migrated to unwrap.plots (FIG_DIR
    re-anchored <repo>/figures/beetles/, _save() via house save_figure,
    _restore_ticks via house tick_style, every public fn gains save=True
    + explicit save_path; figure bodies byte-identical).
  * 1612f9b notebook re-baked end to end (172 s, 0 errors): 25 cells,
    16 kept byte-identical (all model/fit + az.summary diagnostics);
    setup = %matplotlib inline + imports of the 3 sibling scripts via
    cwd-independent path + plots.set_style(); data/post-hoc/plot cells
    one-line + plt.show(); 4 figure cells (8,15,16,24) each 100 dpi
    image/png + "Saved figure to: figures/beetles/<name>.svg|pdf"; saved
    3x SVG + 1x PDF valid; cross-check line printed; 106 pytest passed;
    nb 2750 KB -> 688 KB.
- STASH HAZARD (bitten): beetles_plots.py is a TRACKED file (came in via
  #36), so `git stash` on the study branch stashes the migrated version;
  a failed `git stash pop` + `git stash drop` then LOSES the migration
  (rebuilt it from the unreachable WIP commit via `git fsck
  --unreachable` + `git show <sha>:<path>`). Before dropping any stash:
  `git stash show --stat` and confirm it holds nothing current. Worktree
  state of record right now: on notebooks-to-scripts-beetles, clean tree.
- PHASE 0 DONE:
  * PR #52 (main e571342) landed all 7 annotated branches: #33 circular,
    #34 periodic, #35 von_mises, #36 beetles (definitive), #37
    locust_direction, #45 locust_density (E2E), #46 bats (E2E) - each with
    its annotated markdown + existing *_plots.py (all still single-DPI;
    two-res migration happens per-study). Conflict resolutions documented
    in the merge commits (.gitignore superset x3; README E2E rows from
    #45; datasets.py = main, the #46-side goldschtein_bats typo was NOT
    the real data-repo dir name; keep main spelling).
  * Data-cell fixes in #52: beetles c6 -> load_beetle_security_light()
    (verified: 196 rows, 100 OFF/96 ON, 10 beetles); locust_direction c2
    -> load_directional() + rediscretize_directional() (verified: 2435
    steps / 52 animals = historical step table). Intro md updated. Baked
    outputs predate the rewrites but are data-identical; re-baked by each
    study PR.
  * PR #40 MERGED (main 57e5b0a): zebra = annotated + E2E unwrap-data data
    cell + scripted plot split (notebooks/zebra_plots.py now on main;
    c15/c16 = statistics inline + plots.plot_summary/plot_ppc calls) +
    predict(kind="response_params") kappa fix via pm.compute_deterministics
    (bambinos/bambi#1084). Merge-conflict resolution: .gitignore=main;
    README=main rows + combined zebra row ("fully rerun and validated");
    zebra notebook = #40's version (main's was still pre-split inline).
- GATE LEARNED: main is PUSH-PROTECTED (GH013 "changes via pull request
  only"): direct pushes to main are REJECTED. Pattern: push a branch,
  `gh pr create --base main --head <branch>`, then `gh pr merge N --merge`
  (fresh PRs merge CLEAN - the stale-mergeable race only bit the long-
  closed PRs). docs/* branches are NOT protected (push directly, as the
  handoff mirror does).
- JQ GOTCHA: --jq expressions with backslashes (.\(...)) get mangled in
  this shell ("unexpected token"). Use single-quoted jq without
  backslashes: --jq '.mergeable + " / " + .mergeStateStatus'.
- ZEBRA STATE (post-#40): annotated md; c04 = E2E (load_zebra_bearings
  blocks from unwrap-data; raw = load_zebra_bearings_raw, 3,897,536
  frames, available for provenance reproductions); c09 fit UNSEEDED
  (tune=2000, target_accept=0.8) - STILL RAISE WITH USER on the E2E redo
  (baked outputs shift every run; propose adding an explicit random_seed).
  zebra_plots.py = single-DPI, its own set_style(), no save/save_path.
  Remaining zebra work = the E2E study PR: two-res migration +
  save/save_path kwargs + raise the seeding question + full re-bake.
- BEETLES STATE (post-#52): annotated 25-cell, already thin+scripted
  (plot cells are one-liners); E2E data cell; beetles_plots.py is
  single-DPI with its own set_style() (FIG_DIR under notebooks/figures).
  Template PR = main remaining job for the first study (see work order 1).
- All 8 *_plots.py on main: single DPI=300, own set_style(),
  FIG_DIR=Path(__file__).parent/"figures"/<study> (beetles/zebra) or
  parent.parent/"figures"/<study> (others) - ALL pre-date the two-res
  unwrap.plots convention; every study PR migrates its script to
  unwrap.plots set_style/save_figure + `save` (default on) / `save_path`
  kwargs, with saved files .gitignored (figures/**, notebooks/figures/
  already ignored).
- Upstream bugs kept verbatim (do NOT silently fix): beetles old c19
  az.summary reads idata_treat (idata_nolight_km seems intended); old
  c38 per-draw loop's first except-branch writes the intercept failure
  into the light slot.
- Per-prior noncentered (bambi >=0.19) deferred to the PyMC 6 migration;
  model-level noncentered=False stays where needed (beetles, zebra).
- UPGRADE PR (was the "terminal task"; PULLED FORWARD per user - it is
  the second PR of the beetles template pair, landing right after the
  polish PR, and BEFORE the other 5 study PRs):
  * WHY: the beetles model has DIVERGENCES (user); Rhat is meaningless
    while they exist. Fix = per-prior `noncentered` kwarg (bambi
    >=0.19) with noncentered=True on the NORMAL priors (the kappa
    random-effect group priors, e.g. softplus('1|Individual_kappa');
    the mu random effects are VonMises and stay centered) + new stack.
    The notebook gains a diagnostics pass that reports the divergence
    count FIRST (sample_stats["diverging"]) and only trusts Rhat at
    zero divergences.
  * PINS TO BUMP (requirements.txt, read 2026-10-31): currently
    bambi>=0.15.0, pymc>=5.18.0,<6, arviz>=0.20.0,<1,
    pytensor>=2.25.5, numpy>=1.26, pandas>=2.1, scipy>=1.11.2,
    matplotlib>=3.8,<3.11 (the mpl cap exists only because arviz 0.x
    uses matplotlib.style.core, removed in mpl 3.11). Target:
    bambi>=0.19 (or the version with per-prior noncentered - VERIFY
    which release introduced it), pymc>=6,<7, arviz 1.x, pytensor>=3,
    drop the matplotlib cap if arviz 1.x no longer needs style.core.
  * VENV BEFORE UPDATING: back it up (pip freeze >
    /tmp/opencode/venv_frozen_pre_upgrade.txt, and a full copy for
    rollback). Current (2026-10-31): bambi 0.17.2, pymc 5.28.5,
    arviz 0.23.4 (split arviz-base/plots/stats 0.8.x), pytensor 2.38.3,
    numpy 2.4.6, pandas 3.1.0, scipy 1.18.1. A prior upgrade probe was
    ROLLED BACK (pin clash) - resolve it deliberately this time.
  * IN-PACKAGE RISK: `from pymc.distributions.transforms import
    LogExpM1` in unwrap/prior.py - path moves under the new backend;
    check and move. Verify the `pm.compute_deterministics` signature in
    pymc 6 (used by beetles_posthoc.ppc_bootstrap_stats and the zebra
    #1084 workaround). Check whether bambinos/bambi#1084 (predict
    kind="response_params" kappa bug) is fixed upstream; if yes,
    workarounds drop in the LATER study PRs only (keep verbatim in
    already-merged PRs).
  * CI: .github/workflows/notebooks.yml executes the
    circular_regression notebook - it must stay green on the new stack
    IN THIS PR (touch minimally if needed; its full rework comes in
    its own study PR).
  * Re-bake the beetles notebook after (outputs shift: sampling
    defaults; expect the divergence count to reach zero with
    noncentered normal priors).
- UPGRADE PROGRESS (in flight on branch upgrade-bambi-pymc6;
  venv backup at /home/jake/.venvs/unwrap-env-pymc5-backup + freeze at
  /tmp/opencode/venv_frozen_pre_upgrade.txt):
  * INSTALLED NEW STACK: bambi 0.21.0, pymc 6.3.2, pytensor 3.3.6,
    arviz 1.3.0 (split arviz-base/plots/stats 1.3.x), py3.12. All four
    require PYTHON >= 3.12 (pytensor <3.15) -> tests.yml matrix cut to
    3.12-3.14 with the minimum job on 3.12 (done); setup.py
    python_requires >=3.12 (done); requirements.txt bumped (done:
    bambi>=0.21.0, pymc>=6.3.2,<7, arviz>=1.3.0,<2, pytensor>=3.3.6,<4,
    mpl cap dropped - the <3.11 cap was only for arviz 0.x).
  * PACKAGE FIXES (done, pytest 106 green): unwrap/data.py
    sample_prior_predictive(samples=1) -> (draws=1) [pymc 6 rename];
    unwrap/circular.py: az.hdi(hdi_prob=) -> (prob=),
    az.circular_hdi(...) -> az.hdi(..., circular=True [returns (lo,hi),
    (-pi,pi] for 1-D]), NEW circular_summary(idata, var_names,
    circ_var_names, *, hdi_prob=0.94, skipna, group) reproducing arviz
    0.23's circ_var_names rows (first 5 cols on the circle:
    scs.circmean/circstd high=pi low=-pi, circular hdi, batched
    circular MCSE copied from 0.23 _mc_error; columns = the 0.23 layout
    mean/sd/<p>% lower|upper/mcse/ess_bulk/ess_tail/r_hat by
    IdxLabeller row labels via az.xarray_var_iter).
    tests/test_model.py uses circular_summary.
  * NOTEBOOK CELLS SWAPPED (circular_regression c11 + its c02 import;
    beetles c12/c13 + c02 import): az.summary(..., circ_var_names=...)
    -> circular_summary(...). The OTHER four notebooks' circ cells
    (zebra, locust x2, bats) are left for their own study PRs - they
    will crash on the new stack until then (documented).
  * ARVIZ 1.3 API MAP (for the remaining ports): summary: hdi_prob ->
    ci_prob, kind='stats' = mean/sd/ci ONLY (mcse lives in kind='all'
    as mcse_mean/mcse_sd), default CI is the ETI (ci_kind='hdi'
    required for HDI), row labels "mu[0]" style (IdxLabeller).
    rhat/ess -> DataTree (values via r[name].values; az.ess
    method='bulk' - ess_method kwarg gone). plot_posterior REMOVED
    (locust_density uses it). plot_trace signature (dt, ...).
    InferenceData groups are DataTree: .values/.sizes/[var] still work
    (proven by tests + bakes).
  * EMPIRICAL DIVERGENCE CHECK (done; seed 42, bake-size fits, OMP=1):
    centered kappa RE priors -> 1 divergence, rhat 1.0032;
    noncentered=True on the two Normal kappa RE priors ONLY (mu RE
    priors are VonMises, stay centered) -> 0 divergences, rhat 1.0020.
    Fits saved /tmp/opencode/divfit_{centered,noncentered}_seed42.nc.
  * NEW-STACK LANDMINES FOUND + FIXED IN THIS PR (all probed
    empirically, not guessed):
    - arviz 1.3 has NO InferenceData class: model.fit() and
      az.from_netcdf() return a plain xarray DataTree (groups are
      child nodes). [var], .values, .sizes and
      int(idata.sample_stats.diverging.sum()) all work on it.
      (az.from_netcdf was also the 0.23 loader name.)
    - pymc 6: BOTH pm.compute_deterministics AND
      pm.compute_log_likelihood need ALL free RVs of the graph in the
      posterior. Noncentered kappa RE: bambi stores the transformed
      centered effects + *_sigma and drops the per-group *_offset
      factors, so both calls KeyError until the offsets are
      reconstructed. FIX: beetles_posthoc.graph_free_rvs(idata,
      model) returns the fit tree with a posterior group carrying the
      reconstructed offsets - graph computes effect = sigma * offset
      (sigma (chain,draw), offset per group (chain,draw,n_ind)), so
      offset = effect / sigma with sigma expanded along the group
      axis, 0 where sigma == 0 (exact in float). c10 rebinds:
      idata_treat = pm.compute_log_likelihood(posthoc.graph_free_rvs(
      idata_treat, model_treat), model=model_treat.backend.model).
      NOTE: model.backend is None until AFTER .fit().
    - az.hdi(hdi_prob=) -> az.hdi(prob=) (beetles_plots._mean_hdi);
      az.summary(hdi_prob=) -> ci_prob= and 1.x's default CI is the
      ETI, HDI needs ci_kind="hdi" (beetles c21, circular_summary);
      az.plot_trace(compact=) GONE -> PlotCollection:
      pc = az.plot_trace(...); pc.show() (circular_regression_plots);
      az.plot_posterior REMOVED (locust_density only -> its study PR).
    - circular_summary now matches the arviz 0.23.4 frame EXACTLY
      (verified against the backup venv on the saved noncentered
      fit): columns mean, sd, hdi_<p>%, hdi_<100-p>%, mcse_mean,
      mcse_sd, ess_bulk, ess_tail, r_hat; only the first five stats
      of a circ row (mean, sd, hdi lo/hi, mcse_mean) are
      circularized - mcse_sd/ess/r_hat stay linear - with round(3)
      on the circularized columns like 0.23's auto rounding.
  * BEETLES CELLS DONE: c02 SEED=42 (mirrors circular_regression);
    c10 = noncentered=True on the two kappa RE priors + import pymc
    as pm + fit(draws/tune 2000, target_accept=0.9, random_seed=SEED,
    progressbar=False) + compute_log_likelihood via
    posthoc.graph_free_rvs (rebound) + `print("Divergences:", ...)`;
    c09/c11 markdown document the noncentering and the
    divergence-first R-hat rule; c21 ci_prob/ci_kind; c12/c13
    circular_summary.
  * UPGRADE PR #55 (open, branch upgrade-bambi-pymc6, 4 commits):
    1) package+CI (pins, py3.12, circular.py, data.py, test_model.py),
    2) study scripts (graph_free_rvs, az.hdi(prob=), PlotCollection
    plot_trace, id_data (1,n)), 3) notebooks (beetles c02/c10/c11/c12/
    c13/c21 + circular c11/c14 + plots script), 4) requirements floors
    pandas>=2.2.0 scipy>=1.12.0. Pre-PR local bakes PASSED: beetles
    76 s 0 errors (Divergences: 0; c12 0.23-layout table; 4 figure
    cells + 4 saved files 100/300 dpi; kappa_mu 5.3335 seed 42;
    cross-check line from inside plot_ppc), circular 11 s 0 errors
    (Divergences: 0; 3 figure cells 100 dpi). CI round 1: notebooks
    execute PASS + all latest jobs pass; ubuntu 3.12 MINIMUM failed
    building pandas 2.1.0 from source (no cp312 wheel) -> fixed by
    the floors (verified cp312 manylinux wheels exist for pandas
    2.2.0 / scipy 1.12.0 / numpy 1.26.0 / matplotlib 3.8.0).
    EQUIVALENCE CHECK (done): circular_summary vs arviz 0.23.4 on the
    saved noncentered fit - circularized columns identical after
    rounding (mean/sd/hdi_3%/hdi_97%/mcse_mean/mcse_sd); only
    ess_tail/r_hat differ (arviz 1.3's own estimators: 1.000 vs 1.002,
    max delta 285 in ess_tail) - noted in the PR body.
  * CIRCULAR STUDY PR #56 (MERGED = main d8fbfde; 2 commits:
    circular_regression_plots.py / notebook). circular_regression
    plots join the house convention: set_style() (unwrap.plots
    machinery, FIG_DIR-anchored) + save (on by default) + explicit
    save_path on all 3 figure fns, all vector SVG under
    figures/circular_regression/ (gitignored). Stack landmine found in
    THIS PR: PlotCollection.savefig with a repo-relative save_path
    resolves against the KERNEL cwd (notebooks/), not the repo root ->
    trace save is anchored on FIG_DIR (FIG_DIR / Path(save_path).name,
    mkdir parents, repo-relative print via _plots.repo_rel);
    savefig.dpi=300 flows from the house rc. Notebook 20 -> 16 cells
    (dup markdown pairs + dup H1 from the merged explain branches
    consolidated; Figure-style section added; uniform
    unwrap.__file__-anchored setup; Step-3 markdown gets the
    divergence-first rule; stale "outputs cleared" sentence gone).
    Statistical cells verbatim (verified by diff vs the pre-study
    snapshot /tmp/zebra_no_wrong/circular_pre_studyPR.ipynb). Re-bake:
    11s 0 errors, Divergences: 0 (c09), 3 figure cells 100 dpi inline
    + 3 valid 300 dpi SVGs (trace/predictions/ppc). CI fully green
    incl. the notebooks execute job + 3.12 minimum.
  * PERIODIC STUDY PR #57 (MERGED = main ac2858e; 3 commits:
    periodic_not_circular_plots.py / notebook / README row).
    House-convention plots: set_style() via unwrap.plots (was ad-hoc
    rcParams), save (on) + save_path, vector PDF under
    figures/periodic_not_circular/ (was a 300 dpi PDF next to the
    script + plt.show() in the script), ppc draws + median/HDI moved
    INSIDE plot_beta_fit_ppc (beetles' plot_ppc pattern; cell = one
    call + plt.show()). KEY FINDING: this notebook had never run end
    to end on a modern stack (README said 'full rerun pending') -
    arviz 0.23 read a 2-D hdi array as (draw, shape) WITH A
    FutureWarning that it would change; 1.3 reads it as (n, samples)
    and dropped the d= kwarg -> az.hdi(y.T, prob=...) (time rows,
    draws last) gives the across-draws-per-time-point band the figure
    always claimed; also the doc'd 50% inner HDI band was never
    drawn - now is (flagged in the PR body). Notebook 6 -> 11 cells
    (uniform unwrap.__file__ setup; unused pandas import dropped;
    SEED; version dict; Figure-style section; giant cell -> simulate/
    model(verbatim)/fit; fit seeded + progressbar=False +
    'Divergences: 0' printed first; empty trailing cell dropped).
    Re-bake: 18s 0 errors, 100 dpi inline + valid 300 dpi PDF.
    Snapshot before edit: /tmp/zebra_no_wrong/periodic_pre_studyPR.*.
    PR body kept at /tmp/opencode/pr_periodic_body.md; gate
    /tmp/opencode/verify_periodic.py. CI green (tests 3.12-3.14 +
    min, build, notebooks execute).
  * VON MISES POLICY STUDY PR #58 (MERGED = main 24c2aea; 3 commits:
    script / notebook / README row). Lightest study: pure MC POLICY
    DEMO - no MCMC fit, no external data, no diagnostics - so the MC
    policy simulation IS the analysis and stays inline VERBATIM (both
    blocks verified vs snapshot /tmp/zebra_no_wrong/vmp_pre_studyPR.*);
    only the figure is in the script. No SEED cell (sampling already
    seeded per column default_rng(12345+i)); no divergence print (no
    fit). Plots: set_style() via house machinery + study extra rc
    (STIX mathtext, font 7, hairline, svg.fonttype path); save (on) +
    save_path; single 300 dpi SVG under figures/von_mises_policy/
    (was BOTH .pdf and .svg into notebooks/figures/ + plt.show() in
    the script); figure no longer hardcodes dpi=300 (would have
    forced 300 dpi INLINE under the house tiers). Notebook
    3 -> 11 cells (uniform unwrap.__file__ setup, cwd-walk gone;
    version dict; Figure-style section; policy/config + simulation +
    figure-call split). Re-bake 6s 0 errors, 100 dpi inline + valid
    SVG; visual check OK (top-row lobes sharpen toward goal with
    kappa; bottom row random walk -> straight). get_cmap (deprecated)
    still exists in mpl 3.10.9 - not a trap for now. PR body
    /tmp/opencode/pr_vmp_body.md; gate /tmp/opencode/verify_vmp.py.
  * ZEBRA E2E STUDY PR #59 (MERGED = main 6beabf9; 4 commits:
    preprocess+posthoc / plots / notebook / README). THE FIRST E2E
    STUDY: data via zebra_preprocess.load_zebra_model_table() ->
    unwrap.datasets.load_zebra_bearings() -> ensure_clone()
    (~/.cache/unwrap-data), first use shallow-clones, no datasets/ dir
    in the repo (the old relative parquet path could NOT run from this
    checkout). ~700 lines of inline post-hoc (old c13/c15/c16) ->
    zebra_posthoc: posterior_directions / summary_stats /
    ppc_bootstrap_stats (compute_deterministics bambi#1084 workaround
    KEPT verbatim, prints "response params from PyMC graph OK (kappa
    all positive)"). plots: house set_style/_save/save_path,
    FIG_DIR=figures/zebra (repo), plot fns now take (idata [,model],
    df_blocks, save_path=) and compute via zebra_posthoc internally
    (one-line cells); old saves were .pdf into notebooks/figures/ +
    plt.show() in-script + hard dpi=300 on the figure. Notebook
    18 -> 17 cells; fit seeded (random_seed=SEED=42) +
    "Divergences: 0" printed first; c10 az.summary(circ_var_names=)
    -> circular_summary (0.23 layout kept, c12 ess_bulk/r_hat filter
    unchanged). KEY SILENT-BUG FIX: bambi 0.21 posterior now carries
    the RE HYPER nodes (1|Track_kappa, kappa_*|Track_sigma) next to
    the RE value nodes sharing the factor dim -> the old "robust to
    naming" |Track pattern match grabbed the kappa hyperprior
    (mu_g_name etc.) -> all mu estimates shifted ~O(1); lookups now
    exclude *_kappa/_sigma suffixes (documented in the posthoc
    module). Probe/smoke before bake: /tmp/opencode/probe_zebra.py +
    smoke_zebra.py (both green; Divergences: 0; Dataset obs_dim
    __obs__ via compute_deterministics works on the CENTERED fit -
    no graph_free_rvs needed here). The tight_layout polar UserWarning
    in the summary figure cell is PRE-EXISTING (old bake c15 stderr),
    kept verbatim. Re-bake: full 4x(2000+1000) fit on the new stack;
    snapshot before edit /tmp/zebra_no_wrong/zebra_pre_studyPR.*.
    PR body /tmp/opencode/pr_zebra_body.md; gate
    /tmp/opencode/verify_zebra.py.
  * AFTER MERGE: next = locust_direction -> locust_density -> bats
    (each E2E on the new stack, unwrap.datasets.load_directional /
    load_density_order / load_bats_* accessors). Their notebooks
    still use az.summary(circ_var_names=...) / az.plot_posterior /
    az.plot_trace(compact=) -> port with the API map above + the
    graph_free_rvs pattern if they get noncentered + the
    ci_prob/ci_kind and hdi prob= renames + the az.hdi 2-D-AXIS
    change from #57 + the ZEBRA RE-HYPER-PATTERN-MATCH trap (any
    "|X" pattern matching on the posterior: exclude *_kappa/_sigma
    nodes). E2E check upstream first: bambinos/bambi#1084 (predict
    kind="response_params" kappa angle-wrap bug) - the
    compute_deterministics workaround pattern is now in
    zebra_posthog + beetles_posthoc.
- VENV: /home/jake/.venvs/unwrap-env - PyMC6 stack (bambi 0.21.0,
  pymc 6.3.2, pytensor 3.3.6, arviz 1.3.0, py3.12); unwrap is
  EDITABLE-installed from the work repo (pip check clean). Rollback:
  full copy /home/jake/.venvs/unwrap-env-pymc5-backup + freeze
  /tmp/opencode/venv_frozen_pre_upgrade.txt.
- RE-BAKE RECIPE (proven, /tmp/opencode/rerun_nb.py <nb-path>): nbclient
  NotebookClient + KernelManager.argv override to
  [/home/jake/.venvs/unwrap-env/bin/python, -m, ipykernel_launcher, -f,
  conn_file]; resources={"metadata": {"path": <repo>/notebooks}} (kernel
  cwd); LAUNCH env: PYTHONPATH=<worktree>, MPLBACKEND=Agg,
  OMP_NUM_THREADS=1. GOTCHAS: serialized stream outputs store `text` as a
  LIST of lines (join before substring checks); FIG_DIR must be
  ABSOLUTE: Path(__file__).resolve().parent.parent / "figures" /
  "<study>" (a relative string resolves against kernel cwd=notebooks/
  and the file lands in notebooks/figures/ while the printed repo_rel
  still looks right - it lied); the setup cell MUST start with
  `%matplotlib inline` - under Agg a bare plt.show() emits NO cell
  output and the bake reports images=0 with zero errors (beetles bake
  #1 hit this: 217 s, 65 KB nb, no figures in the notebook).
- VERIFY GATES: (a) clean end-to-end run, 0 errors; (b) every figure
  cell has image/png AND a "Saved figure... figures/<study>/..."
  stream line; (c) no plotting/preproc residue in notebook source (no
  `def ` / plt.figure/plt.subplots/sns./contourf/fill_between); (d)
  pytest tests/ green; (e) ipynb <~1 MB; (f) inline images really are
  ~100 dpi (PIL Image.open(...).info["dpi"]); (g) saved figure files
  valid + nb size noted in the PR.
- Helper scripts in /tmp/opencode: build_beetles_nb2.py (builder
  pattern: KEEP list + kept(i) + SETUP magic - copy per study),
  rerun_nb.py <nb-path> (generic re-bake runner; rerun_beetles.py is
  the older fixed-path ancestor), check_bp.py (plots-script residue
  scan), smoke_bp.py (Agg import + save=False/save-default smoke),
  verify_bake3.py (bake gate checker: errors, per-figure image+save
  pairs, cross-check line, executed-summary counts, PIL inline dpi,
  saved-file validity; kept-cells byte-identity via the builder's KEEP
  list against a pre-refactor snapshot), pr_beetles2_body.md (PR #53
  body - the per-study body template), pr_phase0_body.md (Phase 0
  body), /tmp/zebra_no_wrong/beetles_old36.ipynb (beetles snapshot;
  per study: `git show origin/main:notebooks/<x>.ipynb > <snap>` first).
- Work repo /tmp/opencode/wt-unified: on main (= origin/main = 0e6e43b,
  clean tree; the beetles branch was deleted post-merge and its old
  remote ref is replaced). The OLD origin/notebooks-to-scripts-beetles
  (pre-Phase-0 base @ 0b03a22) was superseded by force-push; the
  re-baselined annotated-base work is what's on main now. User's live
  repo /home/jake/projects/unwrap UNTOUCHED (docs/notebooks-explained @
  277c7ae + 6 uncommitted WIP files; PR branches are authoritative).
- Handoff mirror: this file <-> origin/docs/notebooks-to-scripts-handoff.
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

0. PHASE 0 - DONE (2026-10-31): #52 merged (7 annotated branches +
   beetles/locust_direction data-cell E2E rewires), #40 merged (zebra
   scripted split + predict fix), #48-#51 closed with comments. main is
   push-protected - future main changes go via PR (pattern documented in
   STATE). a85d8a5 (save_figure svg kwarg + test) and #51's
   cwd-independent import pattern still to salvage into the first new PRs
   (policy PR, circular PR respectively).

1. **beetles** (TEMPLATE) - off #36 base: migrate beetles_plots.py to the
   two-res unwrap.plots (save/save_path kwargs); extract preprocessing +
   kappa-mu post-hoc into beetles_preprocess.py / scripts; thin to markdown
   + one-line calls; re-bake; verify gates; 1-2 commits; one PR.

2. circular_regression - off #33 base. CI-executed tutorial: the PR must
   keep .github/workflows/notebooks.yml green (repo-root kernel, image/png
   assert); inspect the duplicated step headings first.

3. periodic_not_circular - off #34 base.

4. von_mises_policy - off #35 base.

5. zebra_escape (E2E) - #40 NOW MERGED: main's zebra is annotated + E2E
   data cell + scripted plot split (zebra_plots.py on main) + predict
   fix. Remaining: two-res migration of zebra_plots.py + save/save_path
   kwargs, RAISE the unseeded c09 fit with the user (suggest explicit
   random_seed), full re-bake.

6. locust_direction_mixture_model (E2E) - off #37 base; raw data via
   load_directional() + preprocess.rediscretize_directional (both in
   main via #47, ported from the closed #44 work).

7. locust_density_order (E2E) - off #45 base (already annotated + E2E);
   migrate locust_density_plots.py to the two-res system + save/save_path;
   re-bake; verify gates.

8. goldshtein_bats (E2E) - off #46 base (already annotated + E2E); the
   user's ask: markdown + refactored plotting + pre/post-processing into
   scripts; largest notebook (96 cells) - biggest job, last.
