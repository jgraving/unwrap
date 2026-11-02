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


## STATE - updated 2026-11-02 (DEDUP PHASE COMPLETE: the shared helper families that were duplicated across the study scripts now live ONCE in the flat `unwrap` package, and every study's scripts are rewired onto the package twins - foundation PRs #69-#71 + the 8 study rewire PRs #72-#79 all merged, main = e02c856; every rewire proved behaviour-neutral by full E2E re-bake + house gate + byte-compare of the baked notebook vs PR-base main; all rewire branches + worktrees deleted both sides; qa/ mirror + this doc refreshed below) - NEXT CHAT STARTS HERE

## DEDUP PHASE (COMPLETE 2026-11-02; main = e02c856)

Goal: the study scripts each carried private copies of the same
circular-statistics / posterior bookkeeping / figure-export helpers.
Move each EXACT family into the flat `unwrap` package (repo root, no
sub-packages) and rewire the study scripts onto it, one focused PR per
study. House rule: `CI` smoke-imports `unwrap.prior` / `unwrap.data`
names - those files never move; notebook-facing modules have no api.md
entries so no docs changes beyond this handoff.

House flow per rewire PR (all 8 followed it):
1. worktree from `origin/main` (main moves with every merge)
2. rewire edit; `py_compile`; pyflakes vs the pre-rewire files with
   `:line:col:` stripped (expected delta = only warnings whose import
   was dropped)
3. runtime identity check where numerics swap (bit-identical or
   ulp-class max-diff vs the pre-rewire module, logged in the PR body)
4. commit 1 (code rewire) -> clear the notebook's code outputs
   (`json.dumps(nb, indent=1) + "\n"`) -> FULL E2E re-bake
   (`rerun_nb.py`, background, no timeout; NEVER two bakes at once -
   RAM) -> the study's `verify_*` gate -> BYTE-COMPARE of the baked
   `.ipynb` vs the PR-base main copy (every remaining diff must be on
   the accepted-transients list below, and be explained)
5. commit 2 (`<study>: re-bake the notebook ...`), push, PR, CI checks,
   merge, delete remote branch + worktree + local branch, `git fetch`

Accepted re-bake transients (documented per PR; nothing else is
allowed): NUTS wall-time lines; warning paths pointing at this
worktree (incl. rewire-shifted line numbers); pandas Styler `0x...`
addresses; saved pandas table/widget random IDs (`T_xxxxx`);
InferenceData repr `created_at`/`sampling_time`; the wall-clock sampler
stat arrays (`perf_counter_start`, `perf_counter_diff`,
`process_time_diff` - differ by construction); the InferenceData var
enumeration order (which also moves the plain repr's 12-of-18 sample_stats
fold window); random widget element UUIDs; arviz `PlotCollection at
0x...` repr addresses. AND (periodic only): the PPC figure's inline
PNG churn - `plot_beta_fit_ppc` samples its posterior predictive draws
with an UNSEEDED `pm.sample_posterior_predictive`, so the median/HDI
bands shift a few pixels on EVERY bake; proven intrinsic by re-baking
the untouched scripts on the same base (3893 differing px vs 3880 on
the rewire, both of 94070, same bbox).

Merged, oldest first (merge commit on main in parentheses):
- #69 locust arrow family (9f03ab8)
- #70 `unwrap.inference` (graph_free_rvs family) (4856395)
- #71 core additions (b7af719; full python suite = 142 tests)
- #72 bats (f0e1f44)
- #73 zebra (3c444fb)
- #74 locust-direction (56e3336)
- #75 locust-density (6ca54fe; posthoc +10/-37, plots +1/-6). Kept
  deliberately local: `_sorted_unique` (the package `sorted_groups`
  twin coerces ids to str, but the call sites filter the float-typed
  `order`/`density` columns by equality - str ids would blank every
  panel), `_circ_dist` (signed), `circular_hdi`/`_linear_hdi_mean`/
  `_subsample_draws` (distinct semantics). Bake 1156 s, 14/14 figs;
  compare: 43 cells, 0 source/PNG diffs; c27 InferenceData repr diffs
  fully accounted for (timestamps + wall-clock stats + var order +
  widget UUIDs; the full HTML var-item multiset matches).
- #76 beetles (85941ce; net -167 loc). Dropped the 27-line in-file
  `rho_circular` duplicate (-> `resultant_length_angles`, bit-identical
  50/50), the ive-ratio `_bessel_ratio_mrl` (3.3e-16 twin), the scipy
  mean-of-pdf `circular_kde` (-> `von_mises_kde(data, points, kappa)`,
  5.8e-16). Kept local: `circ_hdi_centered` (center-then-linear HDI;
  the package `circ_hdi` is arviz-based), `natural_sort_key` (alive at
  its one `sorted(key=...)` reference - grep for `(` misses bare-name
  references), `hierarchical_bootstrap_group` (different indivs
  signature), scipy `circmean`/`circstd` high/low idiom, `vonmises.fit`;
  plots keep `_ensure_style` (font.size 11 + hidden top/right spines on
  top of the house `DEFAULT_RC` - the figures depend on it),
  `circ_median`, `_to_360`, `_circ_hdi_centered`, `_mean_hdi`,
  `_errorbar_mean_hdi`, `_vonmises_interval_mass` (scipy cdf, no twin).
  Bake 84 s; compare: 25 cells, 0 source/PNG/widget diffs; only diff =
  NUTS wall-time 52 -> 50 s.
- #77 circular_regression (f1c95e8). One-line `_save` -> `save_by_path`
  swap (net -4 loc). `plot_trace` KEEPS its manual
  `pc.savefig(FIG_DIR / ...)`: `save_figure`'s SVG branch passes
  `bbox_inches="tight", pad_inches=0.01`, which would re-crop the saved
  trace SVG - documented in the commit message. Bake 12 s; compare:
  16 cells, 0 source/PNG/text diffs; 1 widget diff = the PlotCollection
  repr address.
- #78 periodic_not_circular (519d27f). Same one-line swap; see the
  unseeded-PPC PNG note above. Bake 18 s.
- #79 von_mises_policy (e02c856; current main). Same one-line swap.
  Bake 6 s; compare: ZERO diffs of any kind (11 cells).

Package facts the rewires relied on (verify_*.py gates + rerun_nb.py
mirrored into qa/ above, all worktree-aware: run from the worktree cwd
or pass the root as argv[1]): `save_by_path(fig, fig_dir, save_path,
*, transparent=False, dpi=300)` routes by suffix into
`save_figure(fig, fig_dir, p.stem, ...)`; `save_figure` PDF branch =
`savefig(dpi, bbox_inches="tight", pad_inches=0.01)`, SVG branch =
`savefig(bbox_inches="tight", pad_inches=0.01)` (no dpi); `A_of_kappa`
= vectorized scaled bessel ratio; `wrap_pi`/`circ_mod` =
`(x+pi)%2pi-pi`; `von_mises_kde(data, points, kappa)` wraps data via
circ_mod, `c = 1/(2*pi*i0(kappa))`, `.mean(axis=1)`;
`vonmises_logpdf(theta, mu, kappa)` stable via log-i0; arviz 1.3
InferenceData DataTree; `DEFAULT_RC` = Arial/dpis/facecolor/type-42
only.

Process fixes learned (cost real time; keep them):
- `find(lines, prefix)` in the edit scripts STRIPS each line, so
  prefixes with leading whitespace NEVER match - write anchors
  unprefixed (cost one full beetles re-run).
- after a block deletion, re-normalize blank-line adjacency around the
  survivors (`two_blanks`/`blanks_before`); a banner sits straight
  above its first def (0 blanks).
- recon greps for call sites must also catch BARE-NAME references
  (`key=natural_sort_key`), not just `name(`.
- dead-name asserts must allow the name to remain in comments (beetles
  `directional_stats` in the N.B. comments).
- pyflakes comparison must strip `:line:col:` (numbers shift after
  edits) and diff on message bodies.
- `gh pr edit --body-file` silently no-ops; use
  `gh api -X PATCH repos/jgraving/unwrap/pulls/N -F body="$(cat f)"`.
- NEVER a full-directory `git checkout --` in a baking worktree (it
  once clobbered a bake); partial `git checkout -- notebooks/` before a
  clean script re-run is the accepted restore.

The OLD final state (pre-dedup) is preserved below, starting with
"STATE - updated 2026-11-02 (COMPLETE: ...".

## STATE - updated 2026-11-02 (COMPLETE: every notebook on main has the house treatment - 5 research studies via study PRs (#64 bats = e4f74ae, #65 beetles = 8a0539f, zebra/locust x2 earlier) + periodic PR #66 + policy PR #67 + README intro PR #68; main = 3df54e4; #43 closed as superseded; all work branches + worktrees deleted both sides; this doc + qa/ final sync below)
- BOARD: main = f817823. Open PRs: only #43 (user's draft, unrelated).
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
  (zebra E2E study PR, see ZEBRA note below). #60 MERGED
  (zebra kappa REs non-centered + log-likelihood cross-check,
  see ZEBRA NONCENTERED note below). #61 MERGED
  (locust_direction study PR: 4 commits - posthoc/plots scripts,
  thin-cell E2E notebook, full re-bake, README row; see
  LOCUST_DIRECTION note below). Work repo is on main, clean.- USER DIRECTIVES 2026-10-31 (template polish; beetles is THE template -
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
  * ZEBRA NONCENTERED (#60, MERGED = main 734887c; 2 commits:
    posthoc / notebook+README). User follow-up to #59: the study PR
    kept noncentered=False everywhere (vM group priors); user wants the
    noncentering the stack upgrade was meant to allow. bambi 0.21
    implements non-centered RE ONLY for Normal priors with a tensor
    sigma (VonMises RE prior would raise NotImplementedError - verified
    in the source: priors/prior.py + backend/terms.py build()) so the
    fix is PER-PRIOR (the beetles convention): noncentered=True on the
    two kappa RE Normal priors, model level stays noncentered=False.
    Key mechanics (probes /tmp/opencode/probe_zebra_nc{,2,3}.py +
    smoke_zebra_nc.py, all green): (a) NUTS free RVs gain
    kappa_*|Track_offset (z ~ N(0,1)); bambi's stored POSTERIOR keeps
    the same 12 names as the centered fit - kappa_*|Track is now a
    Deterministic(effect = sigma*offset) - offsets are NOT in it;
    (b) pm.compute_deterministics / pm.compute_log_likelihood need
    ALL free RVs -> posthoc.graph_free_rvs(idata, model) (ported from
    beetles_posthoc) reconstructs offset = effect/sigma exactly;
    (c) fit cell now: idata = pm.compute_log_likelihood(
    posthoc.graph_free_rvs(idata, model), model=model.backend.model)
    after model.fit(...) - the PyMC6 replacement for
    idata_kwargs={'log_likelihood': True} - this DOES leak the two
    offset vars into idata.posterior, so (d) c10 diagnostics cell uses
    an explicit var_names whitelist (the beetles convention, 12 model
    params) and (e) zebra_posthoc RE pattern lookups now also exclude
    *_offset by suffix (stored effect value wanted, not the raw
    factor); (f) ppc_bootstrap_stats runs compute_deterministics on
    graph_free_rvs(idata) + NEW hard cross-check of the per-obs
    response params against the stored log-likelihood (rtol 1e-10,
    prints 'response params from PyMC graph OK; log-likelihood
    cross-check passed'). Re-bake: 4x(2000+1000) 242 s, Divergences:
    0, gate PASS, 106 tests green. PR body
    /tmp/opencode/pr_zebra_nc_body.md.
  * MERGED #61 (locust_direction study PR), main = 529d364.
    4 commits: posthoc+plots scripts (583b1f9) / thin-cell
    re-template (0e927f4, nbformat-valid, kernelspec normalized) /
    re-bake (b5fe42e) / README row (680ac53). Re-bake 376 s TOTAL
    (bimodal fit 194.6 s, unimodal 70.8 s per iopub timestamps -
    the fits did run full 4x(2000+1000)); Divergences: 0 first
    for both; bimodal diag filter EMPTY, unimodal diag filter 19
    non-empty rows (low-ESS kappa_log, same as the original bake
    - do NOT "fix" it); 12 figures 100 dpi inline + 300 dpi PDFs
    under figures/locust/; nb 1.8 MB; gate
    /tmp/opencode/verify_locust.py PASS; pytest 106 passed; CI
    green; merged --merge.
    LOCUST ARVIZ 1.3 / PYMC 6 API FINDINGS (probe + bake):
    (a) az.loo: NO scale= kwarg (log-scale is now the only
    output); ELPDDataLOO attrs: elpd/se/p/pareto_k/elpd_i
    (+n_samples, n_data_points, scale(str), warning, good_k);
    __getitem__ = getattr; loo_i/elpd_loo/p_loo names are GONE.
    Per-obs SE NEVER existed in any arviz (proven against the
    0.23.4 sdist in /tmp/opencode/arviz_old/ - it only returned
    loo_i), so the original's defensive se_i=zeros fallback was
    always taken: zero-length per-animal ELPD bars are the
    faithful port.
    (b) pm.compute_log_likelihood(idata, model=<RAW pm.Model>,
    progressbar=False) works on the locust models - the
    log_likelihood group var is directions_obs (2435 obs);
    this replaces idata_kwargs={"log_likelihood": True}.
    (c) az.hdi(hdi_prob=) -> az.hdi(prob=).
    NBFORMAT GOTCHA (study-PR rebuilds): re-emitting cells by
    hand in JSON - every cell NEEDS a "metadata" key (a missing
    one -> AttributeError in nbformat's strip_transient,
    surfaced as "missing an expected key: metadata"); and
    rerun_nb.py now FORCES kernelspec python3 (the locust nb
    carried a stale custom kernelspec "locust" -> NoSuchKernel
    at bake start). Rebuild script kept at
    /tmp/opencode/locust_rebuild_nb.py (needs the original
    restored from /tmp/zebra_no_wrong/ before re-runs).
    LOCUST_DIRECTION PROBE FINDINGS (2026-11-01,
    /tmp/opencode/probe_locust_dir.py - SEGFAULTED, exit 139):
    (1) SPARSE animal_id TRAP: rediscretize_directional(
    load_directional()) gives 52 animals but animal_id is SPARSE
    (0..53, missing {1,44}) -> the notebook model indexes
    mu[individuals] with RAW ids (max 53 vs num_individuals 52)
    -> C-level out-of-bounds -> segfault 'All initialization
    points failed'. FIX for the study PR: reindex animal_id to
    dense 0..51 in the data cell (e.g.
    df['animal_id'] = pd.factorize(df['animal_id'], sort=True)[0])
    BEFORE any model (posthoc c18/c27 assume dense 0..N-1 too).
    (2) The locust_direction data cell is ALREADY E2E (c02:
    load_directional + rediscretize_directional) -> NO preprocess
    script needed; migration = posthoc (~1100 lines: c18 loo-by-
    individual, c21 population density via pm.compile, c22/23
    bootstrap histogram+KDE with numba/joblib, c27 individual
    PPC 52-panel, c30-c35 hypothesis tests) + plots house
    convention (741 lines, has DUPLICATED mlines/mpatches
    imports) + az.summary(circ_var_names=) x2 (c12/c16 ->
    circular_summary; c12 filters _unsorted, c16 excludes _raw)
    + az.hdi prob= rename (5 sites, all 1-D -> no axis
    transpose) + seed fits + Divergences-first.
    (3) API probed GREEN on this stack: pm.sample still accepts
    return_inferencedata + idata_kwargs={'log_likelihood': True};
    pm.compile([angles], logp) / pm.logp(dist, value).eval() /
    pm.Mixture.dist(w=, comp_dists=) all present in pymc 6.3.2;
    az.loo pointwise works; circular_summary on 2-D circ vars
    (mu_pop shape (2,), mu (52,2)) needs a probe check. Their notebooks
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
    NEXT = locust_density study PR (E2E via load_density_order),
    then bats (load_bats_tracks / load_bats_info /
    load_bats_lake_outline) - same recipe: probe -> posthoc +
    plots scripts -> thin-cell rebuild -> smoke -> full re-bake ->
    gate -> pytest -> push + PR + merge (RE-BAKE RECIPE + VERIFY
    GATES above; helpers in /tmp/opencode are reusable:
    rerun_nb.py (generic, now forces kernelspec python3),
    verify_*.py gate pattern, pr body template
    /tmp/opencode/pr_locust_dir_body.md).
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
- Work repo /tmp/opencode/wt-unified: on main (= origin/main = f817823,
  clean tree after the #62/#63 merges; both PR branches deleted local +
  remote; the wt-dirfix worktree for the #63 fix was removed after
  merge). A bats worktree /tmp/opencode/wt-bats (tmp-stack-42) exists
  from earlier probing. The OLD origin/notebooks-to-scripts-beetles
  (pre-Phase-0 base @ 0b03a22) was superseded by force-push; the
  re-baselined annotated-base work is what's on main now. User's live
  repo /home/jake/projects/unwrap UNTOUCHED (docs/notebooks-explained @
  277c7ae + 6 uncommitted WIP files; PR branches are authoritative).
- Handoff mirror: this file <-> origin/docs/notebooks-to-scripts-handoff.

## LOCUST DENSITY STUDY PR #62 (MERGED = part of main f817823)
- Branch locust-density-study, 3 commits: scripts (preprocess+posthoc+
  plots) / thin-cell notebook + full re-bake / README row. GATE PASS
  (verify_locust_density.py: 0 errors, 14 figure cells each 1 inline
  ~100 dpi png + "Saved figure to:" line, valid PDFs >1 kB, nunique
  440, "Density × Order" best model, N=18641, both diagnostics
  filters Empty). pytest 106. ipynb 1.97 MB (c30 summary-table HTML is
  the bulk; the merged locust_direction sibling is 1.83 MB - same
  class, accepted).
- PREPROCESS: notebooks/locust_density_preprocess.py - rediscretize at
  6 body lengths on the shared unwrap.preprocess arc-length primitives;
  carries order/density through, renumbers animal_id densely 0..439.
- POSTHOC: notebooks/locust_density_posthoc.py - 11 entry points:
  screening_table (per-animal mean resultant heading), fit_table,
  attach_log_likelihood (the LL path, below), model_sweep (5 seeded
  candidate fits -> PSIS-LOO delta-ELPD + "Best model:" print),
  alignment_stats / elpd0_stats / polar_pdf_stats / coefficient_stats /
  data_cell_bootstrap / model_cell_ppc / individual_ppc. Magic
  constants kept verbatim (grids 50/80, HDIs .68/.9/.95, seeds 123/2025,
  B 1000/600/800). Bake pins: p_loo 224.56, N/p_loo 83.01, fraction
  k>0.7 0.0021, sweep Delta-ELPD Order-only -8.86 (SE 5.85) /
  Density -75.7.
- **FINAL-FIT LOG-LIKELIHOOD PATH (non-obvious, validated):**
  `pm.compute_log_likelihood` CANNOT work on the RE sigma-prior spec -
  `transform=Softplus()` on the Normal sigma prior inserts a
  `kappa_1|animal_id_sigma_log_exp` value variable that does NOT
  exact-match the non-centered posterior name (bambi 0.21). Production
  path in `attach_log_likelihood`: fit with `omit_offsets=False` (keeps
  `kappa_1|animal_id_offset` free RV in posterior - the 2nd bambi
  quirk), then `pm.compute_deterministics(posterior, ["mu","kappa"],
  model=model.backend.model)` for response params straight from the PyMC
  graph (bambi `predict(kind="response_params")` wraps kappa through the
  mean link - broken for kappa, bambinos/bambi#1084), manual von Mises
  logpdf with stable `logI0 = log(i0e(kappa)) + kappa`, then
  `idata["log_likelihood"] = xr.Dataset({"y": loglik})` (fit result is
  a DataTree - plain assignment, no `add_groups`). Then
  `az.loo(idata, var_name="y")`. Screening/sweep models have no RE and
  no transformed sigma prior -> plain `pm.compute_log_likelihood`
  (PyMC6 replacement for `idata_kwargs={"log_likelihood": True}`).
- **az.hdi 2-D AXIS GOTCHA (broke density smoke runs 3-4):** for a RAW
  2-D numpy array, `az.hdi` treats the LAST axis as the samples axis.
  Grid-shaped predictions (draws, n_points) must be transposed before
  the call, else the HDI comes back transposed (fill_between ValueError:
  'x' size 50 vs 'y1' size 40). The 3 density sites hand `align.T` /
  `logp0.T` / `pdfs.T`; 1-D sites are fine. The direction study's 2-D
  sites already transposed.
- LINK IMPORT GOTCHA: `unwrap.link.bambi.<name>` exports are Bambi
  `Link` OBJECTS (for `bmb.Model(link=...)`), not callables. NumPy link
  functions live in `unwrap.link` (module) / `unwrap.link.numpy`.
  Posthoc uses NumPy `circular_modulo` for evaluating draws, the
  `BAMBI_LINKS` registry for building sweep models; calling a Link
  object on an array is a TypeError.
- PyMC 6.3 `pm.sample`: NO `include_transformed` kwarg (removed);
  `idata_kwargs` still exists but the final fit drops it (LL off is the
  default; transformed draws always stored).
- Bake ~1 h this run (16-core box, final fit 4 chains x (2000 tune +
  1000 draws), N=18641 / 440 animals).

## DIRECTION ELDP LAYOUT FIX #63 (MERGED = main f817823)
- One-line bug in locust_direction_plots.plot_elpd_by_individual:
  panel b's "b." label was drawn with `transform=axA.transAxes` (a copy
  of panel a's line). "b." rendered at panel A's top-left (stacked on
  "a.") yet stayed a child artist of axis b. `plt.tight_layout()`
  includes every child artist in each axis's bbox, so it saw panel b as
  spanning x~=0.06 (where "b." rendered) to the right edge; the solver
  collapsed panel a to the left ~22% and pinned panel b to the right
  edge -> ~60% white space between the panels + overlapping labels.
  Fix: `transform=axB.transAxes`; notebook re-baked (both fits seeded
  identically, all other outputs unchanged). LESSON: an a./b. label
  with the WRONG axis transform silently corrupts tight_layout -
  bisect by re-adding pieces and printing `ax.get_position()`.

## BATS STUDY PR (IN PROGRESS) - branch goldshtein-bats-study pushed (6 commits off main f817823)
- Branch origin/goldshtein-bats-study: 2ea9fed preprocess module;
  965c9a4 posthoc module; 47f1d7b plots module (27 plot_* figure fns,
  two-res unwrap.plots, save/save_path, FIG_DIR=repo/figures/bats);
  d1d85e8 slim notebook (107 cells: 38 markdown byte-identical to the
  annotated base; data + 2 model cells + evaluation stay inline; every
  figure cell = comment + one-line plots.plot_*(save_path=...) +
  plt.show(); 37 figures). Then smoke/stack fixes: f1a1187 + 5db2255
  (imports/helpers the split cells lost), 82892b2 (vM log-likelihood
  via PyMC graph - see below), ac736dc (arviz 1.x elpd_i rename).
- Slim cell map (new indices): c01 setup (imports + sys.path.insert +
  pre/posthoc/plots + plots.set_style() + SEED=123 + version dict);
  c03/c04 loaders verbatim; c06 pre.rediscretize_latlon_parallel
  (60 m, gap_factor=10, gap_time_s=120, time_unit="days") + step checks;
  c08 pre.merge_tag_metadata; c11 pre.add_dist_and_angle_to_destination
  (harmonic lake field, roost point 33.112373/35.584135, METRIC_EPSG
  32636) + df_harm2 alias; c29/c34/c38 thin posthoc compute (kappa_df /
  align_df / df_track); c40 vM Bambi model - fit ported:
  `model.fit(tune=2000, draws=2000, target_accept=0.9, random_seed=SEED,
  progressbar=False)` and the per-observation log-likelihood attached
  via `pm.compute_log_likelihood(posthoc.graph_free_rvs(idata, model),
  model=model.backend.model, progressbar=False)` — bambi 0.21's own
  Model.compute_log_likelihood applies the vonmises *mean* link to EVERY
  submodel (hardcoded in ConditionalModelTerm.predict:
  `if family.name=="vonmises": invlink = angle(exp(1j*x))`), so the
  kappa linear predictor would be wrapped into an angle (bambinos/
  bambi#1084; the zebra study PR used the identical workaround,
  zebra_posthoc.graph_free_rvs + pm.compute_log_likelihood);
  c42-c45 diagnostics verbatim (az.summary, two empty-filter checks,
  full summary); c51 sigmoid_normal_percentiles thin compute; c52 vMU
  PyMC model verbatim (defines fit_vmum_model_full; pm.sample with
  idata_kwargs={"log_likelihood": True} is valid in pymc 6.3.2 and was
  working in the base bake - do NOT port it); c54 convergence flags
  verbatim; c56/c59/c64 PSIS-LOO via posthoc (models dict: vMU first);
  c69/c82/c93/c98 population + taglevel + pgoal computes via posthoc;
  c101-102 (compact pgoal) and c104-106 (faceted tag curves) REUSE
  c98's out_dist_k/out_steps_k/df_track_with_p/resp/treat_levels (E2E
  order matters; deterministic, no recompute).
- Model cells {40, 52} are the ONLY cells allowed top-level `def`s;
  the gate (below) whitelists exactly those.
- Gate: /tmp/opencode/verify_bats.py (mirrors verify_locust_density.py,
  house pattern): no error outputs; source audit per non-model cell
  (banned: def / sns. / contourf / savefig / get_cmap / out_pdf /
  out_svg / savepath / save_svg_also / az.hdi / pm.sample; figure
  cells end with plt.show() and call only plt.show() of plt; other
  cells no plt at all); 37 figure cells each with exactly 1 inline
  PNG at 100 dpi + "Saved figure to: figures/bats/<stem>" line + valid
  PDF on disk; c01 version dict (pymc+bambi) + set_style echo
  "Figures will be saved to: figures/bats"; c06 stream "tags: 67";
  c43/c44 diagnostic filters "Empty DataFrame"; c45 summary has r_hat;
  c54 "RHAT > 1.01: none ✅" + "ESS bulk or tail < 400: none ✅";
  c56 "Best model: von Mises-Uniform"; c93/c98 print the p_goal head.
  EXPECT_SAVED has the 37 stems; the two faceted p_goal cells save one
  plain PDF with the 4 treatments (4-per-page paging => single page),
  per-page _p{i} suffixes tolerated but unreachable.
- SMOKE: /tmp/opencode/smoke_bats_nb.py (clone slim nb; after c04 keep
  2 tags per treatment; fit 200/200; all max_draws=200; N_BOOT_OBS=50;
  N_PPC_DRAWS=200; PPC_BOOT_PER_DRAW=1; individual_ppc n_cols=4; run
  in the venv via nbclient; purpose = shake the notebook->scripts
  wiring across all 27 figure fns before the hours-long full bake).
  STATUS: PASS (2026-11-02) - 37/37 figures saved, 0 errors, ~2 min.
  Full re-bake RUNNING: rerun_nb.py on the slim nb (timeout=16200 per
  cell, MPLBACKEND=Agg, OMP/OPENBLAS=1, PYTHONPATH=worktree), log
  /tmp/opencode/bats_bake.log; expected ~hours (2x 2000+2000 fits,
  PPC B=1000 obs boot + 4000 draws, individual PPC 67 tags x 8 cols,
  pgoal curves 2000 draws x 2000 grid). Probe notes: vM fit at
  2000+2000 ~ 7x the 200+200 probe time (~69 s at full data =>
  ~8-15 min expected); vMU mixture likely the long pole.
- Bugs found+fixed by smoke/stack (all pushed):
  1. plots module lost imports when the cells were split:
     Polygon/MultiPolygon (shapely), prep (shapely.prepared), math;
     plot_treatment_coefficients_mu_kappa needed `post =
     idata_treatment.posterior` (was cell-local in the base).
  2. posthoc module lost the per-cell helper families the base
     redefined in EVERY cell (one shared copy added): bessel_ratio
     (robust I1/I0 + fallback), vonmises_pdf_scalar,
     vonmises_interval_mass (+scipy-cdf guard), circ_mean, wrap_0_2pi,
     quantile_band, format_treatment_nicely. (Lessons: pyflakes the
     modules after the split; the base's self-contained cells hide
     cross-cell locals.)
  3. vM `kappa > 0` ParameterValueError at c40 — confirmed NOT a fit
     problem (probe_vm_fit2.py: FIT OK at 200/200 on BOTH the smoke's
     2-tags-per-treatment subset and the full 67-tag data): it was
     bambi's Model.compute_log_likelihood applying the vonmises mean
     link to the kappa submodel (bambinos/bambi#1084) - fixed via
     pm.compute_log_likelihood + posthoc.graph_free_rvs on the model
     graph (the zebra pattern), which also passes az.loo
     (elpd -750.8, max k 0.61 on the 8-tag probe).
  4. arviz 1.x PSIS-LOO api in _loo_pointwise_elpd: loo.loo_i ->
     loo.elpd_i (ELPDDataLOO; pareto_k unchanged; the merged locust
     module already did the same rename).
  5. The delta-ELPD call sites (total-by-treatment + by-individual)
     unpacked 2 values from _loo_pointwise_elpd; the base cells had
     local 2-value copies but the module shares the 3-tuple
     (loo_i, pareto_k, loo) helper - both call sites now drop the
     loo object.
  6. tqdm_joblib was defined twice at module level (once per section,
     mirror of the base's self-contained cells); the PPC-section copy
     LACKED the @contextmanager decorator so the later, undecorated
     generator shadowed the pgoal-section copy and c87 died with
     "'generator' object does not support the context manager
     protocol". Duplicated-helper audit: 13 shadowed families total;
     all pairs byte-identical EXCEPT tqdm_joblib (decorator) and
     hdi_1d/hdi_along_samples (type hints/docstring only) - i.e. no
     behavior changes from shadowing; the module docstring now states
     the per-section structure instead of claiming one copy.
- NEXT: finish the smoke (PPC/pgoal tail: bootstrap_observed,
  posterior_predictive_summaries, individual_ppc,
  append_pgoal_post_to_df, population/tag kernel curves + last 10
  figure cells) -> full re-bake /tmp/opencode/rerun_nb.py
  (2x 2000+2000 fits + PPC, hours; run with timeout=16200 per cell,
  PYTHONPATH=worktree) -> verify_bats gate -> commit baked nb + push
  -> gh pr create -R jgraving/unwrap --base main --head
  goldshtein-bats-study (title per the house study-PR pattern) ->
  checks green -> merge -> delete branch -> README row (turn the
  "full workflow runs from this checkout" into the house "fully
  rerun and validated from this checkout (MCMC diagnostics pass; ...)"
  wording; draft in the PR-body scratch: fixed 60 m re-discretization
  with gap segmentation; roost-vs-lake destination assignment; vMU
  von Mises-Uniform mixture fit vs the vM baseline on the
  goal-relative bearings; PSIS-LOO total / by treatment / per animal
  (best model von Mises-Uniform); population + tag-level mixture
  summaries with the directed-component density and the +/-45deg
  mass around the goal; hierarchical-bootstrap PPC by treatment plus
  per-animal bootstrap-vs-posterior checks; per-step goal-component
  responsibility with kernel-smoothed distance/progress curves; vM
  per-observation log-likelihood from the PyMC graph working around
  the bambi vonmises kappa link bug, bambinos/bambi#1084) -> sync
  this doc + docs branch -> BEETLES (see section below; remaining =
  svg->pdf save paths, re-bake, gate, README row).
- Gotchas (bats-specific):
  * Base cells are self-contained (helper redefined per cell) - a
    def-only module build loses them; pyflakes (pip install pyflakes
    in the venv) catches the undefined names cheaply.
  * The shell tool `timeout` is MILLISECONDS (7200 = 7.2 s). Long
    background jobs: omit the timeout (default none).
  * arviz 1.3: InferenceData.groups is a TUPLE, not a method.
  * Kappa xlabels use "Concentration " + rf"($\kappa$)" to dodge the
    Py3.12 `\k` invalid-escape SyntaxWarning (builder audits per-cell
    compile).
  * figures/bats/ is gitignored (repo .gitignore has top-level
    figures/); only the notebook (with baked 100 dpi inline) + the
    three scripts go into the PR.
  * Builder: /tmp/opencode/build_bats_nb.py regens the slim nb from
    the base backup /tmp/opencode/goldshtein_bats_base.ipynb.bak
    (12.1 MB, 96-cell annotated base; audit asserts the figure-cell
    map, banned tokens, markdown byte-identity). Keep builder and
    hand-edits in sync (verified idempotent after the c40 change).
  * E2E data: unwrapp venv /home/jake/.venvs/unwrap-env, data from the
    user-level ~/.cache/unwrap-data clone (no re-download).
## PERIODIC + POLICY HOUSE TREATMENT (user asked: give ALL notebooks the house treatment)
- Both were already thin (scripts merged in the earlier template wave;
  11 cells each; set_style + version-dict setup; single figure cell each).
- periodic_not_circular: zero source edits needed (already saves .pdf,
  already prints Divergences first). Re-bake 17 s, 0 errors, NUTS boot +
  Divergences: 0 (progressbar=False suppresses the "Sampling..." line -
  gate checks the NUTS[sampler] line instead). Gate verify_periodic.py
  PASS. PR #66 (branch periodic-house-treatment, worktree
  /tmp/opencode/wt-periodic-house).
- von_mises_policy: svg->pdf (cell c09 + module default), re-bake 6 s,
  0 errors, figure 100/300 two-res clean. Gate verify_policy.py PASS.
  pytest 106 passed on branch. PR #67 (branch policy-house-treatment,
  worktree /tmp/opencode/wt-policy).
- BOTH PRs also flip their README rows to the house validated wording.
- NOTE: /tmp/opencode/wt-periodic (pr/explain-periodic) is the user's
  explain-series worktree - do not touch; the new one is
  wt-periodic-house.
- After both merged: README intro micro-PR - replaced the stale
  "retain saved results ... not evidence that a fresh run succeeds"
  sentence with validated wording (PR #68 -> main 3df54e4).
  PERIODIC PR #66 -> main 6620ade; POLICY PR #67 -> main 565e81d
  (the two README-row edits conflicted on adjacent rows; resolved by
  merging main into the policy branch - keep both validated rows).
- PR #43 (feat/circular-rhat, unwrap.diagnostics joint circular R-hat)
  CLOSED as superseded: user confirmed the per-study house pattern
  (Divergences-first + circular_summary + plain ArviZ r_hat/ESS filters)
  is the diagnostic treatment and everything converged/mixed; the module
  was never merged and no notebook uses it.

## BEETLES STUDY PR - MERGED (PR #65 -> main 8a0539f; branch + worktree deleted both sides)
- 2 commits on the branch: (1) 28159df svg->pdf save paths (3 cells c08/c15/c16
  + the 3 module signature defaults; _save suffix router + svg fonttype rc
  untouched) + cleared the stale pre-upgrade outputs; (2) 0d42d14 the re-bake
  + README row (house validated wording) + fixed the stale "zebra and beetle
  read preprocessed tables" sentence -> zebra only.
- Re-bake: 77 s total, 0 errors, fit 4x(2000+2000)/accept 0.9 seeded in 52 s,
  Divergences: 0. SCIENCE BIT-IDENTICAL to the prior bake: kappa_mu
  intercept concentration 0.5927, and the per-draw kappa_mu_* summary table
  matches cell-for-cell (seeded determinism confirmed). PPC cross-check
  "response params from PyMC graph OK" passes. 4/4 figures 100 dpi inline +
  300 dpi PDF under figures/beetles/, no svg. pytest 106 passed locally.
- Gate: /tmp/opencode/verify_beetles.py (mirrors verify_bats.py; 4 figure
  cells {8,15,16,24}, EXPECT_SAVED stems, model cell {10} white-list, c02
  version dict, c04 set_style echo, c06 head columns, c10 fit+Divergences,
  c12/c13 r_hat, c19/c21 kappa_mu, c24 cross-check, stray pdf/svg audit).
  BEETLES GATE: PASS.
- PR #65 title: "beetles: land the study on the house template (E2E
  re-bake, two-res PDF figures)"; body /tmp/opencode/beetles_pr_body.md.
- Already largely templated on main (via #53/#54 + #55 stack):
  25 cells (setup c02 house-style + SEED=42; c06 E2E data cell =
  preprocess.load_beetle_model_table() from
  unwrap.datasets.load_beetle_security_light; c10 vM Bambi fit verbatim
  house style: seeded 2000+2000/accept 0.9, non-centered kappa REs,
  per-obs log-likelihood via
  pm.compute_log_likelihood(posthoc.graph_free_rvs(...)) (the
  bambi#1084 workaround is ALREADY in), "Divergences:" printed first);
  c12/c13 circular_summary diagnostics; 4 figure cells c08/c15/c16/c24.
- REMAINING: (1) figure save paths are .svg in c08/c15/c16 (c24 is
  .pdf) - align to the current house gates which check 300 dpi
  %-PDF- under figures/beetles/ (density/direction pattern); the
  module _save already routes by suffix, so it is a 3-line cell edit;
  (2) full re-bake on the current stack (the 24 baked outputs are from
  pre-upgrade bakes) + verify gate (write verify_beetles.py mirroring
  verify_locust_density.py: 4 figure cells, Divergences: 0, summary
  filters, PDF heads) + README row ("full rerun pending" -> "fully
  rerun and validated from this checkout (...)"). Much lighter than
  bats: one 2000+2000 fit (no mixture, no PPC-heavy section beyond the
  single PPC figure).
- Data facts: load_beetle_model_table() adds binary Light, wraps
  Heading deg->rad to (-pi, pi] via circular_modulo, category columns
  the model expects; E2E via the user-level ~/.cache/unwrap-data clone.

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
