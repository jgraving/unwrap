## STATE - updated 2026-11-04 (R EXAMPLES PHASE: P0-P3 ALL MERGED (P3 @ 9b37395 incl. full-budget bake outputs + final markdown); WORDING PR #94 MERGED (8a636b6); PY c25 FORMULA FIX #95 MERGED (58f5a5b); ONLY WORK REMAINING = VECTORIZATION PR (post-P3, bats-only))

### CI DECISION (2026-11-03, user: "Do not run ci on GitHub for notebooks. They need to be locally tested only. Write tests but make it so that GitHub doesn't run them we run them locally.")
- `.github/workflows/notebooks.yml` is DELETED (it ran the Python tutorial + every R notebook on every PR, 240-min job).
- New `tests/test_notebooks.py`: executes every committed notebook (notebooks/*.ipynb python3 + notebooks_r/*.ipynb ir) end to end; asserts >=1 inline PNG per notebook. Skips under GITHUB_ACTIONS (hard guard) AND requires local opt-in UNWRAP_TEST_NOTEBOOKS=1 (unit suite stays fast; the GitHub tests workflow installs .[test] only, so nbclient is absent there too -> importorskip). UNWRAP_NOTEBOOKS=stem,stem subset; UNWRAP_NOTEBOOK_TIMEOUT per-cell seconds (default 14400); kernel-missing = skip.
- R budget flag renamed UNWRAP_R_CI -> UNWRAP_R_QUICK (+ ci_run -> quick_run) across all 4 R notebooks: no CI job anymore, it's a local two-chain smoke budget. Markdown in notebooks + docs/getting-started-r.md, docs/development.md (new "Notebook execution (local only)" section), docs/languages.md updated.
- Local battery = bake (write-back) then UNWRAP_TEST_NOTEBOOKS=1 pytest over the 4 R notebooks (+ Python subset as wanted). 16 cores available: fits already use cores=CHAINS; do NOT run multiple IRkernel bakes in parallel (shared /tmp/opencode/Rlib -> self-install cell race, see gotcha 15).

## R EXAMPLES PHASE (IN PROGRESS 2026-11-03; branch `r-examples-p0` — P0 tutorial baked, PR pending)

### Objective (user)
R versions of ALL example notebooks, side-by-side with Python on the
docs site. Governing requirement: **they need to work** (real baked
outputs + green CI). User: "Just do what you think is best … I won't
know until I see it"; follow precedent for combined R/Python packages.

### Format decision (settled; do not relitigate — revisited 2026-11-03 with user, stands)
Jupyter `.ipynb` + IRkernel (`ir`), NOT Quarto/RMarkdown. Reasons: the
site is mkdocs + mkdocs-jupyter, which renders ipynb natively (proven:
R syntax highlighting `hl-R`, strict build clean at /tmp/opencode/rprobe
before it was cleaned up); per-language Jupyter notes is the
marginaleffects/polars/arrow precedent; Quarto would rebuild the whole
site (nav, template, Colab JS, link validation). Side win: the Colab
button JS only matches `/notebooks/<name>.ipynb`, so R pages get no
Colab toolbar with zero template change.

### P0 deliverables (branch `r-examples-p0`, pushed)
- `notebooks_r/circular_regression.ipynb` — BAKED (real outputs; seed 42,
  4 chains x 1000 / 500 warmup, adapt_delta 0.9, backend cmdstanr,
  `unwrap_von_mises` custom family from `R/unwrap_family.R`) with 3
  INLINE PNGs (ragg) + 3 vector PDFs in `figures/circular_regression/`
  (gitignored; only the notebook is committed).
- `sync_docs_notebooks.py` generalized (notebooks/ and notebooks_r/ →
  docs mirror); `mkdocs.yml` nav pairs (Python/R per study: the tutorial
  pair is in); `docs/getting-started-r.md` (R setup page); `.gitignore`
  (`docs/notebooks_r/`).
- `.github/workflows/notebooks.yml` new `execute-r` job (45 min):
  setup-r → pip `cmake` insurance → pak 0.11
  `repo_add(universe = stan-dev r-universe, ppm = noble) +
  pkg_install(cmdstanr, brms, posterior, loo, circular, ggplot2, ragg,
  IRkernel)` → `cmdstanr::install_cmdstan()` → `IRkernel::installspec()`
  → python nbclient `kernel_name="ir"` with `UNWRAP_R_CI=1` (2x800/400)
  → uploads executed notebook artifact. First real test of the install
  recipe happens on the PR.
- Authoring + tooling (scratch, not repo content): /tmp/opencode/rnb_build.py
  (17-cell source of truth — PATCH IT + rerun before every bake),
  bake_r.py (full budget, writes back), dry_ci.py (CI shape, copy, no
  write-back), verify_r.py (battery), mini_bake.py/parse_cells (verbatim
  cell parse/execute loop — how the trailing-comma bug was found),
  smoke_all.R (Rscript driver, cached fit ~minutes).
- Equivalence vs Python tutorial (baked, seed 42, n=200, mu=wrap(x),
  kappa=softplus(2+0.5x); truth 0/1/2/0.5 pre-softplus):
  |param|post mean|truth|rhat|ESS|
  |b_Intercept|-0.075|0|1.0001|2333|
  |b_x|1.055|1|1.0041|2307|
  |b_kappa_Intercept|2.278|2|1.0048|2329|
  |b_kappa_x|0.598|0.5|1.0009|2072|
  (R-hat/ESS on the 2D [iteration x chain] draws — a flat 2000-vector
  is a pseudo-single-chain and UNDERSTATES R-hat: gave 0.9995/0.9998/
  0.9996/1.0002 instead of the honest 2D values above.)
  All within ~1 SD; divergences: not retained in this brms/cmdstanr
  build (cell probes fit$fit@sim$samples and reports honestly).

### 2026 R-ecosystem gotchas (each cost real time; keep when writing R cells)
1. `cmdstanr` is OFF CRAN — Stan team r-universe
   https://stan-dev.r-universe.dev; install w/ BOTH repos so CRAN deps
   resolve (end-user snippet lives in setup cell + R setup page).
2. pak 0.11 API: `repo_add(...)` + `pkg_install(...)` (NOT
   `pak::pak(repos=)` / `set_repos`). PPM repo gives binaries.
3. RcppParallel configure imports the python `cmake` MODULE →
   `python3 -m pip install cmake --break-system-packages` (CI has the
   insurance line).
4. cmdstanr 0.9: `find_cmdstan_dir()` gone → `cmdstan_path()` (probe
   wrapped in tryCatch in the setup cell).
5. circular 0.5.2 sampler revamp: `rvonmises(n, mu, kappa,
   control.circular = list())`; mu AND kappa length 1 ONLY (no
   broadcasting); plain-numeric mu → coercion warning per call; the
   warning-free shape: `vmu <- function(m) circular::circular(m,
   units="radians", template="none", modulo="2pi")` +
   `vm_ctrl <- list(units="radians", template="none", modulo="2pi")` →
   default [0, 2pi) support, wrap with `mod_circular` (R/unwrap_functions.R)
   to the model's (-pi, pi] convention. `modulo="pi"` FOLDS support to
   [0, pi] — do not use. The PPC cell is a 0.5.x-forced draw-by-observation
   nested loop (2000x200 ≈ 2-4 min full budget).
6. R TRAILING COMMA: `f(a, b,)` parses (in call position) and forwards
   the empty slot into `...` → `vapply` dies with "unused argument
   (alist())"; in bare `(a, b,)` position it is a parse error. Scan
   pattern: line ending `,` + next line starting `)`.
7. Unbraced R if/else: `else` on its own line after an unbraced if is a
   parse error — brace the block.
8. brms 2.23 brmfit is an S3 LIST — no @diagnostics/@rescor slots; this
   brms+cmdstanr combo does NOT retain the per-iteration divergence flag
   anywhere in the fit (verified by recursive walk) → probe + honest
   report; R-hat/ESS carry convergence.
9. Custom-family dpar coefficients are named `b_kappa_Intercept` /
   `b_kappa_x` (brms infix `b_<dpar>_<coef>`), not `kappa_*`; prior
   matching via class= + dpar= verified against the fit (all 4 match).
10. posterior 1.7: the `ess()` generic is no longer exported →
    `posterior::ess_bulk()` (rhat still exported).
11. `as_draws_df` is WIDE (one column per variable; no `value`); long
    trace format must be built with base R (do.call rbind).
    `nrow(vector)` is NULL — vector arithmetic is not a matrix: index
    draws directly in loops.
12. Headless inline figures: R built without cairo/X → IRkernel's
    default png device can't rasterise; `ragg::agg_png` +
    `IRdisplay::display_png` (display the file; keep `dev.off()`
    invisible or IRkernel text-dumps its return value). ragg is a
    PPM binary dep — add wherever the stack is installed.
13. ggplot2 4.x SORTS line data by the x aesthetic (verified on 4.0.3:
    a polyline with data order x = 0,1,0.5,2 is built as 0,0.5,1,2).
    Old-ggplot2 "lines connect data rows in order" is dead. So:
    trace-style plots must use x = iteration (monotonic within a line;
    harmless under the sort) with one line per chain via colour/group;
    "x = value, y = chain-level" traces render as spaghetti walls.
    Related classic-R trap that amplified the mess: `outer(X, Y, FUN)`
    passes the WHOLE X vector to FUN (not scalars) — a custom FUN
    closing over a longer vector returns a vector that is silently
    recycled into the result array. Use the bare `outer(X, Y)`
    elementwise product and build [grid x draws] explicitly.
14. Rootless R detail (machine-local, not repo): PPM-binary .so files
    carry NEEDED libR.so → R core must be built `--enable-R-shlib`;
    non-PIC static libbz2.a in the rtoolchain dir broke the libR.so
    link → fix was re-pointing the dangling `libbz2.so` symlink at the
    system `.so.1.0`. The ir kernelspec carries machine-local env
    (R_LIBS_USER + LD_LIBRARY_PATH → /tmp/opencode/R...).
 15. Do NOT run multiple IRkernel notebook bakes in parallel on this
    machine: the ir kernelspec points every kernel at the SAME
    /tmp/opencode/Rlib, and the notebooks' self-install cell then
    races (observed: two parallel arrow source builds; one kernel left
    the lib with a broken bit64 -> "there is no package called
    'arrow'"). Bake serially (one bake_p1.py call with all names), and
    pre-install any missing dep into Rlib in one process first
    (pak with the noble PPM repo has arrow/binaries). The parallelism
    the box offers lives INSIDE a fit (brm cores = CHAINS, 4 NUTS
    chains concurrent) - that is where the 16 cores get used.
 16. Installing extra R packages with pak + PPM binaries: the URL
    archaeology is a trap. Probed/DEAD ends this session:
    packagemanager.posit.co /ubuntu/noble/latest (404), /cran/__latest__
    (+ bin/linux/ubuntu/noble/contrib/4.x — 404), /cran/__linux__/noble/
    latest (200 index but SOURCE-only: zero File: entries), pak
    repo_add("ppm" = "noble") (treated as hostname), /cran/2024-04-01
    (no noble R-4.5 binaries for the packages needed). The CI recipe in
    .github/workflows/notebooks.yml (repo_add ppm = "…/cran/__linux__/
    <codename>/latest") predates our R-4.5 rootless machine. WORKING
    ROUTE (no root in this container, sudo blocked by no-new-privileges):
    CRAN SOURCE + userland toolchain: `apt-get download gfortran-13
    gcc-13 cpp-13 libgcc-13-dev libgfortran-13-dev libgfortran5
    libquadmath0 libudunits2-dev libudunits2-0 libabsl-dev && dpkg -x
    *.deb /tmp/opencode/rtoolchain/`; gfortran-13 then needs
    -fno-use-linker-plugin (no liblto_plugin) plus the gcc-13 internal
    libs from the same debs; Makeconf now has FC/FFLAGS/FCLIBS lines
    (backup Makeconf.bak-noFC). Makeconf LDFLAGS already carries
    rtoolchain rpath. Then plain install.packages(c(...)) from
    cloud.r-project.org with PATH=+rtoolchain/usr/bin and
    PKG_CONFIG_PATH=+rtoolchain/usr/lib/x86_64-linux-gnu/pkgconfig.
 17. A bake driver's silence is NOT death: bake_p1.py prints
    "===== BAKE <name> =====" once per notebook and only prints
    BAKE-WRITE-OK/FAIL when the whole notebook finishes; nbclient
    streams cell output nowhere in between. This session restarted a
    "dead" locust bake that was actually running fine (both wrote the
    same .ipynb; killed the duplicate by explicit PID). Check the log
    tail + notebook mtime + pgrep before assuming a bake died.

### Process notes
- Co-agent (other PC) was briefed on the R bake but went offline before
  starting; the prompt I wrote for him is STALE (pre cmdstanr-off-CRAN,
  pre pak 0.11) — do not reuse. Codex (cloud, code-only, no R) takes
  P1: draft the 3 mechanical ports (beetles, locust density, zebra) on
  branch `r-examples-p1` off `r-examples-p0`, UNBAKED notebooks, draft
  PR targeting r-examples-p0; I (main agent) run every cell locally and
  fix. Prompt ready in this chat (2026-11-03); user will ship it at P0
  green.
- Bake loop that worked: patch rnb_build.py → rerun it (regenerates
  pristine notebook) → parse-scan (trailing commas) → mini_bake/smoke
  (verbatim cell text, keep.source=TRUE) → bake_r.py full → verify_r.py
  → dry_ci.py (UNWRAP_R_CI=1, copy, no write-back) → strict docs build
  (NO_MKDOCS_2_WARNING=1, sync then mkdocs build --strict) → YAML check
  → commit (ADDITIVE only — branch was shared) → push.

### P0 PR #88 CI fixes (2026-11-03, three recipe bugs found via runner logs)
1. `r-lib/actions/setup-r@v4` does not exist → `@v2`.
2. `setup-r` installs a BARE R: bootstrap pak itself (PPM binary first,
   CRAN source fallback); the 2026 runner image lacks the curl/ssl
   dev headers (pak's `curl` dep) → apt-install
   `libcurl4-openssl-dev libssl-dev` as insurance.
3. ubuntu-latest MIGRATED OFF NOBLE in 2026 (runner now reports the
   `resolute` codename): a hardcoded `__linux__/noble` PPM URL
   silently loses its binaries → resolve `VERSION_CODENAME` from
   /etc/os-release in the recipe.
4. The pak recipe's `repo_add(..., ppm = ..., )` TRAILING COMMA: R
   forwards the empty argument into pak's `...` spec list →
   "argument is missing, with no default". (Base R tolerates trailing
   commas in `install.packages(...)` because it never reads the dots;
   pak's repo_add does. Scan R blocks in YAML/markdown too — the
   notebook-only parse scan misses them.)
Codex (P1) state: branch `r-examples-p1` off 181fb47, draft PR #89
(its nav/CI-loop commit 28daaff is now MOOT: the main-merge into the
local branch superseded it, and the CI loop itself was deleted per the
user's local-only decision). Static review PASSED: priors are EXACT
matches to the Python bambi models (von Mises coefficient vectors via
custom Stan `target +=` contributions + period bounds +
LogNormal(log 30, .5) kmu hyperpriors); helper cells solid (true-
minimum HDI, scaled-Bessel A1, P0-sampler-convention `simulate_angles`,
`save_figure` = ggsave + show_figure); trace/polar cells use
ggplot2-4.x-safe orientations; data loads from the same
jgraving/unwrap-data repo (sparse clone into repo_root/unwrap-data,
gitignored).

## P1 merge sequence (2026-11-03, after the CI decision)
Branch `r-examples-p1-local` in /tmp/opencode/wt-r-p1:
- `1f52b31` merge origin/main (conflicts resolved: notebooks.yml,
  getting-started-r.md; P0 notebook in tree == main's plot-fixed
  version; nav = Python/R pairs).
- `f66b639` lb/ub numeric period literals (brms verbatim-source
  injection gotcha).
- `eca5c68` CI DROPPED: notebooks.yml deleted, tests/test_notebooks.py
  added (GITHUB_ACTIONS hard-skip + UNWRAP_TEST_NOTEBOOKS opt-in),
  UNWRAP_R_CI -> UNWRAP_R_QUICK rename + docs (development/languages/
  getting-started-r) repointed at the local test.
- BAKE 1 (bake_chain.log) finished ALL-BAKES-OK including locust:
  **locust_density_order PASSES the gate clean** (0 error cells, 15
  PNGs; screening c14 + final c23 tables all rhat ~1.00 / ESS >= 600;
  gate c24 empty). Beetles + zebra FAILED the gate on the angular
  coefficients -> diagnosed, fixed, re-baking (see below).
  NOTE: the first chain LOOKED dead (nbclient buffers; the driver only
  prints per-notebook at completion) - do NOT restart a bake while the
  original is still running (a duplicate was running and killed by
  PID once the original finished; always check the log + file mtimes).
- WRAP-CLASS PATHOLOGY (beetles + zebra) - root cause + fix:
  bounded (-pi,pi] angular components + likelihood pinning their SUMS
  mod 2pi + per-individual heading means on BOTH sides of the +/-pi
  seam (beetles: OFF circ means +2.27 for some individuals, -2.33 for
  others) => many 2pi-equivalent posterior "wrap classes" separated by
  ~20-nat barriers. cmdstanr default init (all chains -2.4) and even a
  +-0.6 "anchor" jitter (crosses the seam) put chains in DIFFERENT
  classes; each chain converges inside its class but sine R-hat/ESS
  across chains fail (1.2-2.1 / ESS 1-15 on all 22 angular params;
  Python twin passes only because its 4 chains land in one class).
  Priors/model/data verified identical to Python (generated Stan clean;
  196-row parity; kmu posteriors 2.23/3.31 vs Python 2.225/3.29; Light
  1.254 vs 1.26; all 20 REs matching). VERIFIED FIX: same-wrap-class
  seeding - per-coefficient anchors from circular means (intercept =
  base-condition mean; slope = circ-mean difference; per-individual =
  difference of condition means), mod_circular wrap + runif(-0.08,0.08)
  jitter, set.seed(SEED+k) per chain, kmu at 30. Probe (4x4000/2000,
  notebook's own diagnose_fit): ALL clear (1.0004-1.0021 / ESS
  1167-1647). Patched into beetles c12 + zebra c11 fit cells (+ comment
  + stopifnot guards). Python twin intercept row (2.08 sd 1.04) is the
  one non-matching quantity: weakly-identified 23-D wrap-ridge transient
  spread; document in notebook, do not chase.
- ZEBRA VERDICT (2026-11-03): the 8000-iter re-bake (bake_zb3) FAILED
  IDENTICALLY to the 3000-iter one (71 failing rows: M17 bias0/bias1
  rhat ~1.19/1.20 ESS 5.5/3.9, kmu rhat 1.1859/1.1645 ESS 4.0/4.4,
  family-wide 1.01-1.09; kmu means 9.442 -> 9.437 after 2.7x draws).
  ESS FLAT in N = bottleneck-trapped, not under-sampled. The budget was
  the wrong lever. The zebra .ipynb on disk currently holds this FAILED
  bake - never commit it; c03 budget (8000/3000) + md c10 note edits are
  kept and stay in the final patch. Next: (A) ORIGIN-class init probe
  running now (/tmp/opencode/diag_zo.R -> diag_zo.log, 4x8000/3000, all
  angulars ~0 +/-0.25 jitter + kmu 30, full budget - mirrors PyMC's
  effective init and tests whether the data-anchored per-track init was
  placing chains in a bad region; live cells in /tmp/opencode/zebra_live/,
  fit cell truncated before the brm call, diag cell truncated before the
  posterior-extract call). If (A) passes -> patch zebra c11 to origin
  init, keep bounded priors as written, re-bake ONLY zebra
  (bake_p1.py zebra_escape), check_gate. If (A) fails -> (B) UNBOUNDED
  angular coefficients: exact twin geometry - delete the three
  uniform(-pi(), pi()) bound lines in c09, KEEP every
  set_prior("target += unwrap_von_mises...") line unchanged (re-verified
  1:1 vs the twin dict: beta0[1] VM(0, 1e-6), beta0[2] VM(0, 2),
  bias0 vect VM(0, kmu_intercept), bias1 vect VM(0, kmu_slope), kmu
  LogNormal(log30, 0.5) - so the faithful port needs NO new prior lines;
  the user's "small kappa like 1e-5" rule applies to any coefficient
  lacking an explicit VM prior, and here every one has one). Origin
  init (+/-0.25, kmu 30) replaces the data-anchor init (drop seam
  guards - no bounds to straddle); rewrite the md c10 wrap paragraph to
  "unconstrained reals + periodic VM priors (Python circ_mod geometry),
  chains converge in the same mode because the init is shared".
  Periodic VM on an unbounded real makes the total posterior an
  improper countable mixture of identical proper classes; NUTS stays in
  a single mode per init and every quantity of interest is circular
  (mod 2*pi), so all circular outputs are well-defined - exactly PyMC's
  geometry (its rhat 1.00). (C) nutpieR - see NUTPIE STATE below.
  (D) nimbleHMC last resort.
  STATUS (2026-11-03, later): (B) promoted to PRIMARY (the bounded box
  was a deviation from the twin in the first place: PyMC circ_mod has
  no +-pi bounds). UNBOUNDED-STANCODE VERIFIED via
  /tmp/opencode/zebra_live/unbounded_stancode.R: parameters block =
  vector[K_beta0/K_bias0/K_bias1] UNBOUNDED; prior block has NO default
  Gaussian on the angulars (only the six target += VM/lognormal lines +
  kappa/sd lines); zero uniform lines. So (B) = exact twin geometry.
  Probe (B) RUNNING as /tmp/opencode/diag_zu.R -> diag_zu.log
  (4x8000/3000, origin init +/-0.25, kmu 30; first launch died on my own
  grepl bug - escaped brackets under fixed=TRUE - fixed, relaunched).
  NOTEBOOK PATCH PRE-WRITTEN and syntax-checked (R parse OK, 0 trailing
  commas, anchor asserts): /tmp/opencode/patch_zebra_unbounded.py -
  whole-cell replacement of c09 (drop the three uniform bound lines,
  rewrite comment), c10 (md: circ_mod geometry + per-class convergence
  + sine R-hat invariance), c11 (stan check asserts UNBOUNDED vectors +
  no uniform_lpdf; origin init mk_init, n_track guard; same brm call).
  Apply ONLY after the diag_zu probe passes the gate, then re-bake
  zebra via bake_p1.py zebra_escape and check_gate.py.
- ZEBRA VERDICT FINAL (2026-11-03): diag_zu probe PASSED -
  FAILING ROWS 0 of 80, max rhat 1.0003, min ESS 4570 (kmu now ESS in
  the thousands; angular per-chain wrapped means agree to 0.005-0.01;
  old worst M17 was rhat 1.20/ESS 5.5). Patch APPLIED to the
  notebook. FIRST POST-PATCH BAKE FAILED WITH c13 "object 'fit' not
  found" - ROOT CAUSE FOUND (NOT model, NOT kernel): my patch wrote
  the cell sources as splitlines(keepends=False) LISTS; nbformat
  joins list source lines with "" (EMPTY string), so each patched
  cell flattened to ONE line whose first line is a COMMENT = the
  entire cell silently commented out (no-op, no error, 0.0s: c09
  defined nothing, c11 fit never ran). This is why exec showed
  c02-c13 with exec_counts but c13 missing fit, and why Rscript
  probes with the same code worked. Fixed with
  /tmp/opencode/fix_flatten.py: re-joined lines, stored sources as
  PLAIN STRINGS (nbformat passes those verbatim), verified
  byte-identical to intended text + R-parse of the actual on-disk
  joined form + 0 trailing commas. RULE FOREVER: write .ipynb cell
  sources as plain strings (or lists with keepends=True), NEVER
  splitlines() without newlines; and when validating patched cells,
  R-parse the text AFTER the on-disk join/serialization, not the
  in-memory intended string. GOTCHA 18 (2nd full-bake death): the patch
  had also set execution_count/outputs on MARKDOWN cell c10 - forbidden
  additional properties; the bake's final nbformat.validate() failed
  AFTER the real fit completed (50 min lost, nothing written). Fixed
  (removed both keys) + ALWAYS nbformat.validate the file before
  launching a bake. Third full bake = bake_zb4b.log... (relaunch) nbclient 0.11 gotchas while debugging:
  no client.initialize(); client.setup_kernel() is a contextmanager;
  kc.execute() sync wrapper returns the msg-id string (probe via a
  temp nbformat code cell through execute_cell instead).
- DEFERRED OPTIMIZATION (user OK'd "later", 2026-11-03): R notebook
  post-fit stage is single-threaded base R (PPC/diagnostic loops +
  ragg figures) vs numpy in the twins; sampling itself is already 4
  parallel C++ cmdstan processes. Candidate follow-up PR: vectorize
  the PPC/diagnostic loops (results-identical, faster), and/or brms
  future.apply where the pattern fits - DO NOT mix into the faithful
  port PRs. Could be a codex brief.
- USER RULE (angular priors): NEVER uniform(-pi,pi) bounds/flat
  overrides for circular coefficients - use periodic VON MISSES with
  small kappa (1e-5-ish) as the weak baseline. (Zebra's faithful port
  already satisfies this: 1e-6/2/kmu-based VM targets = the twin dict.)
- NUTPIE STATE: no Linux binaries on andytimm.r-universe (bin/linux
  listing empty; r-4.5-x86_64 precompiled dir 404s; R errors on
  type="binary" on Linux). Source build needs rustc >= 1.85: rustup
  installed in $HOME/.rustup, toolchain+install chain running in
  background -> /tmp/opencode/rustup_install.log (first attempt failed
  on RUSTUP_HOME mismatch between the rustup install and the R script
  env; fixed with consistent exported env). Keep as option-C insurance
  ONLY if (B) fails; using nutpie requires hand-written Stan (nutpie
  takes a path to a .stan file, not a brms object) + notebook surgery
  (nutpie draws -> manual draws array; re-implement every downstream
  brmfit call). API: nutpie_sample(model, data, num_draws, num_warmup,
  num_chains, cores, seed, target_accept), nutpie_diagnostics.
- IN FLIGHT (resolved -> see ZEBRA VERDICT): the 8000-iter re-bake
  (bake_zb3.log) finished and FAILED the gate identically; 4 chains,
  8000 iter, 3000 warmup (5000 retained). The 3000-iter bake failed the
  gate with UNDER-MIXING, not class splits: rhat 1.00-1.21 family-wide,
  kmu_intercept/slope ESS 3.7/5.0 (still descending the LogNormal prior
  ridge, mean 9.4 vs prior 30), failing angulars ESS 5-120. Beetles
  re-bake is already CLEAN on disk (init fix + md annotation + baked
  outputs: rhat 1.0000-1.0009, ESS 849-1944, gate empty, 9 PNGs) -
  commit pending alongside zebra. Zebra md c10 documents the bumped
  budget; c03 carries ITER/WARMUP 8000/3000.
  NOTE: the first 8000-iter re-bake was killed at 41 min on a SUSPECTED
  zombie (log frozen at "===== BAKE ... =====", no R processes visible)
  shortly after the user's chat crash - but the "no R processes" check
  was FAULTY (pgrep without -f matches the 15-char comm name, which is
  just "R"; the real cmdline is
  /tmp/opencode/R/.../exec/R --slave -e IRkernel::main() ...), so the
  old fit may have been healthy and slow (2.7x the previous budget).
  Corrected liveness rule: `pgrep -af "IRkernel::main"` for the kernel
  AND `ps -C R` / top for a `model_8+` cmdstan sampling process at
  high CPU. The log only prints per-notebook at completion (bake
  silence is normal - gotcha 17), so log silence must NEVER be the
  kill criterion alone. Restarted as bake_zb3.log; deterministic
  (SEED 42 + data-anchored init) so the result is identical, only the
  wall time is redone.
- GATE CHECKER (use this, the old inline checker was bogus):
  /tmp/opencode/check_gate.py <ipynb...> - R print() wraps wide tables
  into TWO row blocks (block 1 'idx name mean sd rhat', block 2 'idx
  ess_bulk'); parsing only toks[-2]/toks[-1] of block-1 lines mis-read
  sd/rhat as rhat/ess AND silently 'passed' (the earlier lstrip/isdigit
  variant skipped every row). The script re-joins the blocks; failing
  row = rhat >= 1.01 or ess <= 400. Validated: reproduces the known zebra
  bad rows (M17 1.2117/5.6, kmu 1.2041/3.7) exactly and reads the clean
  beetles/locust notebooks as 0 failures.
- PR #89 body DRAFT: /tmp/opencode/pr89-body.md (local-test narrative;
  fill the zebra figure count at final verification).
- REMAINING: verify re-bakes clean (gate tables empty) + commit ->
  local test run (UNWRAP_TEST_NOTEBOOKS=1, 4 R notebooks, + Python
  tutorial if quick) -> push origin r-examples-p1 -> REOPEN #89 (gh
  api -X REOPEN; it was auto-CLOSED when r-examples-p0 was deleted) ->
  PATCH base to main (closed PR base PATCH 422s; do it after REOPEN) ->
  PR body update (local-test narrative, drop "CI" claims) -> merge
  (plain) once tests.yml + docs.yml green (no notebook job anymore) ->
  delete branch both sides + worktree.
- P2: codex PR #90 is UP (r-examples-p2 @ 1ce0a98, draft, base main).
  FULL REVIEW in /tmp/opencode/review-p2.md: 89/89 cells parse, exact
  prior matches (periodic c08; bats c41 + c53 incl. exact uniform mix
  component), UNWRAP_R_QUICK + hard gates, local-only wording. The
  earlier 'LOO bug' is RETRACTED (false alarm, re-verified on the actual
  c48 code: loo::loo gets RAW ll; the per-observation shift feeds only
  relative_eff, which is invariant to it) - the cell is correct, do not
  'fix' it; at bake time just sanity-check delta-ELPD sign/ballpark vs
  the Python committed value. The codex fix prompt was SENT as a
  fresh-session brief (codex's chat crashed): only the same-wrap-class
  init for c39/c44 (reference code embedded) + the explicit LOO
  retraction + budget-vs-init triage note. CODEX PUSHED 5f3eb3d (verified
  clean, 2026-11-03): 1 commit, c39+c44 only; 62/62 cells parse, 0
  trailing commas; anchors = reference-treatment circ mean + wrapped
  offsets, per-tag = tag mean - ref anchor - own-treat offset; kmu_tag 30;
  stopifnot count+seam guards; AND each cell checks its generated Stan
  parameters block (bounded vector names/dims) + make_standata design
  order BEFORE sampling. PR #90 STATUS section updated. Next when P1
  merges: merge main into r-examples-p2 (mkdocs.yml conflict expected),
  bake the 3 notebooks serially (deps pre-installed), verify gates +
  delta-ELPD ballpark vs Python, then re-target #90 to main. P2
  deps readxl/sf/geosphere/patchwork are PRE-INSTALLED into Rlib
  (userland gfortran route, see gotcha 16). Expect mkdocs.yml text
  conflict when merging main post-P1; map = static sf fallback
  (documented deviation, accept).
- P2 BAKE ROOT-CAUSES (2026-11-03, both fixed in r-examples-p2 @ 815c945):
  (A) ENDPOINT-GUARD LOST IN TRANIT. periodic_not_circular c11 guarded the
      simulated Beta response with 1 - .Machine$double.eps, but brms
      serializes the standata JSON at ~15 significant figures, so the max
      Y (1 - 2^-52) ROUNDED TO 1.0 in the data file; Stan's beta_lpdf
      then returned -inf at EVERY initial value (all inits rejected,
      "Log probability evaluates to log(0)"), and brm surfaced only
      "Fitting failed. Unable to retrieve the metadata". R-side checks
      never caught it (the round-trip only happens in the JSON). FIX:
      cap at 1 - 1e-10 (round-trips exactly, ~10 orders below any real
      signal). LESSON: any data cell passed to Stan must guard/round to a
      margin that survives ~15-sig-fig JSON serialization, not just the
      in-memory .Machine$double.eps.
  (B) draws_array [ METHOD DOES NOT DROP. posterior 1.7.1 as_draws_array()
      keeps class draws_array; its [.draws_array forces drop=FALSE, so
      as.matrix(a[,,p]) yields ONE long column (draws x chains flattened,
      e.g. 20000x1), NOT a draws x chains matrix. The P1 & P2 diagnose_fit
      helpers are WRITTEN FOR THAT flattened layout (rhat/ESS consume it).
      P2's copies added a spurious stopifnot(ncol(x)>=2L) that fails
      because the object is 1-column by design. LESSON: never assert a
      multi-column shape on as.matrix(draws_array[,,p]); keep the
      finiteness check only. (subset_draws keeps the class too; the raw
      brm fit RDS in /tmp/opencode/periodic_fit_small.rds reproduces it.)
  CLI NOTE (cmdstan 2.40): the compiled-model argument grammar changed -
  legacy `--key=value` now ERRORS; use bare `key=value` (e.g.
  `data file=x.json method=sample num_warmup=...`), or just go through
  cmdstanr/brm which emit the new form. My direct-binary diagnosis that
  "every data file fails" was this CLI change, not a broken build.
- P2 STATE (2026-11-03, post-bake): worktree /tmp/opencode/wt-r-p2b,
  branch r-examples-p2-bake. ALL THREE NOTEBOOKS BAKED (bake_p2.log:
  ALL-BAKES-OK) + gate-verified (check_gate.py, extended for stream
  outputs - P1/P2 print() tables land in stream, not text/plain):
  periodic rhat 1.0002-1.0024 ESS 1991-2882 (6/6); bats baseline
  rhat 0.9999-1.0019 ESS 1016-9368 (144), mixture rhat 0.9999-1.0029
  ESS 668-9322 (216), PNG 3/5/39, 0 error cells; bats delta-ELPD vs
  Python reference (+118.763 total): R = +118.488 total, per-treatment
  rows within 0.2, pareto_k_max 0.50. Commits on origin/r-examples-p2:
  0923f17 (merge main, mkdocs pair-style), 377211f (endpoint guard),
  815c945 (ncol guard), a3e347b (bakes), b0ccff7 (von_mises_policy R
  markdown rewritten in the paper's nav-vonmises feature-box language -
  utility/gradient/goal-heading/policy derivation display equations +
  section notes in paper symbols; verified mkdocs build --strict locally
  and in rendered HTML, all equations -> arithmatex). PYTHON twin mirror
  MERGED to main as e1e5539 (PR #93, docs/vm-policy-md): same rewrite with
  Python code-name mapping (goal_heading, vonmises_pdf, Z = 2*pi*I0);
  checks all green; worktree + branch cleaned up. NOTE: origin/main is
  now e1e5539 (P2 branch merged origin/main @ 0923f17 predates the
  e1e5539 md change in the PYTHON notebook only - r-examples-p2 does NOT
  touch notebooks/von_mises_policy.ipynb, so #90 has no conflict risk;
  re-verify only if a build fails). PR #90: base main verified,
  body UPDATED via REST (this gh build's `gh pr edit` FAILS on the
  deprecated GraphQL projects mutation; workaround: build
  {"body": ...} -> gh api --input file -X PATCH
  repos/jgraving/unwrap/pulls/90 with GITHUB_REPOSITORY env).
  NEXT: local acceptance pytest (3 P2 notebooks, dir-qualified
  UNWRAP_NOTEBOOKS, log /tmp/opencode/test_p2.log) -> when green:
  gh pr ready 90 -> checks -> gh pr merge 90 --merge --delete-branch ->
  clean up wt-r-p2b + remote branch -> tell user "P3 ready - send now"
  (prompt /tmp/opencode/codex-p3-brief.md verbatim).
   P2 MERGED (2026-11-04): local acceptance GREEN (test_p2.log: 3 passed,
   12 skipped in 8998s - bats + periodic + von_mises_policy in-memory,
   full budget, data root ~/.cache/unwrap-data). PR #90 READY 11-04,
   all checks green (pr90_checks.log: 16 pass / 0 fail / deploy skip),
   MERGED to main as 53e0a3c with --delete-branch; main now has all 7
   R notebooks (origin/main = 53e0a3c, history: 53e0a3c merge #90 <-
   e1e5539 merge #93 <- 3fd38e8 P1). Cleanup DONE: wt-r-p2b removed,
   r-examples-p2-bake branch deleted, remote r-examples-p2 deleted
   (verified 404). REMAINING: (1) fold this file into the repo per the
   runbook - DONE: folded to the log branch
   docs/notebooks-to-scripts-handoff as handoff-r-phase.md
   (fold commit 266fedf, byte-identical verified with cmp);
   (2) tell the user the exact phrase "P3 ready - send now" and ship
   codex-p3-brief.md VERBATIM (it assumes all P0-P2 templates on main -
   now TRUE). DISCOVERY (2026-11-04, pre-phrase): DRAFT PR #91 ALREADY
   EXISTS - "Draft P3 R port: locust direction-mixture study", head
   r-examples-p3 @ 9fdfcb8 (unbaked
   notebooks_r/locust_direction_mixture_model.ipynb, 1 commit), base
   main, branched from main @ 15abf88 (P0-only main, 24 commits STALE:
   predates P1+P2). Codex evidently ran an earlier brief on 11-03.
   When the user shipped codex-p3-brief.md (user said SHIPPED
   2026-11-04; codex is working): codex must CONTINUE from
   PR #91 / r-examples-p3 and REBASE it onto current origin/main
   (merge main first, as P2 did), NOT start a fresh branch.
   P3 ACCEPTANCE PLAN (main agent):
   - Draft audit (done 2026-11-04): the 9fdfcb8 draft is a SUBSTANTIVE
     46-cell port (custom_family locust_mix, both fits, chunked
     explicit_loo, per-individual ELPD, PPC, hypothesis + alignment
     tables). It predates P1/P2 templates -> recheck: same-wrap-class
     init per the bats c39/c44 pattern, diagnose_fit gate at BOTH
     budgets, UNWRAP_R_QUICK switch, nav pair, mk_init/stopifnot
     guards, no trailing commas.
   - Python REFERENCE re-running in background: locust_ref_run.py ->
     locust_ref.log (venv editable unwrap re-pointed to
     /tmp/opencode/wt-mainref origin/main worktree - the old
     /tmp/opencode/wt-unified mapping was dangling after cleanup and
     `import unwrap` died; fixed with pip install -e). Captures:
     gate tables (failing rows) BOTH models, az.loo per model
     (elpd/se/pareto_k_max) + loo_compare, key posterior means.
     CAVEAT IN ADVANCE: the COMMITTED Python unimodal output ALREADY
     has 19 failing gate rows (worst mu_pop rhat 1.05, ESS 110; the
     bimodal committed table is EMPTY). If the R unimodal bake
     fails the strict gate the same way, the levers are (a)
     data-anchored same-wrap-class init on mu_pop/mu_raw, (b) budget
     bump like zebra 4x8000/3000, (c) document-vs-twin deviation as
     last resort - decide after seeing the ref-run numbers.
   - Bake: bake_p3.py written (BAKE_ROOT=wt-r-p3, 6h per-cell budget,
     same validate protocol); deps confirmed in Rlib (circular ggplot2
     ragg IRdisplay patchwork posterior brms cmdstanr loo readxl all
     TRUE). Budget convention: UNWRAP_R_QUICK 2x3000/1500 else
     4x4000/2000 (like bats) - verify the cell mirrors it.
   - After green: same finish sequence (pr ready -> checks -> merge
     --delete-branch -> worktree/branch cleanup), then hand off the
     P3 acceptance evidence and note the phase (P0-P3) is COMPLETE.
   P3 REVIEW (2026-11-04, codex shipped 6da2d63 = d8f4b01 merge-main +
   442a697 notebook + 6da2d63 nav): PASSED static. Worktree
   /tmp/opencode/wt-r-p3 (branch r-examples-p3-review).
   - Model parity vs Python twin: EXACT (both models, every prior; the
     PR STATUS has the full port tables + the sort-encoding explanation,
     satisfying the brief's acceptance).
   - Sort/label-switching encoding: brms custom_family
     vars="first_component" repurposed as a PARAMETER-DEPENDENT ORDER -
     tparameters `real first_component = b_popone[1] <= b_poptwo[1] ? 1
     : 2` (linear order, per brief), likelihood locust_mix_lpdf selects
     unsorted components by the saved index, flips the animal eta sign
     when index 2 leads, log_sum_exp of the two log-weighted VM
     densities with locust_wrap(mean sums). extract_state reads the
     saved first_component from the posterior draws and re-sorts every
     component-indexed quantity per draw (stopifnot vs p1>p2). Verified
     correct.
   - Static probe (p3_stan_probe.py -> p3_stan_probe.log, cells
     c04-c24 verbatim, NO sampling): 0 errors; generated Stan clean
     (one first_component decl, bounded b_popone/b_poptwo/b_devone/
     b_devtwo [1,1,52,52], kmu_one/two lower=0, all target += lines,
     design order 0..51 checked); data parity on the REAL data: 54
     recordings / 52 animals / 2435 steps / 3-204 per animal (exact
     Python counts, stopifnot-enforced in c13).
   - Init: bimodal = kmeans(2, nstart=20) on (cos,sin) cluster anchors
     + per-animal offsets, seam-margin stopifnot, +/-0.08 post-wrap
     jitter, seed SEED+k per chain (same-wrap-class pattern).
     Unimodal = pooled circ_mean + seam-gap re-anchor if any animal
     lands near +/-pi (codex found the pooled mean was seam-close).
   - Budgets: quick 2x4000/2000, FULL 4x8000/4000 (codex bumped above
     bats' 4x4000/2000 on purpose: Python unimodal has 18 marginal
     gate rows - see locust_reference_elpd.txt). PPC/hyp/align
     constants mirror locust_posthoc exactly (1000/1000 boot, 360
     grid, 16 bins, bw 30/15, seeds 123 + 10000+idx); LOO = the bats
     pattern (raw ll -> loo::loo, scaled relative_eff); prints
     ELPD table + delta table for acceptance.
   - Nav pair commit verified (mkdocs.yml, same pair style).
   - Full-budget BAKE IN FLIGHT (bake_p3.log, driver bake_p3.py,
     per-cell timeout raised to 10 h; the bimodal fit cell is the long
     pole, expect hours). Machine note: a llama-server at ~430% CPU
     shares the box (load ~7-8) - wall times stretched but bakes
     proceed.
   - ACCEPTANCE when green: gate tables (bimodal clean strict;
     unimodal may sit on its marginal rows - reference file has the
     criteria + levers), ELPD vs reference (bimodal -3404.1 SE 39.2,
     uni -3501.8 SE 38.7, delta ~+97.6 bimodal better, uni k-hat
     marginal EXPECTED), posterior parity targets (mu_pop ~ +/-1.45,
     mu_kappa ~27/21, weights_pop ~0.5) - all in
     /tmp/opencode/locust_reference_elpd.txt.
   WORDING PR SHIPPED (2026-11-04, PR #94, branch
   docs/wording-paper-language @ a536456, base 53e0a3c, worktree
   /tmp/opencode/wt-mdpolish): markdown-only pass over all 7 Python +
   7 R example notebooks (the R locust-mixture twin EXCLUDED - its
   wording ships in P3). Audit = full md dump of all 15 notebooks
   (md_audit_all.txt) read end to end + the paper's five Example
   sections mined as the language standard (beetles ~1272, zebra ~1359,
   locust order-density ~1443, bats ~1620, locust mixture ~1740 in
   paper.tex). docs/*.md pages scanned: already clean, untouched.
   Changes: (1) dev meta removed (CI/test-harness refs, clone/cache
   plumbing, bambi internals + bug workaround, "In the implementation"
   mapping tables, "The centerpiece.", PyTensor line, "Local review",
   headless-map-fallback, "Notebook execution is local by design" in
   all R Setup cells); (2) stale draft wording ("This draft is
   intentionally unbaked" x3, "unbaked drafts make no convergence
   claim", R bats "(baseline, Bambi)"/"(PyMC)" headings + "section
   heading identifies the Python reference" notes); (3) R twins stand
   alone ("This R port of X.ipynb" -> model description); (4) slop/
   redundancy (locust-mixture TOC + duplicated hypothesis re-derivation
   c26, label-restating PPC bullets c19, doubled weak-ID story in
   locust_density c34, "security-light trial table", "Rmd's
   grouped-data machinery"); (5) Figure-style save_path=/gitignored
   plumbing -> one-liner in all 8 Python notebooks (vmp included);
   (6) FIXED PRE-EXISTING SITE BUG: locust-mixture c13 display math
   was not blank-line-separated -> arithmatex never processed it (4
   stray $$ on the page) + unescaped set braces in Theta_i/Y_i/Y.
   Local acceptance done: nbformat.validate all, full strict docs
   build (== docs CI job) passes, rendered HTML page-by-page check
   clean (only the 2 MathJax-config $$ literals on every page; 404.html
   has 0, expected). Net diff: 15 files (14 + vmp py), -183 lines.
   Patch tooling: /tmp/opencode/apply_mdpolish.py (exact-string,
   one-match assertions, skip-argv for the already-applied files).
   PR #94 OPEN, checks pending (docs job is the meaningful gate).
   No re-bake (no executed cell changed). When green: plain merge
   --delete-branch per house flow.
   VECTORIZATION PR SCOPED (2026-11-04, plan in
   /tmp/opencode/vectorization-plan.md): only two changes - bats c29
   kde_vm (361-iteration vapply -> matrix colMeans, algebraically
   identical, ~100x) and bats c81 sim_one_tag (per-row scalar
   rvonmises(1,..) vapply -> two rvonmises(n, scalar, ..) calls;
   PROBED: n-draws call == scalar stream identical, 371x faster for
   n=2000; new y vector same distribution, different stream -> PPC
   bands shift by MC noise, disclosed in PR; gates/LOO untouched).
   Everything else (c48/c86/c90, zebra PPC) already vectorized or
   negligible. PR = bats notebook only, base = post-P3 main, rebake
   bats only. Queue it AFTER P3 merges (machine serializes bakes).
   (3) deferred: R post-fit vectorization follow-up PR.
  NOTE: bats c48 LOO cell untouched; reference values in
  /tmp/opencode/bats_reference_elpd.txt.
  P3: locust direction — custom
  2-component circular MIXTURE Stan model, the one real modeling risk
  (validate vs the Python posterior).

## CLOSE-OUT (2026-11-04 session 2: P3 landed, wording landed, P3 R markdown pass done)
  ALL THREE PRS MERGED, main now 9b37395. P0-P3 of the R examples
  phase are COMPLETE; nothing P-related is open or queued except the
  vectorization PR (see plan above; P3 merge unblocks it - machine is
  free) and the deferred R post-fit slowness PR.
  PR #94 (wording) MERGED @ 8a636b6: markdown-only, 15 files
  (7 py + 7 R + py vmp figure block), +127/-269; raw-JSON patch
  runner (/tmp/opencode/apply_mdpolish.py) preserved each file's own
  serialization (vmp pair = ensure_ascii=True no trailing newline;
  all others ensure_ascii=False + "\n") so the diff is exactly the
  edited markdown cells (0 artifacts); strict docs build + page-by-
  page HTML check green; CI full matrix green.
  PR #91 (P3 R port) MERGED @ 9b37395: head 1f4eb72 = branch 6da2d63
  (static-reviewed model) + baked full-budget outputs + final markdown
  pass (/tmp/opencode/apply_p3r_polish.py, 27 patches: TOC/draft-notes/
  dev-meta removed, R twin stands alone, budgets documented in c19,
  arithmatex c23 display-math fix + escaped braces, Python variable
  names + vestigial "Formula correction" note dropped, paper-aligned
  result sentence in c50). Acceptance (bake = acceptance run; test
  suite is a redundant re-execute subset): 51 cells 0 errors; BOTH
  fits strict-clean (0 failing rows, the in-notebook stopifnot gate);
  ELPD bimodal -3404.488 (SE 39.242) / unimodal -3501.249 (SE 38.763),
  delta +96.761 (SE 16.710) vs Python ref -3404.123/-3501.771 (+97.6);
  honest caveat recorded in PR body: bimodal max Pareto k 0.7041 (1/2435)
  0.004 over 0.70 (Python 0.676) = PSIS boundary noise, comparison
  unaffected; population hypothesis: lateral (±pi/2) elpd -1.285 vs
  forward 0 elpd -2.567, P(lateral>forward)=1 -> vectorial /
  ring-attractor supported.
  PR #95 MERGED @ 58f5a5b: fixed the Python locust-mixture c25
  hypothesis-test display (shared Z over two different-concentration
  components is not the mixture normalizer) -> per-component
  M(y|mu,kappa) form the code actually evaluates (matches c26 + R c42);
  same cell drops the `weights_pop`-style variable names. No executed
  cell changed; strict docs build + HTML arithmatex check green.
  GH QUIRK (cost time twice): `gh api` POST/delete on
  repos/{owner}/{repo}/pulls/N/merge intermittently 404s when cwd is
  inside a worktree; retry from a cwd outside the repo (e.g. /tmp)
  with GITHUB_REPOSITORY set, or use `gh pr merge` from the worktree
  that holds the head branch. `gh pr merge --delete-branch` also
  fails when another worktree has main checked out ("already used by
  worktree") - plain merge + git push --delete works fine.
  WORKTREES: wt-mdpolish, wt-c25, wt-r-p3 removed; local branches
  docs/wording-paper-language, docs/locust-hypot-formula,
  r-examples-p3-review deleted; remote branches r-examples-p3,
  docs/locust-hypot-formula, docs/wording-paper-language gone.
  KEPT: wt-handoff (this file), wt-mainref (detached 53e0a3c - the
  venv's editable unwrap target), house repo /home/jake/projects/
  unwrap on docs/notebooks-explained (READ-ONLY rule unchanged).
  NEXT: vectorization PR per /tmp/opencode/vectorization-plan.md -
  base = current main (9b37395), bats-only: c29 kde_vm matrix
  (~100x), c81 sim_one_tag rvonmises n-draws (371x, stream-identical
  probe done), PPC-bands MC-noise disclosure, full bats rebake
  (bake_p3.py driver reusable with BAKE_ROOT + name; liveness =
  pgrep "IRkernel::main"), check_gate.py vs bats_reference_elpd.txt +
  gate tables + delta-ELPD +118.49.
