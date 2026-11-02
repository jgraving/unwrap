# QA scratch for the notebooks-to-scripts task

Per-study verification gates, re-bake drivers, smoke drivers and
notebook builders used while landing the study PRs. These are working
scaffolding (hardcoded slim-cell maps), NOT shipped with the study
PRs - the GitHub mirror here exists purely so the toolset survives an
offline/loss event. The live copies run from /tmp/opencode in the jake
environment; this mirror is refreshed whenever the handoff doc
(../TASK-unwrap-notebooks-to-scripts.md) is synced.

Gates are worktree-aware: run from a study worktree's cwd, or pass the
root as argv[1] (`python verify_*.py /path/to/worktree`). Each gate
checks the baked notebook: zero error outputs, the version dict, the
study-specific content (fit lines, `Divergences: 0`, best-model /
summary / cross-check lines), and that every saved figure file exists
and is a valid vector file.

- verify_bats.py — goldshtein_bats (the big one; `BATS GATE`)
- verify_beetles.py — beetles (`BEETLES GATE`)
- verify_circular.py — circular_regression (`CIRCULAR STUDY GATE`;
  3 saved SVGs incl. the PlotCollection trace)
- verify_zebra.py — zebra_escape (`ZEBRA GATE`)
- verify_locust_dirfix.py — locust_direction_mixture_model (`LOCUST GATE`)
- verify_locust_density.py — locust_density_order (`LOCUST DENSITY GATE`)
- verify_periodic.py — periodic_not_circular (`PERIODIC GATE`)
- verify_policy.py — von_mises_policy (`POLICY GATE`)
- rerun_nb.py — full re-bake driver (venv jupyter kernel, generous
  timeouts); prints `DONE in Ns: images=K saved-figures=K errors=0`
- smoke_bats_nb.py — reduced-size E2E smoke (2 tags/treatment, 200/200)
- build_bats_nb.py — regenerates the slim bats notebook from the
  annotated base (base backup restorable via
  `git show main:notebooks/goldshtein_bats.ipynb`)
- build_bats_plots.py — regenerates goldshtein_bats_plots.py
