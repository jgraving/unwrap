# QA scratch for the notebooks-to-scripts task

Per-study verification gates, re-bake drivers, smoke drivers and
notebook builders used while landing the study PRs. These are working
scaffolding (absolute /tmp/opencode paths, hardcoded slim-cell maps),
NOT shipped with the study PRs - the GitHub mirror here exists purely
so the toolset survives an offline/loss event. The live copies run from
/tmp/opencode in the jake environment; this mirror is refreshed whenever
the handoff doc (../TASK-unwrap-notebooks-to-scripts.md) is synced.

- verify_locust_density.py / verify_bats.py — per-bake gates (house
  pattern; expect "LOCUST DENSITY GATE: PASS" / "BATS GATE: PASS")
- rerun_nb.py — full re-bake driver (venv kernel, 16200 s/cell)
- smoke_bats_nb.py — reduced-size E2E smoke (2 tags/treatment, 200/200)
- build_bats_nb.py — regenerates the slim bats notebook from the
  annotated base (base backup restorable via
  `git show main:notebooks/goldshtein_bats.ipynb`)
- build_bats_plots.py — regenerates goldshtein_bats_plots.py
