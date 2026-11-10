# Agent instructions

Rules for AI agents and contributors working in this repository.

## Model priors and sampling specifications

Adjustments to a model's priors or to its sampling configuration (draws,
warmup, target acceptance rate / adapt delta, chains, cores, initial
values, seed) must be systematically tested and justified with documentation
before they land in a notebook or a PR:

- Test the full fit at the proposed specification and check the notebook's
  own hard gates (R-hat, effective sample size, divergences) on identifiable
  quantities. An untested specification change does not go in.
- Isolate changes: do not mix a prior change with a sampling-configuration
  change in one untested step. Change one thing, test it, then layer the
  next.
- Document the evidence for every change that lands: the specification
  tested, the seed, the observed R-hat / ESS / divergence numbers, and the
  verdict. Put it in the PR description, and where a value is load-bearing,
  a short factual note next to the constant in the affected notebook cell.
- Do not raise or drop sampling settings (acceptance target, draw budget)
  as a first response to a failing gate, and do not relax them by default.
  Change them only with the failing/passing numbers quoted and a stated
  reason.
- If a gate fails, stop and report the numbers (parameter, R-hat, ESS,
  divergences) before retrying under a different specification. One bake
  per tested change; no silent retry loops.
