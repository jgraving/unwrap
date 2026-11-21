# Use unwrap with a coding agent

Give the agent your research question, a data dictionary, and the unit of
observation. For headings, say how you measured the angles and what counts
as zero. If an animal appears in several rows, include its identifier and any
trial, time, or trajectory columns. These details help the agent choose a
model and interpret its predictions.

Open your analysis project in your coding agent. Copy the prompt below and
replace the bracketed fields. Give it a local file path or a synthetic sample
with the same columns and structure. Start from the
[model chooser](choose-a-model.md), the [worked examples](examples.md),
and the [unwrap source repository](https://github.com/jgraving/unwrap).

## Copyable prompt

```text
Help me answer this question with unwrap (jgraving/unwrap).

Scientific question: [What do I want to estimate or compare?]
Data: [Local file path or a representative synthetic sample]
Unit of observation: [What does one row represent?]
Response: [Column, degrees or radians, zero direction, positive rotation]
Predictors: [Columns and the conditions or values I want to compare]
Sampling design: [Animal and group IDs, repeated trials, time or step order]
Environment: [Python with Bambi/PyMC or R with brms/Stan, and my setup]
Deliverable: [A notebook or script, plus the plots and summaries I need]

1. Read unwrap's model chooser (docs/choose-a-model.md) and the closest
   worked example in notebooks/ (Python) or notebooks_r/ (R).
   Name the example and explain how its sampling design matches mine.
   Ask me about missing details only when they affect the model or the answer.

2. Check whether the outcome is a direction, an unoriented axis, or a
   measurement with a periodic predictor. Confirm units, angular reference,
   wraparound, and the amount of data per condition and per group. State any
   assumptions and preserve the original data.

3. Choose the simplest model that answers my question. Say whether the
   predictors should affect mean direction, concentration, or both. If there
   are repeated measurements, decide which group effects and any time
   dependence the design calls for. Assess separately whether the scientific
   question or predictive checks motivate a mixture of heading patterns.

4. Make the analysis reproducible. Record the environment, seed, data
   transformations, model formula, and priors. Check what the prior predicts,
   then fit the model. Report sampling diagnostics such as R-hat, effective
   sample size, and divergences, and compare posterior predictions with the
   observations. Use a labeled synthetic dataset to exercise the workflow if
   I have only provided a schema. Say which checks and fits actually ran.

5. Report predictions or contrasts that answer my question and include
   posterior uncertainty. Use circular differences for angles, and explain
   directional consistency separately from precision of the estimated mean.
   Flag limitations from the sampling design, fit, or predictive checks, and
   give the commands to rerun the analysis.
```
