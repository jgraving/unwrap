# Choose a circular model

Start with what you measured and what you want to learn. Suppose you recorded
an animal’s heading as it left an arena. You might ask where it tended to
point, how consistently it pointed there, or whether it used more than one
direction. The same dataset can support different questions.

![Schematic heading distributions lead to circular GLMs, circular GLMMs, and mixtures, with worked examples in Python and R. A second panel links to introductory guides and a prompt for coding agents.](assets/unwrap-model-map.svg)

[Open the linked model map](assets/unwrap-model-map.svg)

## Start from your data and question

| Your data and question | Starting model | Worked examples |
| --- | --- | --- |
| **Headings with predictors.** Does a treatment change the preferred heading, directional consistency, or both? | [Circular GLM](math/glms.md). Model mean direction and concentration as functions of predictors. | Simulated circular regression in [Python](notebooks/circular_regression.ipynb) or [R](notebooks_r/circular_regression.ipynb). |
| **Repeated headings within animals or groups.** Do animals differ in their preferred heading or in how they respond to a cue? | [Circular GLMM](math/glmms.md). Group effects account for shared differences among measurements of the same animal or group. Partial pooling lets the animals’ estimates inform one another. | Dung beetle orientation in [Python](notebooks/beetles.ipynb) or [R](notebooks_r/beetles.ipynb). |
| **Several plausible heading patterns.** Could headings reflect directed and undirected travel, or could there be two preferred directions? | [Circular mixture](math/mixtures.md). Give each component its own heading distribution and estimate how they combine. | Bat homing ([Python](notebooks/goldshtein_bats.ipynb), [R](notebooks_r/goldshtein_bats.ipynb)) for directed and uniform components. Locust steering ([Python](notebooks/locust_direction_mixture_model.ipynb), [R](notebooks_r/locust_direction_mixture_model.ipynb)) for competing directions. |

For a single set of headings with no predictors, start with a
[von Mises distribution](math/von-mises.md) to estimate mean direction and
concentration. If the thing you measure is a rate or a count that varies
with time of day, time is periodic but your outcome is not an angle. The
[daily activity example](examples.md#periodic-but-not-circular-a-24-hour-activity-rhythm)
shows how to put the periodicity in the predictor.

The likelihood and the grouping structure answer different questions.
You can add animal effects to a mixture when you observe each animal repeatedly.
Group effects capture shared variation within animals. Closely spaced steps
can additionally depend on preceding steps. Use the biological question and
predictive checks to decide whether another component is warranted.

## Start from what you know

- **I know regression and GLMs.** Start with [Angle basics](angle-basics.md).
  Then see how [circular and positive links](links.md) map linear predictors
  to mean direction and concentration.
- **I know circular statistics.** Start with [Bayesian basics](bayes-basics.md)
  to see what changes after you observe headings. Then read the
  [circular GLM guide](models.md#from-a-biological-question-to-a-glm) to see
  where predictors enter the model.
- **I’m new to statistical modeling.** Start with the
  [biological questions](index.md#start-with-a-biological-question), then
  work through [Angle basics](angle-basics.md) and
  [Bayesian basics](bayes-basics.md).
- **I’m using a coding agent.** Use the [copyable prompt](coding-agents.md#copyable-prompt)
  to describe your question and data, find a relevant example, and check the fit.

## Fit a first model

Run the simulated circular regression tutorial in
[Python](notebooks/circular_regression.ipynb) or
[R](notebooks_r/circular_regression.ipynb). It generates its own data and
shows how to fit the model, interpret its predictions, and compare those
predictions with the observed headings. Follow the
[Python setup](getting-started.md) or [R setup](getting-started-r.md)
instructions to run it locally.
