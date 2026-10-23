%% Circular Dataset Structures (Matlab only)
% Matlab version of "R/Circular Dataset Structures.Rmd": simulated circular
% datasets that provide examples of effects on mean direction, concentration,
% and individual differences. Fitted circular models are used to estimate
% these effects and compare hypotheses about their causes.
%
% This uses the same strategy as MatlabOnly.m: the log posterior of the unwrap
% model (circular modulo link for $\mu$, softplus link for $\kappa$) is written
% out explicitly and sampled with Hamiltonian Monte Carlo (|hmcSampler|,
% Statistics and Machine Learning Toolbox). No Stan installation is needed.
%
% Differences from the R version:
%
% * gradients are analytic (finite differences are too slow for the
%   hierarchical models)
% * leave-one-out cross-validation uses Pareto-smoothed importance sampling
%   (PSIS-LOO), with exact refits for observations with Pareto $k > 0.7$
%   (equivalent to |loo(..., reloo = TRUE)|, instead of moment matching)
% * posterior densities of mean vectors are shown as clouds of posterior draws
% * Rhat is the classic (not rank-normalised) version used in MatlabOnly.m;
%   summary tables show it both before (|Rhat_raw|) and after unwrapping
%   circular parameters (|Rhat|)
%
% *Inputs:* |circular_datasets/*.csv|, exported by "R/Export Circular Datasets.R",
% so that all versions model identical simulated data, and the functions in
% |unwrap_functions/|.
%
% Set |saveFigures = true| before running to save figures to |figures_circular/|.
% Set |quickTest = true| for a fast check with short chains.

if ~exist('saveFigures', 'var'); saveFigures = false; end
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'unwrap_functions'));
dataDir = fullfile(scriptDir, 'circular_datasets');
figDir = fullfile(scriptDir, 'figures_circular');
if saveFigures && ~exist(figDir, 'dir'); mkdir(figDir); end

rng(120810506) % ISBN Batschelet 1981

% Sampler settings (as for brms defaults: 4 chains, 1000 warmup, 1000 draws)
opts = struct('numChains', 4, ...
              'numSamples', 1000, ...
              'burnin', 1000, ...
              'targetAccept', 0.65, ...
              'reloo', true, ... % exact refits for observations with high Pareto k
              'maxRefits', 10); % if more observations have high k, report them instead
if exist('quickTest', 'var') && quickTest
    opts.numSamples = 200; opts.burnin = 200; opts.reloo = false;
end

% Colours
col.kappa = hex2rgb('#1E78B5'); % kappa parameter
col.rho = hex2rgb('#F08024'); % mean vector length
col.sd = hex2rgb('#E74A29'); % SD and mean angle
col.sd2 = hex2rgb('#E57461'); % other SD heuristics
col.obs = hex2rgb('#3E1F51'); % control observations
col.treat = hex2rgb('#006400'); % treatment observations

%% Divergence from home direction
%
% <<../R/icons/ant.png>>
%
% _Desert ants of the genus Cataglyphis are known for their use of path
% integration (dead reckoning) to return to their nest. Disruption of visual
% orientation cues can result in deviations from the expected homewards
% direction._
%
% In navigation experiments, return directions are often compared with the
% expected home direction to determine whether homeward navigation is
% perturbed. Here, the example data ($n = 20$) were generated from a
% distribution with a mean that differs from the expected direction (0°) by
% only 15° ($\kappa = 10$).
%
% A common short-cut is the v-test modification of the Rayleigh test, which
% tests uniformity against clustering in some direction less than 90° from the
% expected direction. The data are very far from uniform and the mean angle is
% close to the expected angle, so the test is significant. Even so, the data
% show a trend away from the expected angle, with 75% of points falling to the
% left of 0°.

dtDiv = readtable(fullfile(dataDir, 'divergence.csv'));
yDiv = dtDiv.y;
fprintf('Proportion left of 0°: %.2f\n', mean(wrapPi(yDiv) < 0));
[vStat, vP] = vTest(yDiv, 0);
fprintf('v-test (expected direction 0°): V = %.3f, p = %.2g (normal approximation)\n', vStat, vP);

%% Fit the model
% Each of $\mu$ and $\kappa$ requires a single parameter estimate (an
% intercept). We expect a mean direction anywhere on the circle and high
% concentration ($\kappa > 1$), so highest prior density is set in that range:
% $\mu \sim N(0, 2\pi)$ and, on the inverse-softplus scale,
% $\kappa \sim N(3, 3)$. Since the mean direction is far from the wrap-point
% ($\pm\pi$), Rhat for $\mu$ can be interpreted on the linear scale.

mDiv = makeModel(yDiv, ...
    'Xmu', ones(numel(yDiv), 1), 'muNames', {'b_Intercept'}, ...
    'muPrior', [1 0 2*pi], ... % normal(0, 2pi)
    'Xk', ones(numel(yDiv), 1), 'kNames', {'b_kappa_Intercept'}, ...
    'kPrior', [3 3]); % normal(3, 3)
fitDiv = fitUnwrapModel(mDiv, opts, 'model_divergence');
printSummary(fitDiv);

%%
% The posterior estimates are summarised as the median mean vector (red
% arrow), compared with the true mean vector (dark purple arrow), and a cloud
% of posterior mean vectors.

figure('Name', 'Divergence model', 'Color', 'w');
tiledlayout(1, 2);
ax = nexttile; circBase(ax, 0);
circPoints(ax, yDiv, col.obs);
drawCloud(ax, getDraw(fitDiv, 'b_Intercept'), softplus(getDraw(fitDiv, 'b_kappa_Intercept')));
circArrow(ax, deg2rad(-15), A1(10), col.obs, 5);
circArrow(ax, circMedian(getDraw(fitDiv, 'b_Intercept')), ...
    A1(softplus(median(getDraw(fitDiv, 'b_kappa_Intercept')))), col.sd, 2);
ax = nexttile;
vertHist(ax, rad2deg(unwrapCircular(getDraw(fitDiv, 'b_Intercept'))), ...
    'mean angle', [-30 15], col.sd, 0);
saveFig(saveFigures, figDir, 'divergence_model');
fprintf('mu left of expected angle: %.1f%% of estimates\n', ...
    100*mean(wrapPi(getDraw(fitDiv, 'b_Intercept')) < 0));

%%
% The fitted model indicates that the mean direction is to the left of the
% expected direction in $>99\%$ of $\mu$ estimates.

%% Test the expected angle hypothesis
% To test the hypothesis that the true mean direction falls in the expected
% direction, we fit another model in which the mean direction is fixed at the
% expected direction (no $\mu$ parameter), and compare the two models.

mExpect = makeModel(yDiv, ...
    'Xk', ones(numel(yDiv), 1), 'kNames', {'b_kappa_Intercept'}, 'kPrior', [3 3]);
fitExpect = fitUnwrapModel(mExpect, opts, 'model_expect');
looDiv = looUnwrap(fitDiv, opts);
looExpect = looUnwrap(fitExpect, opts);
looCompare({'model_divergence', 'model_expect'}, {looDiv, looExpect});
fprintf('ELPD difference (divergence - expected): %s\n', elpdDiff(looDiv, looExpect));

%%
% The model for a divergent mean direction has higher expected log predictive
% density (ELPD) than the model for the expected direction. This ELPD
% difference is $>2\times$ larger than its standard error, so we can state
% with some confidence that the hypothesis of a divergence from the expected
% direction is more likely than the hypothesis of alignment with the expected
% direction.
%
% If the expected direction is not 0°, we align the data with the expected
% angle, then fit both models (expected direction and divergent direction) to
% the aligned data. Here is an example with the expected angle at -15°.

yAligned = yDiv - deg2rad(-15);
fitTrue = fitUnwrapModel(setY(mExpect, yAligned), opts, 'model_true');
fitFalse = fitUnwrapModel(setY(mDiv, yAligned), opts, 'model_false');
sFalse = printSummary(fitFalse, false);
fprintf('model_false mu estimate (deg; mean, SD, 2.5%%, 97.5%%): %s\n', ...
    mat2str(round(rad2deg(table2array(sFalse(1, {'Mean', 'SD', 'Q2_5', 'Q97_5'}))))));
looTrue = looUnwrap(fitTrue, opts);
looFalse = looUnwrap(fitFalse, opts);
looCompare({'model_true', 'model_false'}, {looTrue, looFalse});
fprintf('ELPD difference (expected at -15° - divergent): %s\n', elpdDiff(looTrue, looFalse));

%%
% Now the expected direction model has higher predictive power, even though
% both find a similar mean direction. If the divergence model had higher
% predictive power, but within one standard error of the expected direction
% model, we could not rule out the possibility that the expected direction is
% the true direction, though evidence would be weaker.

%% Change in direction
%
% <<../R/icons/sandhopper.png>>
%
% _Sandhoppers (Talitridae) forage along beaches, following the movement of
% the sea to avoid drying out. Their movement is guided by celestial cues,
% including the moon._
%
% If the expected direction is not known, it can be estimated from control
% data and compared with treatment data. In this example, the treatment
% induces a rightward turn of 30°, changing $\mu$ from 180° to 210°, with
% intermediate concentration ($\kappa = 3$).
%
% The treatment indicator $x$ is used as a predictor of both $\mu$ and
% $\kappa$. Predictors are not centred, so the intercepts describe the control
% condition (equivalent to |0 + Intercept| in |brms|), and the |x|
% coefficients describe the change caused by the treatment. For stability, a
% slightly narrower prior on the intercept, $N(0, 3\pi/2)$, is used.

dtDelta = readtable(fullfile(dataDir, 'delta.csv'));
nD = height(dtDelta);
mDelta = makeModel(dtDelta.y, ...
    'Xmu', [ones(nD, 1) dtDelta.x], 'muNames', {'b_Intercept', 'b_x'}, ...
    'muPrior', [1 0 3*pi/2; 1 0 pi/2], ...
    'Xk', [ones(nD, 1) dtDelta.x], 'kNames', {'b_kappa_Intercept', 'b_kappa_x'}, ...
    'kPrior', [3 3; 0 3]);
fitDelta = fitUnwrapModel(mDelta, opts, 'model_delta');
printSummary(fitDelta);

%%
% Because the mean direction for the control falls at the boundary of
% $(-\pi, \pi)$, estimates for the intercept can appear in both positive and
% negative directions in different chains (|Rhat_raw|). After unwrapping, the
% chains agree (|Rhat|).
%
% Predictions for the control use |b_Intercept| and |b_kappa_Intercept|, and
% for the treatment |b_Intercept + b_x| and |b_kappa_Intercept + b_kappa_x|.

ctrl = dtDelta.x == 0;
plotTreatment(fitDelta, dtDelta.y(ctrl), dtDelta.y(~ctrl), ...
    deg2rad([180 210]), [3 3], col, 'Change in direction');
bx = getDraw(fitDelta, 'b_x');
nexttile;
vertHist(gca, rad2deg(wrapPi(bx)), 'change in mean angle', [-90 90], col.sd, 0);
saveFig(saveFigures, figDir, 'change_in_direction');
ciTurnDelta = quantile(rad2deg(unwrapCircular(bx)), [0.025 0.975]);
fprintf('Rightwards turn: %.1f%% of estimates\n', 100*mean(wrapPi(bx) > 0));
fprintf('Median turn: %.0f°; 95%% CI: %.0f° to %.0f°\n', ...
    rad2deg(circMedian(bx)), ciTurnDelta);

%%
% Uncertainty in the example data results in a wide range of posterior
% estimates for the treatment direction. Nonetheless, most estimates for
% change in $\mu$ (percentage above) are to the right of 0°, suggesting a
% rightwards change in direction.

%% Test the direction change hypothesis
% We fit models with treatment effects on $\kappa$ only, and with no effects
% of treatment, and compare them with the full model.

mKappaDelta = makeModel(dtDelta.y, ...
    'Xmu', ones(nD, 1), 'muNames', {'b_Intercept'}, 'muPrior', [1 0 3*pi/2], ...
    'Xk', [ones(nD, 1) dtDelta.x], 'kNames', {'b_kappa_Intercept', 'b_kappa_x'}, ...
    'kPrior', [3 3; 0 3]);
mNoDelta = makeModel(dtDelta.y, ...
    'Xmu', ones(nD, 1), 'muNames', {'b_Intercept'}, 'muPrior', [1 0 3*pi/2], ...
    'Xk', ones(nD, 1), 'kNames', {'b_kappa_Intercept'}, 'kPrior', [3 3]);
fitKappaDelta = fitUnwrapModel(mKappaDelta, opts, 'model_kappadelta');
fitNoDelta = fitUnwrapModel(mNoDelta, opts, 'model_nodelta');
printSummary(fitKappaDelta);
printSummary(fitNoDelta);
looDelta = looUnwrap(fitDelta, opts);
looKappaDelta = looUnwrap(fitKappaDelta, opts);
looNoDelta = looUnwrap(fitNoDelta, opts);
looCompare({'model_delta', 'model_kappadelta', 'model_nodelta'}, ...
    {looDelta, looKappaDelta, looNoDelta});
fprintf('ELPD difference (delta - no delta): %s\n', elpdDiff(looDelta, looNoDelta));
fprintf('ELPD difference (kappa delta - no delta): %s\n', elpdDiff(looKappaDelta, looNoDelta));

%%
% In this case there is uncertainty about the effect of treatment on
% direction. Although most estimates from the model with a change in direction
% indicate a rightwards turn (percentage and 95% CI above), these differences
% could also be explained well by a null model with no change (ELPD difference
% within its standard error). Since the median estimate is close to the true
% value of 30°, the model's estimates appear to be robust, if not conclusive.
% Conversely, there is no evidence for an effect of treatment on
% concentration, since the model with an effect of treatment on $\kappa$ only
% has lower predictive density than the null model, although this difference
% is small.

%% Change in concentration
%
% <<../R/icons/honeybee.png>>
%
% _Honeybees performing a "waggle dance" on a horizontal observation hive that
% view a green light-spot, orient their dances as if viewing the sun. The
% concentration of a series of dance directions decreases when the brightness
% of this light spot is reduced._
%
% An experimental treatment may also affect concentration independently of
% mean direction. In this example the treatment has exactly the same $\mu$
% (0°) as the control, but $\kappa$ decreases by $>80\%$ (from 3 to 0.5).
% We can use the same model as for the change in direction.

dtConc = readtable(fullfile(dataDir, 'conc.csv'));
mConc = setY(mDelta, dtConc.y);
mConc.Xmu(:, 2) = dtConc.x; mConc.Xk(:, 2) = dtConc.x;
fitConc = fitUnwrapModel(mConc, opts, 'model_conc');
printSummary(fitConc);
fprintf('mu estimates to the left of 0°: %.0f%%\n', 100*mean(wrapPi(getDraw(fitConc, 'b_x')) < 0));

%%
% The model successfully recovers a large change in $\kappa$ caused by the
% treatment, alongside only a negligible potential change in $\mu$ (similar
% posterior density to the left and right of 0°). Since $\kappa$ is estimated
% on the inverse-softplus scale, estimates for control and treatment are
% softplus transformed before calculating their difference.

ctrl = dtConc.x == 0;
plotTreatment(fitConc, dtConc.y(ctrl), dtConc.y(~ctrl), [0 0], [3 0.5], col, ...
    'Change in concentration');
bkx = getDraw(fitConc, 'b_kappa_x');
bk0 = getDraw(fitConc, 'b_kappa_Intercept');
nexttile;
vertHist(gca, softplus(bk0 + bkx) - softplus(bk0), 'change in kappa', [-5 1], col.kappa, 0);
saveFig(saveFigures, figDir, 'change_in_concentration');
fprintf('Reduction in kappa: %.1f%% of estimates\n', 100*mean(bkx < 0));

%%
% Although there is some uncertainty in the particular size of the change in
% $\kappa$, nearly all posterior estimates (percentage above) support a
% treatment effect that decreases $\kappa$.

%% Test the concentration change hypothesis

mConcKappa = setY(mKappaDelta, dtConc.y); mConcKappa.Xk(:, 2) = dtConc.x;
mConcNoKappa = setY(mNoDelta, dtConc.y);
fitConcKappa = fitUnwrapModel(mConcKappa, opts, 'model_conc_kappa');
fitConcNoKappa = fitUnwrapModel(mConcNoKappa, opts, 'model_conc_nokappa');
printSummary(fitConcKappa);
printSummary(fitConcNoKappa);
looConc = looUnwrap(fitConc, opts);
looConcKappa = looUnwrap(fitConcKappa, opts);
looConcNoKappa = looUnwrap(fitConcNoKappa, opts);
looCompare({'model_conc', 'model_conc_kappa', 'model_conc_nokappa'}, ...
    {looConc, looConcKappa, looConcNoKappa});
fprintf('ELPD difference (kappa only - no kappa): %s\n', elpdDiff(looConcKappa, looConcNoKappa));
fprintf('ELPD difference (mu & kappa - kappa only): %s\n', elpdDiff(looConc, looConcKappa));

%%
% In this case, model comparison is less conclusive than the posterior
% estimates. The model with an effect of treatment on $\kappa$ only predicts
% the data better than the model with no effects of treatment, but this
% difference is small relative to its standard error. With this sample size,
% cross-validation therefore provides only weak support for the change in
% concentration indicated by the posterior estimates. As we might have also
% expected, the effect of treatment on $\mu$, present in the original model
% but not the $\kappa$-only model, does not improve predictive density, so
% there is little evidence that the treatment affects $\mu$.

%% High inter-individual correlation
%
% <<../R/icons/moth.png>>
%
% _Tethered Bogong moths can be induced to fly in "flight simulators", which
% allows repeated measurements of flight directions across multiple
% individuals simultaneously. Individual directional biases may differ
% somewhat, but are typically highly correlated, aligning with the
% population's migratory direction._
%
% In datasets that include repeated measurements across individuals, it is
% necessary to account for individual differences before attempting to
% identify any population level trends. In the simplest case, the orientation
% behaviour of different individuals is very similar: all individuals share a
% common goal direction $\beta_0$, but differ in alignment with that direction
% by their individual offsets $\beta_i$. The distribution across $\beta_i$ is
% determined by the hyperparameter $\kappa_\mu$, the concentration across
% individual offsets.
%
% In this example, each individual's headings are highly concentrated
% ($\kappa = 5$) around its own bias $\mu_i = \beta_0 + \beta_i$, and the
% individual biases are moderately concentrated ($\kappa_\mu = 2$).
%
% The hierarchical model estimates $\beta_0$ as a single estimate across
% individuals and $\beta_i$ as an offset for each individual, with the prior
% $\beta_i \sim vM(0, \kappa_\mu)$ and, as expected for high correlation,
% $\log \kappa_\mu \sim N(\log 15, 0.6)$. Individual differences in $\kappa$
% follow a normal distribution (inverse-softplus scale) with standard
% deviation $\sigma_\kappa \sim HalfT(3, 0, 0.5)$, using a non-centred
% parameterisation. Circular parameters are summarised after unwrapping:
% for some individuals, estimates from different chains are separated around
% the circle (|Rhat_raw|), which can make it appear that the model has not
% converged.

dtHigh = readData(dataDir, 'highcorr.csv');
trueHigh = readtable(fullfile(dataDir, 'highcorr_true.csv'));
mHigh = makeIndividualModel(dtHigh.y, dtHigh.ID, ...
    'beta0Prior', [1 0 3*pi/2], 'kappaPrior', [3 3], 'sdScale', 0.5);
fitHigh = fitUnwrapModel(mHigh, opts, 'model_highcorr');
printSummary(fitHigh);

%%
% The plots show each individual's data (dark purple points) with its true
% mean vector (dark purple arrow), the posterior median mean vector (red arrow)
% and a cloud of posterior mean vectors, followed by the individual biases
% (orange points) and the population concentration $\kappa$.

plotIndividuals(fitHigh, dtHigh, deg2rad(trueHigh.mu), 5.0*ones(height(trueHigh), 1), ...
    0, 2.0, 5.0, col, 'High correlation');
saveFig(saveFigures, figDir, 'high_correlation');
fprintf('softplus(b_kappa_Intercept) > true kappa: %.1f%%\n', ...
    100*mean(softplus(getDraw(fitHigh, 'b_kappa_Intercept')) > 5.0));

%%
% The estimates of $\kappa_\mu$ have a slight positive bias, due to a
% combination of low sample size and our prior to higher values. The model
% nonetheless recovers good estimates of $\beta_0$ and $\beta_i$, as well as
% population $\kappa$, the concentration of the average individual, which
% needs to be distinguished from individual differences and the concentration
% across individuals.

%% Low inter-individual correlation
%
% <<../R/icons/beetle.png>>
%
% _When rolling a dung ball, ball-rolling dung beetles can be repeatedly
% returned to the centre of a circular arena on exiting, resulting in repeated
% measurements of exit direction. While a single individual exits
% consistently in the same direction, there is typically no relationship
% between individual directional biases._
%
% In this example each individual's headings have a high concentration
% ($\kappa = 5$), but the concentration across individual directional biases
% is low ($\kappa_\mu = 0.1$). The same model and priors are used.

dtLow = readData(dataDir, 'lowcorr.csv');
trueLow = readtable(fullfile(dataDir, 'lowcorr_true.csv'));
mLow = makeIndividualModel(dtLow.y, dtLow.ID, ...
    'beta0Prior', [1 0 3*pi/2], 'kappaPrior', [3 3], 'sdScale', 0.5);
fitLow = fitUnwrapModel(mLow, opts, 'model_lowcorr');
printSummary(fitLow);
plotIndividuals(fitLow, dtLow, deg2rad(trueLow.mu), 5.0*ones(height(trueLow), 1), ...
    0, 0.1, 5.0, col, 'Low correlation');
saveFig(saveFigures, figDir, 'low_correlation');
ciBeta0Low = quantile(rad2deg(unwrapCircular(getDraw(fitLow, 'b_beta0_Intercept'))), [0.025 0.5 0.975]);
fprintf('beta0 (degrees clockwise) 2.5%%, 50%%, 97.5%%: %s; 95%% CI width: %.0f°\n', ...
    mat2str(round(ciBeta0Low)), ciBeta0Low(3) - ciBeta0Low(1));
fprintf('softplus(b_kappa_Intercept) > true kappa: %.1f%%\n', ...
    100*mean(softplus(getDraw(fitLow, 'b_kappa_Intercept')) > 5.0));

%%
% As for the high-correlation example, the model overestimates $\kappa_\mu$,
% but recovers accurate estimates of $\beta_i$ and population $\kappa$.
% Perhaps unsurprisingly given the set of individual heading biases, the
% population mean direction $\beta_0$ was not successfully recovered, but the
% broad credible interval around the estimate suggests considerable
% uncertainty about the true direction. This may lead us to ask whether or not
% there is a shared population mean direction.

%% Test shared heading bias hypothesis
% We fit a model without $\kappa_\mu$: $\beta_0$ is anchored at 0° and the
% individual biases have a (near) uniform prior, $\beta_i \sim vM(0, \approx 0)$.
% A higher target acceptance ratio gives smaller sampling steps: slower, but
% more robust sampling. Without the information from $\kappa_\mu$ the
% posterior for each $\beta_i$ is periodic and unbounded, so step size tuning
% may warn that it is attempting a large number of leapfrog steps.

mNoCorr = makeIndividualModel(dtLow.y, dtLow.ID, ...
    'beta0Prior', [], 'kappaPrior', [3 3], 'sdScale', 0.5, 'betaiKappa', 1e-16);
optsSlow = opts; optsSlow.targetAccept = 0.97;
fitNoCorr = fitUnwrapModel(mNoCorr, optsSlow, 'model_nocorr');
printSummary(fitNoCorr);
looLow = looUnwrap(fitLow, opts);
looNoCorr = looUnwrap(fitNoCorr, optsSlow);
looCompare({'model_lowcorr', 'model_nocorr'}, {looLow, looNoCorr});
fprintf('ELPD difference (kappa_mu - no kappa_mu): %s\n', elpdDiff(looLow, looNoCorr));

%%
% Although this dataset was simulated with a von Mises distribution with very
% low concentration ($\kappa_\mu = 0.1$), the distribution of individual
% headings does not detectably differ from uniform. This is plausible for
% animals without a shared goal direction, and at this level of concentration
% a considerably larger number of individuals would be needed to recover any
% group-level bias. Nonetheless, both models recover the concentration of the
% average individual's headings ($\kappa \approx 5$).

%% Variable individual parameters
%
% <<../R/icons/fly.png>>
%
% _Tethered vinegar flies (Drosophila) can be induced to fly in a "virtual
% reality" LED arena, in which their wing motion causes visual cues to rotate.
% A green light spot induces orientation behaviour, with a weak population
% level bias towards the light, but highly variable performance across
% individuals._
%
% In this example, there is a moderate concentration across individual biases
% ($\kappa_\mu = 1$), and the average individual's headings are somewhat
% concentrated ($\alpha_0 = 1.5$), but the concentration within each
% individual's headings varies ($\sigma_\kappa = 0.5$), including a
% $\approx 0.1\%$ chance of completely uniform heading distributions
% ($\kappa = 0$). Despite the modest inter-individual correlation, 4
% individuals diverge from the population mean direction by more than 90°, and
% individual concentrations vary from $\kappa = 1.14$ to $\kappa = 2.66$.
%
% We use an unbiased (near uniform) von Mises prior for $\beta_0$, a prior
% expecting lower average concentration, $\alpha_0 \sim N(1, 3)$, and a wider
% prior for individual differences in concentration,
% $\sigma_\kappa \sim HalfT(3, 0, 1)$.

dtVar = readData(dataDir, 'var.csv');
trueVar = readtable(fullfile(dataDir, 'var_true.csv'));
mVar = makeIndividualModel(dtVar.y, dtVar.ID, ...
    'beta0Prior', [2 0 1e-16], 'kappaPrior', [1 3], 'sdScale', 1.0);
fitVar = fitUnwrapModel(mVar, opts, 'model_var');
printSummary(fitVar);
plotIndividuals(fitVar, dtVar, deg2rad(trueVar.mu), trueVar.kappa, ...
    0, 1.0, 1.5, col, 'Variable individuals');
saveFig(saveFigures, figDir, 'variable_individuals');
fprintf('beta0 > 0°: %.1f%% of estimates\n', ...
    100*mean(unwrapCircular(getDraw(fitVar, 'b_beta0_Intercept')) > 0));

%%
% The model accounts for individual differences in both directional bias
% ($\beta_i$) and concentration ($\alpha_i$), and the hyperparameters
% ($\kappa_\mu$ and $\sigma_\kappa$) are recovered quite accurately.
% Individuals with very low concentration and large deviations from the
% population mean are estimated to be slightly closer to the population mean
% than their true values (shrinkage). The population mean direction is
% recovered, with a wide posterior distribution for $\beta_0$ centred close
% to the true value.

%% Individual heading changes
%
% <<../R/icons/bird.png>>
%
% _Birds exhibiting migratory restlessness hop repeatedly in their intended
% flight direction, allowing this to be recorded repeatedly for the same
% individuals across multiple nights. Individuals may differ somewhat in
% their hopping directions, but are generally aligned provided they detect
% the same compass cues. A treatment that affects their perception of these
% cues causes the entire population to change their hopping direction._
%
% In this example, the treatment induces a population-level rightwards turn
% of 30°. The average individual's headings are concentrated ($\kappa = 3$)
% but vary greatly between individuals ($\sigma_\kappa = 2.5$), and heading
% biases are highly concentrated across individuals ($\kappa_\mu = 4$). There
% are large differences between individual responses: one individual turns to
% the left after treatment, and individual concentrations range from
% $\kappa = 1.4$ to $\kappa = 8.8$.
%
% The treatment $x$ predicts $\beta_0$ (|beta0 ~ 1 + x|), $\beta_i$
% (|betai ~ 0 + ID + x:ID|) and $\kappa$ (kappa ~ 0 + Intercept + x, with
% individual effects on both the intercept and x).
% Individual biases in the control follow $\beta_i \sim vM(0, \kappa_\mu)$,
% and in the treatment $\beta_i + \beta_{i,x} \sim vM(0, \kappa_{\mu,x})$,
% where $\kappa_{\mu,x} = \log(\exp(k_\mu + k_{\mu,x}) + 1)$ and
% $k_{\mu,x} \sim N(0, 1)$. Models with effects of treatment on $\kappa$ only
% ('kappa model') and with no effects of treatment ('null model') are fitted
% for comparison, all with a high target acceptance ratio.

dtHd = readData(dataDir, 'hd.csv');
trueHd = readtable(fullfile(dataDir, 'hd_true.csv'));
mHd = makeIndividualModel(dtHd.y, dtHd.ID, ...
    'x', dtHd.x, 'muTreatment', true, 'kappaTreatment', true, ...
    'beta0Prior', [1 0 3*pi/2; 1 0 pi/2], 'kappaPrior', [3 3; 0 3], ...
    'sdScale', [3 3]);
mHdKappa = makeIndividualModel(dtHd.y, dtHd.ID, ...
    'x', dtHd.x, 'muTreatment', false, 'kappaTreatment', true, ...
    'beta0Prior', [1 0 3*pi/2], 'kappaPrior', [3 3; 0 3], 'sdScale', [3 3]);
mHdNull = makeIndividualModel(dtHd.y, dtHd.ID, ...
    'beta0Prior', [1 0 3*pi/2], 'kappaPrior', [3 3], 'sdScale', 3);
fitHdTreat = fitUnwrapModel(mHd, optsSlow, 'model_hd_treatment');
fitHdKappa = fitUnwrapModel(mHdKappa, optsSlow, 'model_hd_kappa');
fitHdNull = fitUnwrapModel(mHdNull, optsSlow, 'model_hd_null');
printSummary(fitHdTreat);
printSummary(fitHdKappa);
printSummary(fitHdNull);
for v = {'b_kappa_x', 'cor_ID__kappa_Intercept__kappa_x'}
    fprintf('%s 95%% CI: %.2f to %.2f\n', v{1}, quantile(getDraw(fitHdTreat, v{1}), [0.025 0.975]));
end

%%
% Parameter recovery is quite successful for $\kappa$ and $\kappa_\mu$, given
% the minimal sample size. For |b_kappa_x| there is substantial posterior
% density for both increases and decreases in concentration in the treatment
% condition (95% CI above), and the range of estimates for |kappa_mu_x|
% largely overlaps with that for |kappa_mu| in the control condition. The
% correlation between individual differences in concentration in the control
% and in their response to treatment also has substantial posterior density
% for both correlation and anti-correlation, so there is no clear evidence of
% a relationship between variation in the two conditions (correctly, since
% none was simulated).

%% Test hypothesis of treatment effect

looHdNull = looUnwrap(fitHdNull, optsSlow);
looHdKappa = looUnwrap(fitHdKappa, optsSlow);
looHdTreat = looUnwrap(fitHdTreat, optsSlow);
lcHd = looCompare({'model_hd_null', 'model_hd_kappa', 'model_hd_treatment'}, ...
    {looHdNull, looHdKappa, looHdTreat});
fprintf('ELPD difference (mu & kappa - kappa): %s\n', elpdDiff(looHdTreat, looHdKappa));
fprintf('ELPD difference (mu & kappa - null): %s\n', elpdDiff(looHdTreat, looHdNull));

elpdPlot = [looHdNull.elpd, looHdKappa.elpd, looHdTreat.elpd];
sePlot = [looHdNull.se, looHdKappa.se, looHdTreat.se];
figure('Name', 'Model comparison', 'Color', 'w'); hold on;
cols = [col.obs; col.kappa; col.treat];
for i = 1:3
    errorbar(i, elpdPlot(i), sePlot(i), 'o', 'Color', cols(i, :), ...
        'MarkerFaceColor', cols(i, :), 'MarkerSize', 10, 'LineWidth', 2);
end
% difference between treatment and kappa models, plotted relative to the kappa model
diffTK = looHdTreat.elpd - looHdKappa.elpd;
seTK = seDiff(looHdTreat, looHdKappa);
errorbar(4, looHdKappa.elpd + diffTK, seTK, 'o', 'Color', col.rho, ...
    'MarkerFaceColor', col.rho, 'MarkerSize', 10, 'LineWidth', 2);
yline(looHdKappa.elpd, 'Color', [0.5 0.5 0.5]);
xticks(1:4); xticklabels({'null model', 'kappa model', 'mu & kappa model', 'difference'});
xlim([0.5 4.5]); ylabel('expected log predictive density');
title(sprintf('ELPD difference (mu & kappa - kappa): %.1f ± %.1f', diffTK, seTK));
saveFig(saveFigures, figDir, 'hd_model_comparison');

%%
% The results of model comparison are quite clear: the model for the
% hypothesis of an effect of treatment on both $\mu$ and $\kappa$ predicts the
% data more robustly than the other two hypotheses.

%% Plot model predictions
% Predictions are calculated separately for each condition for each
% individual: control (top row) and treatment (bottom row).

plotHeadingChange(fitHdTreat, dtHd, trueHd, col);
saveFig(saveFigures, figDir, 'hd_predictions');

%%
% We can also compare the model parameters with the population parameters:
% individual biases in the control and treatment, the average individual, and
% the change in population mean angle.

b0 = getDraw(fitHdTreat, 'b_beta0_Intercept');
b0x = getDraw(fitHdTreat, 'b_beta0_x');
kmu = getDraw(fitHdTreat, 'kappa_mu');
kmux = getDraw(fitHdTreat, 'kappa_mu_x');
k0 = softplus(getDraw(fitHdTreat, 'b_kappa_Intercept'));
figure('Name', 'Population parameters', 'Color', 'w');
tiledlayout(1, 4);
ax = nexttile; circBase(ax, 0); title(ax, 'individual biases (control)');
circPoints(ax, deg2rad(trueHd.mu), col.rho);
circArrow(ax, 0, A1(4.0), col.obs, 5);
drawCloud(ax, b0, kmu);
circArrow(ax, circMean(b0), A1(median(kmu)), col.sd, 2);
ax = nexttile; circBase(ax, 0); title(ax, 'individual biases (treatment)');
circPoints(ax, deg2rad(trueHd.mu_delta), col.rho);
circArrow(ax, deg2rad(30), A1(4.0), col.treat, 5);
drawCloud(ax, b0x, kmux);
circArrow(ax, circMean(b0x), A1(median(kmux)), col.sd, 2);
ax = nexttile; circBase(ax, 0); title(ax, 'average individual');
circArrow(ax, 0, A1(3.0), col.obs, 5);
drawCloud(ax, b0, k0);
circArrow(ax, circMean(b0), A1(median(k0)), col.sd, 2);
ax = nexttile;
vertHist(ax, rad2deg(unwrapCircular(b0x)), 'change in mean angle', [-180 180], col.sd, 0);
saveFig(saveFigures, figDir, 'hd_population');
fprintf('Rightwards turn: %.1f%% of estimates; 95%% CI: %.0f° to %.0f° (true turn 30°)\n', ...
    100*mean(wrapPi(b0x) > 0), quantile(rad2deg(unwrapCircular(b0x)), [0.025 0.975]));

%%
% The estimates of the population mean angle favour a right turn (percentage
% above), but there is a lot of uncertainty about the exact size of the turn.
% This is influenced by the wide spread of individual mean angles and turn
% angles, which the model recovered well.

%% Summary
%
% * Circular GLMs can be used to model effects on mean direction $\mu$ and
%   concentration $\kappa$.
% * The posterior density of parameter estimates can be compared with
%   expected values; for example, an expected angle (such as the 'home
%   direction') or concentration.
% * Where expected values of $\mu$ and $\kappa$ are not known, fitted
%   estimates can be compared between data for control and treatment
%   conditions as an offset parameter between the two.
% * Grouping factors, such as individual biases, can be accounted for in a
%   hierarchical model, and the properties of the population of individuals
%   estimated as hyperparameters for the concentration of individual mean
%   directions ($\kappa_\mu$) and standard deviation of individual
%   concentrations ($\sigma_\kappa$).
% * Competing hypotheses about effects within a dataset can be formalised as
%   circular GLMs with differing parameters, and their plausibility contrasted
%   directly via model-comparison methods such as expected log predictive
%   density.
