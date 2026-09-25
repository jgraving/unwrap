%% set up data
% Define the predictor matrix and response vector for the circular regression model.
% The first column acts as an intercept indicator and the second column contains
% the predictor values used in the linear mean and concentration models.
X = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;... Intercept identifier
    -1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1...
     0    0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9 1.0]'; % predictor variable

% Store the observed circular response values as a column vector.
Y = [5.65 6.15 5.86 6.13 5.99 6.12 5.96 0.12 6.13 6.27...
     6.21 0.07 6.19 0.10 0.25 0.04 0.30 0.35 0.37 0.40 0.36]';

% Define the sample size and the number of population-level effects.
N = length(Y); % sample size
K = 2; % number of population level effects

% Reuse the same design matrix structure for the concentration model.
K_kappa = K;
X_kappa = X;

%% --- Preprocessing (Transformed Data) ---
% Remove the intercept column and center the remaining predictor values.
Kc = size(X, 2) - 1;
means_X = mean(X(:, 2:end));
Xc = X(:, 2:end) - means_X;

% Repeat the same preprocessing for the concentration predictor matrix.
Kc_kappa = size(X_kappa, 2) - 1;
means_X_kappa = mean(X_kappa(:, 2:end));
Xc_kappa = X_kappa(:, 2:end) - means_X_kappa;

% Package all model inputs into a struct passed to the log-posterior function.
data_struct = struct('N', N, 'Y', Y, 'Xc', Xc, 'Kc', Kc, ...
                     'Xc_kappa', Xc_kappa, 'Kc_kappa', Kc_kappa);

%% --- Initialize Sampler ---
% Set the parameter vector size and a zero-valued initial guess.
numParams = Kc + 1 + Kc_kappa + 1;
startPoint = zeros(numParams, 1);

% Configure multiple chains for basic MCMC diagnostics.
numChains = 4;
samplesPerChain = 2000;
samples = cell(numChains, 1);
acceptanceRatio = zeros(numChains, 1);

% Run one tuned HMC chain at a time.
for c = 1:numChains
    chainStart = startPoint + 0.1 * randn(numParams, 1);

    sampler = hmcSampler(@(p) vonMisesLogPosterior(p, data_struct), chainStart);

    % Tune the sampler before drawing samples.
    sampler = tuneSampler(sampler, ...
        'StepSizeTuningMethod', 'dual-averaging', ...
        'MassVectorTuningMethod', 'iterative-sampling');

    % Draw samples from the tuned sampler.
    [chain, ~, accratio] = drawSamples(sampler, ...
        'NumSamples', samplesPerChain);

    % Store the chain and its acceptance ratio.
    samples{c} = chain;
    acceptanceRatio(c) = accratio;
end

% Combine all chains into a single 3-D array for post-processing.
sampleArray = cat(3, samples{:});

%% --- MCMC Diagnostics & Plots ---
% Assign readable labels for each model parameter.
[paramMean, paramSD, rhat, ess, hdiLower, hdiUpper] = mcmcSummaries(sampleArray);

% Assign readable labels for each model parameter.
paramNames = [compose("b%d", 1:Kc), "Intercept", ...
    compose("b_kappa%d", 1:Kc_kappa), "Intercept_kappa"];

summaryTbl = table(paramNames(:), paramMean, paramSD, rhat, ess, hdiLower, hdiUpper, ...
    'VariableNames', {'Parameter', 'Mean', 'StdDev', 'Rhat', 'ESS', 'HDI_Lower', 'HDI_Upper'});
disp(summaryTbl);

chainTbl = table((1:numChains)', acceptanceRatio(:), ...
    'VariableNames', {'Chain', 'AcceptanceRatio'});
disp(chainTbl);

% Plot trace plots for each parameter across chains.
figure('Name', 'MCMC Trace Plots', 'Color', 'w');
for i = 1:numParams
    subplot(numParams, 1, i);
    hold on;
    for c = 1:numChains
        plot(sampleArray(:, i, c), 'DisplayName', sprintf('Chain %d', c));
    end
    ylabel(paramNames{i});
    grid on;
    if i == 1
        title('Trace Plots');
        legend('Location', 'eastoutside');
    end
    if i == numParams
        xlabel('Iteration');
    end
end




%% --- Plot Predictions ---
[numSamples, ~, numChains] = size(sampleArray);
% Stack chains (dim 3) under each other before flattening, keeping parameters in columns.
pooledSamples = reshape(permute(sampleArray, [1 3 2]), numSamples * numChains, numParams);

xObs = X(:, 2);
xGrid = linspace(min(xObs), max(xObs), 200)';
xGridC = xGrid - means_X;

muDraws = zeros(numel(xGrid), size(pooledSamples, 1));

for s = 1:size(pooledSamples, 1)
    p = pooledSamples(s, :)';
    b = p(1:Kc);
    intercept = p(Kc + 1);
    muDraws(:, s) = intercept + xGridC * b;
end

muLinearMean = mean(muDraws, 2);
muCircularMean = angle(mean(exp(1i * muDraws), 2));
yWrapped = mod(Y + pi, 2 * pi) - pi;

% Wrap each draw's line to [-pi, pi) and break it where it wraps,
% then join all draws into one NaN-separated line for fast plotting.
drawsWrapped = mod(muDraws + pi, 2 * pi) - pi;
drawsWrapped([false(1, size(drawsWrapped, 2)); abs(diff(drawsWrapped)) > pi]) = NaN;
xDraws = repmat([xGrid; NaN], 1, size(drawsWrapped, 2));
yDraws = [drawsWrapped; NaN(1, size(drawsWrapped, 2))];

figure('Name', 'Model Predictions vs Data', 'Color', 'w');
hold on;
hDraws = plot(xDraws(:), yDraws(:), '-', 'Color', [0.5 0.5 0.5 0.02], 'LineWidth', 0.5);
hData = scatter(xObs, yWrapped, 45, 'filled', 'MarkerFaceAlpha', 0.7);
hLinear = plot(xGrid, muLinearMean, 'k--', 'LineWidth', 1.5);
hCircular = plot(xGrid, muCircularMean, 'r-', 'LineWidth', 2);
grid on;
xlabel('Predictor');
ylabel('Circular Response');
title('Observed Data and Posterior Predictions');
legend([hData, hDraws, hLinear, hCircular], ...
    {'Observed data', 'Posterior draws', 'Latent mean', 'Circular mean'}, 'Location', 'best');


%% --- Helper Functions ---
% Evaluate the log posterior density for the sampler.
% A finite-difference gradient is returned when requested.
function [logp, grad] = vonMisesLogPosterior(params, data)
    params = params(:);
    logp = calculateLP_logic(params, data);

    if nargout > 1
        grad = finiteDiffGrad(@(p) calculateLP_logic(p, data), params);
    end
end

% Compute the log posterior from the priors and von Mises likelihood.
function lp = calculateLP_logic(p, data)
    K = data.Kc;
    Kk = data.Kc_kappa;

    % Unpack the regression coefficients and intercept terms.
    b = p(1:K);
    Intercept = p(K + 1);
    b_kappa = p(K + 2 : K + 1 + Kk);
    Intercept_kappa = p(K + 1 + Kk + 1);

    % Gaussian priors on the slope parameters.
    lp = sum(-0.5 * (b ./ (pi / 2)).^2 - log((pi / 2) * sqrt(2 * pi)));
    lp = lp + sum(-0.5 * (b_kappa ./ 0.5).^2 - log(0.5 * sqrt(2 * pi)));

    % Prior on the concentration intercept.
    mu_Ik = log1p(exp(2.0));
    lp = lp + (-0.5 * (Intercept_kappa - mu_Ik).^2 - log(sqrt(2 * pi)));

    % Weak circular prior term for the location intercept.
    lp = lp + 1e-16 * cos(Intercept);

    % Linear predictors for circular mean and concentration.
    mu = Intercept + data.Xc * b;
    k_lin = Intercept_kappa + data.Xc_kappa * b_kappa;
    kappa = log1p(exp(k_lin));

    % Circular likelihood contribution.
    mu_wrapped = mod(mu + pi, 2 * pi) - pi;

    term1 = kappa .* cos(data.Y - mu_wrapped);
    log_i0 = log_besseli0_stable(kappa);
    term2 = log(2 * pi) + log_i0;

    lp = lp + sum(term1 - term2);
end

% Stable evaluation of log(I0(x)) for the von Mises normalizing constant.
function log_i0 = log_besseli0_stable(x)
    log_i0 = log(besseli(0, x, 1)) + abs(x);
end

% Central finite-difference gradient used when an analytic gradient is not provided.
function grad = finiteDiffGrad(fun, x)
    h = 1e-6;
    grad = zeros(size(x));

    for i = 1:numel(x)
        xp = x;
        xm = x;
        xp(i) = xp(i) + h;
        xm(i) = xm(i) - h;
        grad(i) = (fun(xp) - fun(xm)) / (2 * h);
    end
end

function [paramMean, paramSD, rhat, ess, hdiLower, hdiUpper] = mcmcSummaries(sampleArray)
    [numSamples, numParams, ~] = size(sampleArray);

    pooled = reshape(permute(sampleArray, [1 3 2]), [], numParams);
    paramMean = mean(pooled, 1)';
    paramSD = std(pooled, 0, 1)';

    chainMeans = squeeze(mean(sampleArray, 1));
    chainVars = squeeze(var(sampleArray, 0, 1));

    W = mean(chainVars, 2);
    B = numSamples * var(chainMeans, 0, 2);
    varHat = ((numSamples - 1) / numSamples) * W + (B / numSamples);

    rhat = sqrt(varHat ./ W);
    rhat(W == 0) = 1;

    ess = zeros(numParams, 1);
    hdiLower = zeros(numParams, 1);
    hdiUpper = zeros(numParams, 1);

    for p = 1:numParams
        x = pooled(:, p);
        ess(p) = estimateESS(sampleArray(:, p, :));
        [hdiLower(p), hdiUpper(p)] = highestDensityInterval(x, 0.95);
    end
end

function ess = estimateESS(x)
[numSamples, ~, numChains] = size(x);
lagMax = min(numSamples - 1, 1000);

meanRho = zeros(lagMax, 1);

for c = 1:numChains
    xc = x(:, 1, c);
    xc = xc - mean(xc);

    denom = sum(xc.^2);
    if denom == 0
        ess = numSamples * numChains;
        return
    end

    rho = zeros(lagMax, 1);
    for lag = 1:lagMax
        rho(lag) = sum(xc(1:end-lag) .* xc(1+lag:end)) / denom;
    end

    meanRho = meanRho + rho;
end

meanRho = meanRho / numChains;

tau = 1;
for k = 1:2:(lagMax - 1)
    pairSum = meanRho(k) + meanRho(k + 1);
    if pairSum <= 0
        break
    end
    tau = tau + 2 * pairSum;
end

ess = numSamples * numChains / tau;
end

function [hdiLow, hdiHigh] = highestDensityInterval(x, credMass)
x = sort(x(:));
n = numel(x);

if n == 0
    hdiLow = NaN;
    hdiHigh = NaN;
    return
end

m = max(1, floor(credMass * n));
if m >= n
    hdiLow = x(1);
    hdiHigh = x(end);
    return
end

widths = x(m + 1:end) - x(1:end - m);
[~, idx] = min(widths);

hdiLow = x(idx);
hdiHigh = x(idx + m);
end