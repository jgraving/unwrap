%% set up data
%enter the predictor as a matrix and response as a vector
X = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;... Intercept identifier 
    -1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1...
    0    0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9 1.0]';% predictor variable
Y = [5.65 6.15 5.86 6.13 5.99 6.12 5.96 0.12 6.13 6.27...
    6.21 0.07 6.19 0.10 0.25 0.04 0.30 0.35 0.37 0.40 0.36];% response variable (angles)
N = 21;% sample size
K = 2;% number of population level effects
%Kc = 1;% %number after centring
%typically these will be the same for kappa
K_kappa = K;
%Kc_kappa = Kc;
X_kappa = X;

%% --- Preprocessing (Transformed Data) ---
Kc = size(X, 2) - 1;
means_X = mean(X(:, 2:end));
Xc = X(:, 2:end) - means_X;

Kc_kappa = size(X_kappa, 2) - 1;
means_X_kappa = mean(X_kappa(:, 2:end));
Xc_kappa = X_kappa(:, 2:end) - means_X_kappa;

data_struct = struct('N', size(Y,1), 'Y', Y, 'Xc', Xc, 'Kc', Kc, ...
                     'Xc_kappa', Xc_kappa, 'Kc_kappa', Kc_kappa);

%% --- Initialize Sampler ---
numParams = Kc + 1 + Kc_kappa + 1;
startPoint = zeros(numParams, 1);

% vonMisesLogPosterior samples from lp_internal
% lp_internal defines the model structure and priors
% Create the sampler object
sampler = hmcSampler(@(p) vonMisesLogPosterior(p, data_struct), startPoint);

%% --- Tuning and Sampling ---
[sampler, rolledOut] = estimateStepsize(sampler, 'NumStepSizeIterations', 10);
[sampler, rolledOut] = estimateMassMatrix(sampler, 'NumSamples', 500);

%% Draw samples (The Matlab equivalent to Stan's 'sampling' or 'draws')
chains = 4;
samplesPerChain = 1000;
[chains, acceptanceRatio] = drawSamples(sampler, 'NumSamples', samplesPerChain, 'NumChains', chains);


