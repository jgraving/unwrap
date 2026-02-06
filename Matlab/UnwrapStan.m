%% UnwrapStan.m
%---------------------------------------------------------------
%       USAGE: UnwrapStan
%
%      AUTHOR: James Foster                 DATE: 2026 01 16
%    
% DESCRIPTION: Runs an example circular regression using the circular
%              modulo and softplus links. 
%---------------------------------------------------------------
% Before running
% make sure you have installed:
% - https://github.com/brian-lau/MatlabStan
% - https://github.com/brian-lau/MatlabProcessManager/
% - cmdstant (recommended https://mc-stan.org/cmdstanr/)
%
% hard code the cmdstan paths and versions in:
% - stan_home.m (e.g. d = '/Users/username/.cmdstan/cmdstan-2.37.0';)
% - StanModel.m
% %%BY SETTING
%       function ver = stan_version(self)
%             ver = [2, 37, 0]%for cmdstan-2.37.0
%       end
% %%AND
%           if isempty(p.Results.stan_version)
%             self.stan_version = [2 37 0];
% %             self.stan_version = self.get_stan_version();
%          else
%             self.stan_version = p.Results.stan_version;
%          end
% %% TO CORRESPOND WITH YOUR cmdstan VERSION
%%

%% write code
unwrap_stan_code = {
    '// generated with brms 2.23.0'
    'functions {'
    '   real log_expm1(real x) {'
    '     return log(expm1(x));'
    '   }'
    '   vector log_expm1(vector x) {'
    '     return log(expm1(x));'
    '   }'
    '  real unwrap_von_mises_lpdf(real y, real mu, real kappa) {'
    '    return von_mises_lpdf(y | modulo_circular(mu), kappa);'
    '  }'
    '  real unwrap_von_mises_rng(real mu, real kappa) {'
    '    return von_mises_rng( modulo_circular(mu) , kappa);'
    '  }'
    '  real modulo_circular(real y) {'
    '    return fmod(y + pi(), 2*pi()) - pi();'
    '  }'
    '}'
    'data {'
    '  int<lower=1> N;  // total number of observations'
    '  vector[N] Y;  // response variable'
    '  int<lower=1> K;  // number of population-level effects'
    '  matrix[N, K] X;  // population-level design matrix'
    '  int<lower=1> Kc;  // number of population-level effects after centering'
    '  int<lower=1> K_kappa;  // number of population-level effects'
    '  matrix[N, K_kappa] X_kappa;  // population-level design matrix'
    '  int<lower=1> Kc_kappa;  // number of population-level effects after centering'
    '}'
    'transformed data {'
    '  matrix[N, Kc] Xc;  // centered version of X without an intercept'
    '  vector[Kc] means_X;  // column means of X before centering'
    '  matrix[N, Kc_kappa] Xc_kappa;  // centered version of X_kappa without an intercept'
    '  vector[Kc_kappa] means_X_kappa;  // column means of X_kappa before centering'
    '  for (i in 2:K) {'
    '    means_X[i - 1] = mean(X[, i]);'
    '    Xc[, i - 1] = X[, i] - means_X[i - 1];'
    '  }'
    '  for (i in 2:K_kappa) {'
    '    means_X_kappa[i - 1] = mean(X_kappa[, i]);'
    '    Xc_kappa[, i - 1] = X_kappa[, i] - means_X_kappa[i - 1];'
    '  }'
    '}'
    'parameters {'
    '  vector[Kc] b;  // regression coefficients'
    '  real Intercept;  // temporary intercept for centered predictors'
    '  vector[Kc_kappa] b_kappa;  // regression coefficients'
    '  real Intercept_kappa;  // temporary intercept for centered predictors'
    '}'
    'transformed parameters {'
    '  // prior contributions to the log posterior'
    '  real lprior = 0;'
    '  lprior += normal_lpdf(b | 0.0,pi()/2);'
    '  lprior += unwrap_von_mises_lpdf(Intercept | 0,1e-16);'
    '  lprior += normal_lpdf(b_kappa | 0.0,0.5);'
    '  lprior += normal_lpdf(Intercept_kappa | log1p_exp(2.0),1.0);'
    '}'
    'model {'
    '  // likelihood including constants'
    '    // initialize linear predictor term'
    '    vector[N] mu = rep_vector(0.0, N);'
    '    // initialize linear predictor term'
    '    vector[N] kappa = rep_vector(0.0, N);'
    '    mu += Intercept + Xc * b;'
    '    kappa += Intercept_kappa + Xc_kappa * b_kappa;'
    '    kappa = log1p_exp(kappa);'
    '    for (n in 1:N) {'
    '      target += unwrap_von_mises_lpdf(Y[n] | mu[n], kappa[n]);'
    '  }'
    '  // priors including constants'
    '  target += lprior;'
    '}'
    'generated quantities {'
    '  // actual population-level intercept'
    '  real b_Intercept = Intercept - dot_product(means_X, b);'
    '  // actual population-level intercept'
    '  real b_kappa_Intercept = Intercept_kappa - dot_product(means_X_kappa, b_kappa);'
    '}'
};

%% set up data
%enter the predictor as a matrix and response as a vector
vonmises_data = struct('X', [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;... Intercept identifier 
                           -1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1...
                           0    0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9 1.0]',... predictor variable
                     'Y', [5.65 6.15 5.86 6.13 5.99 6.12 5.96 0.12 6.13 6.27...
                           6.21 0.07 6.19 0.10 0.25 0.04 0.30 0.35 0.37 0.40 0.36],... response variable (angles)
                      'N', 21,... sample size
                      'K', 2,... number of population level effects
                       'Kc', 1); %number after centring
%typically these will be the same for kappa
vonmises_data.K_kappa = vonmises_data.K;
vonmises_data.Kc_kappa = vonmises_data.Kc;
vonmises_data.X_kappa = vonmises_data.X;
                       
%% fit the model
fit = stan('model_code',unwrap_stan_code,...
           'data',vonmises_data,...
           'stan_home', '/Users/jamesfoster/.cmdstan/cmdstan-2.37.0');

%% summarise the fit
%Rhat should be <1.01
print(fit);

% b[1]                %change in angle (slope) as function of X (in radians)
% Intercept           %angle at X = 0 (in radians)
% b_kappa[1]          %change in inverse_softplus(kappa) as a function of X
% Intercept_kappa     %inverse_softplus(kappa) at X = 0

%N.B. MatlabStan does not currently offer tools to unwrap estimates
%but the traceplot should help the user to identify periodic chains
fit.traceplot

%% Inspect estimates

%Extract samples (returns a struct)
samples = fit.extract('permuted', true); 

%Plot the distribution of the population-level effects 
figure;
%mu
subplot(4,2,1);
boxplot(samples.Intercept); % Boxplot of all 'b' coefficients
title('Posterior Distribution of Intercept');
grid on;
subplot(4,2,2);
histogram(samples.Intercept, 'Normalization', 'pdf');
title('Posterior Density of Intercept');
xlabel('Value'); ylabel('Density');
subplot(4,2,3);
boxplot(samples.b); % Boxplot of all 'b' coefficients
title('Posterior Distribution of b');
grid on;
subplot(4,2,4);
histogram(samples.b, 'Normalization', 'pdf');
title('Posterior Density of b');
xlabel('Value'); ylabel('Density');
%kappa
subplot(4,2,5);
boxplot(samples.Intercept_kappa); % Boxplot of all 'b' coefficients
title('Posterior Distribution of kappa Intercept');
grid on;
subplot(4,2,6);
histogram(samples.Intercept_kappa, 'Normalization', 'pdf');
title('Posterior Density of kappa Intercept');
xlabel('Value'); ylabel('Density');
subplot(4,2,7);
boxplot(samples.b_kappa); % Boxplot of all 'b' coefficients
title('Posterior Distribution of kappa b');
grid on;
subplot(4,2,8);
histogram(samples.b_kappa, 'Normalization', 'pdf');
title('Posterior Density of kappa b');
xlabel('Value'); ylabel('Density');