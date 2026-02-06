function [logp, grad] = vonMisesLogPosterior(params, data)
    % 1. Ensure params is a column vector dlarray
    params = params(:);
    if ~isa(params, 'dlarray')
        params = dlarray(params);
    end

    % 2. Define the objective function as a local variable
    % This returns ONLY the scalar log-posterior (lp)
    target = @(p) calculateLP_logic(p, data);

    % 3. dlfeval handles the gradient extraction
    [lp_dl, grad_dl] = dlfeval(target, params);
    
    % 4. Convert back to standard doubles for hmcSampler
    logp = extractdata(lp_dl);
    grad = extractdata(grad_dl);
end

% Keep this in the SAME file, at the bottom
function lp = calculateLP_logic(p, data)
    % --- Unpack ---
    K = data.Kc;
    Kk = data.Kc_kappa;
    
    b = p(1:K);
    Intercept = p(K+1);
    b_kappa = p(K+2 : K+1+Kk);
    Intercept_kappa = p(K+1+Kk+1);

    % --- Priors using normpdf ---
    lp = sum(log(normpdf(b, 0, pi/2)));
    lp = lp + sum(log(normpdf(b_kappa, 0, 0.5)));
    
    mu_Ik = log1p(exp(2.0));
    lp = lp + log(normpdf(Intercept_kappa, mu_Ik, 1.0));
    
    % Flat-ish circular prior for intercept
    lp = lp + (1e-16 * cos(Intercept));

    % --- Likelihood ---
    mu = Intercept + data.Xc * b;
    k_lin = Intercept_kappa + data.Xc_kappa * b_kappa;
    kappa = log1p(exp(k_lin)); 
    
    mu_wrapped = mod(mu + pi, 2*pi) - pi;
    
    % Von Mises density calculation
    term1 = kappa .* cos(data.Y - mu_wrapped);
    % Modified Bessel function of the first kind via besselj
    term2 = log(2 * pi * besselj(0, 1i * kappa)); 
    
    lp = lp + sum(term1 - real(term2));
end