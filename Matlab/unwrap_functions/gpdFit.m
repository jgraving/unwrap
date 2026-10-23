function [k, sigma] = gpdFit(x)
    % Generalised Pareto fit (Zhang & Stephens 2009, with the loo package's weak prior)
    x = sort(x(:));
    N = numel(x);
    prior = 3;
    M = 30 + floor(sqrt(N));
    jj = (1:M)';
    xstar = x(floor(N/4 + 0.5));
    theta = 1/x(N) + (1 - sqrt(M./(jj - 0.5)))/prior/xstar;
    lTheta = zeros(M, 1);
    for j = 1:M
        a = -theta(j);
        kj = mean(log1p(a*x));
        lTheta(j) = N*(log(a/kj) - kj - 1);
    end
    wTheta = exp(lTheta - logSumExp(lTheta));
    thetaHat = sum(theta.*wTheta);
    k = mean(log1p(-thetaHat*x));
    sigma = -k/thetaHat;
    k = (N*k + 10*0.5)/(N + 10); % shrink towards 0.5
    if ~isfinite(k); k = Inf; end
end
