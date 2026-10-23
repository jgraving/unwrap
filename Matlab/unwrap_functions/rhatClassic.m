function rhat = rhatClassic(X)
    % X: numSamples x numChains (as in MatlabOnly.m)
    n = size(X, 1);
    W = mean(var(X, 0, 1));
    B = n*var(mean(X, 1));
    varHat = ((n - 1)/n)*W + B/n;
    rhat = sqrt(varHat/W);
    if W == 0; rhat = 1; end
end
