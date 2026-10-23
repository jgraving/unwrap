function [V, p] = vTest(x, mu0)
    % v-test (Rayleigh test with a specified mean direction), normal approximation
    n = numel(x);
    R = abs(mean(exp(1i*x)));
    V = R*cos(circMean(x) - mu0);
    u = V*sqrt(2*n);
    p = 1 - normcdf(u);
end
