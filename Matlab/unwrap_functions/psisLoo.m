function [elpd_i, k] = psisLoo(LL)
    % Pareto-smoothed importance sampling LOO (Vehtari et al. 2017, as in the loo package)
    [S, N] = size(LL);
    M = ceil(min(0.2*S, 3*sqrt(S))); % tail length
    elpd_i = zeros(N, 1);
    k = zeros(N, 1);
    for i = 1:N
        lw = -LL(:, i);
        lw = lw - max(lw);
        [lwS, ord] = sort(lw);
        cutoff = lwS(S - M);
        [k(i), sigma] = gpdFit(exp(lwS(S - M + 1:S)) - exp(cutoff));
        if isfinite(k(i))
            p = ((1:M)' - 0.5)/M;
            lwS(S - M + 1:S) = log(gpdQuantile(p, k(i), sigma) + exp(cutoff));
            lw(ord) = lwS;
        end
        lw(lw > 0) = 0; % truncate at the largest raw weight
        lw = lw - logSumExp(lw);
        elpd_i(i) = logSumExp(lw + LL(:, i));
    end
end
