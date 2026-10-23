function tbl = printSummary(fit, show)
    % Posterior summary; circular parameters are summarised after unwrapping
    % (Rhat_raw shows the apparent lack of convergence before unwrapping)
    if nargin < 2; show = true; end
    D = fit.draws;
    C = size(fit.chains, 3);
    S = size(D, 1)/C;
    nP = size(D, 2);
    rhatRaw = NaN(nP, 1); rhat = NaN(nP, 1); ess = NaN(nP, 1);
    Du = D;
    for p = 1:nP
        if fit.isCirc(p)
            Du(:, p) = unwrapCircular(D(:, p));
        end
        rhatRaw(p) = rhatClassic(reshape(D(:, p), S, C));
        rhat(p) = rhatClassic(reshape(Du(:, p), S, C));
        ess(p) = estimateESS(reshape(Du(:, p), S, 1, C));
    end
    q = quantile(Du, [0.025 0.975])';
    tbl = table(fit.names(:), mean(Du)', std(Du)', q(:, 1), q(:, 2), rhat, rhatRaw, round(ess), ...
        'VariableNames', {'Parameter', 'Mean', 'SD', 'Q2_5', 'Q97_5', 'Rhat', 'Rhat_raw', 'ESS'});
    if show
        fprintf('\n%s (circular parameters unwrapped, radians)\n', fit.name);
        disp(tbl);
    end
end
