function lo = looUnwrap(fit, opts)
    m = fit.m;
    LL = pointwiseLogLik(fit.pooled, m);
    [lo.elpd_i, lo.k] = psisLoo(LL);
    kThresh = 0.7;
    bad = find(lo.k > kThresh);
    if ~isempty(bad)
        if opts.reloo && numel(bad) <= opts.maxRefits
            fprintf('%s: %d observation(s) with Pareto k > %.1f, refitting without them\n', ...
                fit.name, numel(bad), kThresh);
            for i = bad(:)'
                keep = true(m.N, 1); keep(i) = false;
                fitI = fitUnwrapModel(subsetModel(m, keep), opts, ...
                    sprintf('%s without obs %d', fit.name, i));
                llI = pointwiseLogLik(fitI.pooled, m);
                lo.elpd_i(i) = logSumExp(llI(:, i)) - log(size(llI, 1));
            end
        else
            warning('%s: %d observation(s) with Pareto k > %.1f', fit.name, numel(bad), kThresh);
        end
    end
    lo.elpd = sum(lo.elpd_i);
    lo.se = sqrt(m.N*var(lo.elpd_i));
end
