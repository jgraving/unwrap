function fit = fitUnwrapModel(m, opts, name)
    m = indexModel(m);
    fprintf('Fitting %s (%d parameters)... ', name, m.nPar);
    tic;
    chains = zeros(opts.numSamples, m.nPar, opts.numChains);
    acc = zeros(opts.numChains, 1);
    for c = 1:opts.numChains
        start = -2 + 4*rand(m.nPar, 1); % as Stan: uniform(-2, 2)
        smp = hmcSampler(@(t) unwrapLogPosterior(t, m), start, ...
            'CheckGradient', c == 1); % check analytic gradient once
        smp = tuneSampler(smp, 'StartPoint', start, ...
            'StepSizeTuningMethod', 'dual-averaging', ...
            'MassVectorTuningMethod', 'iterative-sampling', ...
            'TargetAcceptanceRatio', opts.targetAccept);
        [chains(:, :, c), ~, acc(c)] = drawSamples(smp, 'StartPoint', start, ...
            'Burnin', opts.burnin, 'NumSamples', opts.numSamples);
    end
    fprintf('%.0f s, acceptance ratio %s\n', toc, mat2str(round(acc', 2)));
    fit.name = name;
    fit.m = m;
    fit.chains = chains;
    fit.acceptance = acc;
    [S, ~, C] = size(chains);
    fit.pooled = reshape(permute(chains, [1 3 2]), S*C, m.nPar);
    [fit.draws, fit.names, fit.isCirc] = derivedDraws(fit.pooled, m);
end
