function m = makeIndividualModel(y, ID, varargin)
    % Hierarchical model: mu = beta0 + betai, kappa with individual effects
    p = struct('x', zeros(numel(y), 1), 'muTreatment', false, 'kappaTreatment', false, ...
        'beta0Prior', [1 0 3*pi/2], 'kappaPrior', [3 3], 'sdScale', 0.5, ...
        'betaiKappa', NaN); % NaN: estimate kappa_mu
    for i = 1:2:numel(varargin)
        p.(varargin{i}) = varargin{i + 1};
    end
    N = numel(y);
    x = p.x(:);
    [~, ~, g] = unique(ID);
    J = max(g);
    D = full(sparse(1:N, g, 1, N, J)); % individual indicators
    idNames = compose('b_betai_ID%d', 1:J);
    % population mean direction (beta0)
    if isempty(p.beta0Prior)
        Xmu = zeros(N, 0); muNames = {}; % beta0 ~ 0
    elseif p.muTreatment
        Xmu = [ones(N, 1) x]; muNames = {'b_beta0_Intercept', 'b_beta0_x'};
    else
        Xmu = ones(N, 1); muNames = {'b_beta0_Intercept'};
    end
    % individual biases (betai) and their von Mises priors
    if p.muTreatment
        Zb = [D D.*x];
        bNames = [idNames, compose('b_betai_ID%d:x', 1:J)];
        groups = struct('A', {[eye(J) zeros(J)], [eye(J) eye(J)]}, ...
            'kappaFixed', {NaN, NaN}, 'skmCoef', {[1 0], [1 1]}, ...
            'name', {'kappa_mu', 'kappa_mu_x'});
        skmPrior = [3 log(15) 0.6; 1 0 1.0]; % lognormal on kappa_mu, normal on change
        skmNames = {'softkappamu', 'softkappamux'};
    else
        Zb = D; bNames = idNames;
        if isnan(p.betaiKappa)
            groups = struct('A', eye(J), 'kappaFixed', NaN, 'skmCoef', 1, 'name', 'kappa_mu');
            skmPrior = [3 log(15) 0.6]; skmNames = {'softkappamu'};
        else
            groups = struct('A', eye(J), 'kappaFixed', p.betaiKappa, 'skmCoef', zeros(1, 0), ...
                'name', '');
            skmPrior = zeros(0, 3); skmNames = {};
        end
    end
    % concentration (kappa): no centring of x, so the intercept is the control
    if p.kappaTreatment
        Xk = [ones(N, 1) x]; kNames = {'b_kappa_Intercept', 'b_kappa_x'};
        Zr = [ones(N, 1) x]; reNames = {'Intercept', 'x'};
    else
        Xk = ones(N, 1); kNames = {'b_kappa_Intercept'};
        Zr = ones(N, 1); reNames = {'Intercept'};
    end
    m = makeModel(y, 'Xmu', Xmu, 'muNames', muNames, 'muPrior', p.beta0Prior, ...
        'Xk', Xk, 'kNames', kNames, 'kPrior', p.kappaPrior, ...
        'Zb', Zb, 'bNames', bNames, 'groups', groups, ...
        'skmPrior', skmPrior, 'skmNames', skmNames, ...
        'g', g, 'Zr', Zr, 'reNames', reNames, 'sdScale', p.sdScale(:));
end
