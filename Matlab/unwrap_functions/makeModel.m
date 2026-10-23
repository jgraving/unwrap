function m = makeModel(Y, varargin)
    % Generic unwrap model:
    %  mu    = Xmu*b_mu + Zb*b_i                     (circular modulo link)
    %  kappa = softplus(Xk*b_kappa + sum(Zr .* r_ID)) (softplus link)
    %  y ~ vonMises(mu, kappa)
    % Individual biases b_i have von Mises priors, defined in groups as
    %  A*b_i ~ vonMises(0, softplus(skmCoef*softkappamu)) (or a fixed kappa)
    % Individual effects on kappa (r_ID) use a non-centred parameterisation:
    %  r_ID[j,:] = diag(sd) * L * z_j, z_j ~ N(0,1), L = chol(correlation)
    N = numel(Y);
    m = struct('Y', Y(:), ...
        'Xmu', zeros(N, 0), 'muNames', {{}}, 'muPrior', zeros(0, 3), ...
        'Xk', zeros(N, 0), 'kNames', {{}}, 'kPrior', zeros(0, 2), ...
        'Zb', zeros(N, 0), 'bNames', {{}}, ...
        'groups', struct('A', {}, 'kappaFixed', {}, 'skmCoef', {}, 'name', {}), ...
        'skmPrior', zeros(0, 3), 'skmNames', {{}}, ...
        'g', ones(N, 1), 'Zr', zeros(N, 0), 'reNames', {{}}, 'sdScale', zeros(0, 1));
    for i = 1:2:numel(varargin)
        m.(varargin{i}) = varargin{i + 1};
    end
end
