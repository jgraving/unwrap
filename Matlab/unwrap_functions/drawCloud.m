function drawCloud(ax, mu, kappa)
    % Posterior distribution of mean vectors as a cloud of draws
    idx = randperm(numel(mu), min(1000, numel(mu)));
    rho = A1(kappa(idx));
    scatter(ax, rho.*sin(mu(idx)), rho.*cos(mu(idx)), 6, [0.85 0.2 0.1], 'filled', ...
        'MarkerFaceAlpha', 0.06, 'MarkerEdgeColor', 'none');
end
