function plotIndividuals(fit, dt, muTrue, kappaTrue, mu0, kappaMuTrue, kappaPopTrue, col, ttl)
    % Individual predictions and population parameters
    J = numel(muTrue);
    hasBeta0 = any(strcmp(fit.names, 'b_beta0_Intercept'));
    if hasBeta0
        b0 = getDraw(fit, 'b_beta0_Intercept');
    else
        b0 = zeros(size(fit.draws, 1), 1);
    end
    k0 = getDraw(fit, 'b_kappa_Intercept');
    figure('Name', ttl, 'Color', 'w');
    tl = tiledlayout(3, 4, 'TileSpacing', 'compact');
    title(tl, ttl);
    for i = 1:J
        ax = nexttile; circBase(ax, 0);
        circPoints(ax, dt.y(dt.ID == i), col.obs);
        mu = b0 + getDraw(fit, sprintf('b_betai_ID%d', i));
        kp = softplus(k0 + getDraw(fit, sprintf('r_ID__kappa[%d,Intercept]', i)));
        drawCloud(ax, mu, kp);
        circArrow(ax, muTrue(i), A1(kappaTrue(i)), col.obs, 4);
        circArrow(ax, circMedian(mu), A1(median(kp)), col.sd, 2);
    end
    % population of biases
    ax = nexttile; circBase(ax, 0); title(ax, 'individual biases');
    circPoints(ax, muTrue, col.rho);
    circArrow(ax, mu0, A1(kappaMuTrue), col.rho, 4);
    if any(strcmp(fit.names, 'kappa_mu'))
        kmu = getDraw(fit, 'kappa_mu');
        drawCloud(ax, b0, kmu);
        circArrow(ax, circMean(b0), A1(median(kmu)), col.sd, 2);
    end
    % population kappa
    ax = nexttile;
    vertHist(ax, softplus(k0), 'population kappa', [0 10], col.kappa, kappaPopTrue);
end
