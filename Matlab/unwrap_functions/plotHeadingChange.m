function plotHeadingChange(fit, dt, trueHd, col)
    % Individual heading change predictions (control & treatment)
    J = height(trueHd);
    b0 = getDraw(fit, 'b_beta0_Intercept'); b0x = getDraw(fit, 'b_beta0_x');
    k0 = getDraw(fit, 'b_kappa_Intercept'); k0x = getDraw(fit, 'b_kappa_x');
    figure('Name', 'Individual heading change predictions', 'Color', 'w');
    tiledlayout(2, J, 'TileSpacing', 'compact');
    for i = 1:J
        ax = nexttile; circBase(ax, 0);
        circPoints(ax, dt.y(dt.ID == i & dt.x == 0), col.obs);
        mu = b0 + getDraw(fit, sprintf('b_betai_ID%d', i));
        kp = softplus(k0 + getDraw(fit, sprintf('r_ID__kappa[%d,Intercept]', i)));
        drawCloud(ax, mu, kp);
        circArrow(ax, deg2rad(trueHd.mu(i)), A1(trueHd.kappa(i)), col.obs, 4);
        circArrow(ax, circMedian(mu), A1(median(kp)), col.sd, 2);
    end
    for i = 1:J
        ax = nexttile; circBase(ax, 0);
        circPoints(ax, dt.y(dt.ID == i & dt.x == 1), col.treat);
        mu = b0 + b0x + getDraw(fit, sprintf('b_betai_ID%d', i)) + ...
            getDraw(fit, sprintf('b_betai_ID%d:x', i));
        kp = softplus(k0 + k0x + getDraw(fit, sprintf('r_ID__kappa[%d,Intercept]', i)) + ...
            getDraw(fit, sprintf('r_ID__kappa[%d,x]', i)));
        drawCloud(ax, mu, kp);
        circArrow(ax, deg2rad(trueHd.mu_delta(i)), A1(trueHd.kappa(i)), col.treat, 4);
        circArrow(ax, circMedian(mu), A1(median(kp)), col.sd, 2);
    end
end
