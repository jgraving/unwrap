function plotTreatment(fit, yCtrl, yTreat, muTrue, kappaTrue, col, ttl)
    % Control and treatment predictions for the treatment models
    b0 = getDraw(fit, 'b_Intercept'); bx = getDraw(fit, 'b_x');
    k0 = getDraw(fit, 'b_kappa_Intercept'); kx = getDraw(fit, 'b_kappa_x');
    figure('Name', ttl, 'Color', 'w');
    tiledlayout(1, 3);
    ax = nexttile; circBase(ax, 0); title(ax, 'control');
    circPoints(ax, yCtrl, col.obs);
    drawCloud(ax, b0, softplus(k0));
    circArrow(ax, muTrue(1), A1(kappaTrue(1)), col.obs, 5);
    circArrow(ax, circMedian(b0), A1(softplus(median(k0))), col.sd, 2);
    ax = nexttile; circBase(ax, 0); title(ax, 'treatment');
    circPoints(ax, yTreat, col.treat);
    drawCloud(ax, b0 + bx, softplus(k0 + kx));
    circArrow(ax, muTrue(2), A1(kappaTrue(2)), col.treat, 5);
    circArrow(ax, circMedian(b0 + bx), A1(softplus(median(k0 + kx))), col.sd, 2);
end
