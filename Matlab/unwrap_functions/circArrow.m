function circArrow(ax, mu, rho, colr, lw)
    quiver(ax, 0, 0, rho*sin(mu), rho*cos(mu), 0, 'Color', colr, 'LineWidth', lw, ...
        'MaxHeadSize', 0.4/max(rho, 0.1));
end
