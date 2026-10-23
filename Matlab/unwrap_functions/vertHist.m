function vertHist(ax, x, ttl, lims, colr, refval)
    histogram(ax, x, 'Orientation', 'horizontal', 'FaceColor', colr, 'FaceAlpha', 0.4, ...
        'EdgeColor', 'none', 'BinLimits', lims);
    ylim(ax, lims);
    yline(ax, refval, 'Color', [0.6 0.6 0.6], 'LineWidth', 4);
    title(ax, ttl);
end
