function circBase(ax, refline)
    % Circular axes: zero at the top, clockwise rotation
    hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
    t = linspace(0, 2*pi, 361);
    plot(ax, sin(t), cos(t), 'k-', 'LineWidth', 1);
    for a = (0:3)*pi/2
        plot(ax, [0.95 1.05]*sin(a), [0.95 1.05]*cos(a), 'k-');
    end
    if ~isempty(refline)
        plot(ax, [0 sin(refline)], [0 cos(refline)], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
    end
    xlim(ax, [-1.6 1.6]); ylim(ax, [-1.6 1.6]);
end
