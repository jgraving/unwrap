function circPoints(ax, ang, colr)
    % Stacked points in 5° bins
    binW = 2*pi/72;
    b = mod(round(mod(ang, 2*pi)/binW), 72);
    ub = unique(b);
    for i = 1:numel(ub)
        n = sum(b == ub(i));
        r = 1 + 0.06*(1:n);
        a = ub(i)*binW;
        scatter(ax, r*sin(a), r*cos(a), 15, colr, 'filled');
    end
end
