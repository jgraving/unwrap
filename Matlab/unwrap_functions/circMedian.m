function m = circMedian(x)
    % angle minimising the mean circular distance to all points
    x = wrapPi(x(:));
    cand = x(round(linspace(1, numel(x), min(numel(x), 1000))));
    dev = arrayfun(@(c) mean(pi - abs(pi - abs(x - c))), cand);
    [~, i] = min(dev);
    m = cand(i);
end
