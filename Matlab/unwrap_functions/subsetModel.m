function m = subsetModel(m, keep)
    % Keep a subset of observations (for exact leave-one-out refits)
    m.Y = m.Y(keep); m.Xmu = m.Xmu(keep, :); m.Xk = m.Xk(keep, :);
    m.Zb = m.Zb(keep, :); m.Zr = m.Zr(keep, :); m.g = m.g(keep);
end
