function s = seDiff(l1, l2)
    d = l1.elpd_i - l2.elpd_i;
    s = sqrt(numel(d)*var(d));
end
