function txt = elpdDiff(l1, l2)
    % ELPD difference between two models (first minus second) and its
    % standard error, as text "difference ± SE" (as in loo_compare)
    txt = sprintf('%.1f ± %.1f', sum(l1.elpd_i - l2.elpd_i), seDiff(l1, l2));
end
