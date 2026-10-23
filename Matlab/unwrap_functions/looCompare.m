function tbl = looCompare(names, los)
    elpd = cellfun(@(l) l.elpd, los)';
    se = cellfun(@(l) l.se, los)';
    [~, best] = max(elpd);
    elpdDiff = elpd - elpd(best);
    seD = cellfun(@(l) seDiff(l, los{best}), los)';
    tbl = table(names(:), elpdDiff, seD, elpd, se, ...
        'VariableNames', {'Model', 'elpd_diff', 'se_diff', 'elpd', 'se_elpd'});
    tbl = sortrows(tbl, 'elpd_diff', 'descend');
    fprintf('\nLOO comparison\n');
    disp(tbl);
end
