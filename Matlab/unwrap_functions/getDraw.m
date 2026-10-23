function x = getDraw(fit, name)
    x = fit.draws(:, strcmp(fit.names, name));
end
