function y = logSumExp(x)
    mx = max(x);
    y = mx + log(sum(exp(x - mx)));
end
