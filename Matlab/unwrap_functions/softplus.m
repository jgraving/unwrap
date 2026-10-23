function y = softplus(x)
    y = max(x, 0) + log1p(exp(-abs(x)));
end
