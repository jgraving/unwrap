function q = gpdQuantile(p, k, sigma)
    q = sigma*expm1(-k*log1p(-p))/k;
end
