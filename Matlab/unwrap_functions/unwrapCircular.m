function u = unwrapCircular(x)
    % Remove discontinuities by centring on the mean direction
    mux = circMean(x);
    u = wrapPi(x - mux) + mux;
end
