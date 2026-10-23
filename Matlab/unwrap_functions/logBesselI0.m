function y = logBesselI0(k)
    y = log(besseli(0, k, 1)) + abs(k);
end
