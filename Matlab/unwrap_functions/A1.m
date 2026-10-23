function y = A1(k)
    % Mean vector length of a von Mises distribution
    y = besseli(1, k, 1)./besseli(0, k, 1);
end
