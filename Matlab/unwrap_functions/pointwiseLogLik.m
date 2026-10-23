function LL = pointwiseLogLik(th, m)
    LL = zeros(size(th, 1), m.N);
    for s = 1:size(th, 1)
        LL(s, :) = unwrapLogLik(th(s, :), m)';
    end
end
