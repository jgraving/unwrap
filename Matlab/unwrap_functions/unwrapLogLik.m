function ll = unwrapLogLik(theta, m)
    % Pointwise log likelihood (for LOO)
    theta = theta(:);
    etaMu = m.Xmu*theta(m.iMu) + m.Zb*theta(m.iB);
    etaK = m.Xk*theta(m.iK);
    if m.R > 0
        r = randomEffects(theta, m);
        etaK = etaK + sum(m.Zr .* r(m.g, :), 2);
    end
    kappa = softplus(etaK);
    ll = kappa.*cos(m.Y - etaMu) - log(2*pi) - logBesselI0(kappa);
end
