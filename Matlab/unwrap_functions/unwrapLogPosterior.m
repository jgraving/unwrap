function [lp, grad] = unwrapLogPosterior(theta, m)
    theta = theta(:);
    grad = zeros(m.nPar, 1);
    bmu = theta(m.iMu); bk = theta(m.iK); bi = theta(m.iB); skm = theta(m.iSkm);

    % Linear predictors
    etaMu = m.Xmu*bmu + m.Zb*bi;
    etaK = m.Xk*bk;
    if m.R > 0
        [r, sd, rho, cc, z] = randomEffects(theta, m);
        etaK = etaK + sum(m.Zr .* r(m.g, :), 2);
    end
    kappa = softplus(etaK);

    % von Mises likelihood with the circular modulo link:
    % cos(y - mod(mu)) = cos(y - mu), so the derivatives are unaffected
    d = m.Y - etaMu;
    lp = sum(kappa.*cos(d) - log(2*pi) - logBesselI0(kappa));
    dMu = kappa.*sin(d);
    dEtaK = (cos(d) - A1(kappa)).*logistic(etaK);
    grad(m.iMu) = m.Xmu'*dMu;
    grad(m.iB) = m.Zb'*dMu;
    grad(m.iK) = m.Xk'*dEtaK;

    % Priors on population-level effects
    for p = 1:m.Pmu
        [lpp, gp] = priorLogPdf(bmu(p), m.muPrior(p, :));
        lp = lp + lpp; grad(m.iMu(p)) = grad(m.iMu(p)) + gp;
    end
    for p = 1:m.Pk
        [lpp, gp] = priorLogPdf(bk(p), [1 m.kPrior(p, :)]);
        lp = lp + lpp; grad(m.iK(p)) = grad(m.iK(p)) + gp;
    end

    % von Mises priors on individual biases
    for gI = 1:numel(m.groups)
        grp = m.groups(gI);
        v = grp.A*bi;
        nJ = numel(v);
        if isnan(grp.kappaFixed)
            s = grp.skmCoef*skm;
            kv = softplus(s);
        else
            kv = grp.kappaFixed;
        end
        lp = lp + sum(kv*cos(v)) - nJ*(log(2*pi) + logBesselI0(kv));
        grad(m.iB) = grad(m.iB) - grp.A'*(kv*sin(v));
        if isnan(grp.kappaFixed)
            dkv = sum(cos(v)) - nJ*A1(kv);
            grad(m.iSkm) = grad(m.iSkm) + grp.skmCoef(:)*dkv*logistic(s);
        end
    end
    % Priors on the (inverse-softplus) concentration of individual biases
    for j = 1:m.nskm
        pr = m.skmPrior(j, :);
        if pr(1) == 3 % lognormal on softplus(skm), as in brms: no Jacobian
            v = softplus(skm(j));
            lp = lp - log(v) - log(pr(3)*sqrt(2*pi)) - (log(v) - pr(2))^2/(2*pr(3)^2);
            grad(m.iSkm(j)) = grad(m.iSkm(j)) + ...
                (-1/v - (log(v) - pr(2))/(pr(3)^2*v))*logistic(skm(j));
        else
            [lpp, gp] = priorLogPdf(skm(j), pr);
            lp = lp + lpp; grad(m.iSkm(j)) = grad(m.iSkm(j)) + gp;
        end
    end

    % Individual effects on kappa
    if m.R > 0
        G = zeros(m.J, m.R); % gradient with respect to r_ID
        for c = 1:m.R
            G(:, c) = accumarray(m.g, dEtaK.*m.Zr(:, c), [m.J 1]);
        end
        if m.R == 1
            gZ = sd*G;
            gSd = sum(G.*z);
        else
            gZ = [sd(1)*G(:, 1) + rho*sd(2)*G(:, 2), cc*sd(2)*G(:, 2)];
            gSd = [sum(G(:, 1).*z(:, 1)); sum(G(:, 2).*(rho*z(:, 1) + cc*z(:, 2)))];
            gRho = sum(sd(2)*G(:, 2).*(z(:, 1) - rho/cc*z(:, 2)));
            % LKJ(1) prior: uniform correlation, sampled as atanh(rho) (+ Jacobian)
            lp = lp + log(1 - rho^2);
            grad(m.iCor) = gRho*(1 - rho^2) - 2*rho;
        end
        % standard normal z
        lp = lp - 0.5*sum(z(:).^2) - 0.5*numel(z)*log(2*pi);
        grad(m.iZ) = gZ(:) - z(:);
        % half Student-t(3, 0, s) on sd, sampled as log(sd) (+ Jacobian)
        s = m.sdScale(:);
        lp = lp + sum(log(2) + gammaln(2) - gammaln(1.5) - 0.5*log(3*pi) - log(s) ...
            - 2*log(1 + (sd./s).^2/3) + log(sd));
        grad(m.iSd) = gSd(:).*sd - 4*sd.^2./(3*s.^2 + sd.^2) + 1;
    end
end
