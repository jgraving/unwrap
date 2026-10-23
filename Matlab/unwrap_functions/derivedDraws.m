function [D, names, isCirc] = derivedDraws(th, m)
    % Named draws on the scale reported by brms
    D = [th(:, m.iMu), th(:, m.iK), th(:, m.iB), th(:, m.iSkm)];
    names = [m.muNames, m.kNames, m.bNames, m.skmNames];
    isCirc = [true(1, m.Pmu), false(1, m.Pk), true(1, m.Q), false(1, m.nskm)];
    for gI = 1:numel(m.groups)
        if isnan(m.groups(gI).kappaFixed)
            D = [D, softplus(th(:, m.iSkm)*m.groups(gI).skmCoef(:))]; %#ok<AGROW>
            names = [names, m.groups(gI).name]; %#ok<AGROW>
            isCirc = [isCirc, false]; %#ok<AGROW>
        end
    end
    if m.R > 0
        sd = exp(th(:, m.iSd));
        z = th(:, m.iZ);
        z1 = z(:, 1:m.J);
        if m.R == 1
            r = sd.*z1;
            cr = zeros(size(th, 1), 0);
        else
            rho = tanh(th(:, m.iCor));
            z2 = z(:, m.J + (1:m.J));
            r = [sd(:, 1).*z1, sd(:, 2).*(rho.*z1 + sqrt(1 - rho.^2).*z2)];
            cr = rho;
        end
        rNames = {};
        for c = 1:m.R
            rNames = [rNames, compose('r_ID__kappa[%d,%s]', (1:m.J)', m.reNames{c})']; %#ok<AGROW>
        end
        sdNames = strcat('sd_ID__kappa_', m.reNames);
        corNames = {};
        if m.R == 2; corNames = {'cor_ID__kappa_Intercept__kappa_x'}; end
        D = [D, sd, cr, r];
        names = [names, sdNames, corNames, rNames];
        isCirc = [isCirc, false(1, m.R + m.nCor + m.J*m.R)];
    end
end
