function [r, sd, rho, cc, z] = randomEffects(theta, m)
    % Transform standardised individual effects into effects on kappa
    z = reshape(theta(m.iZ), m.J, m.R);
    sd = exp(theta(m.iSd)); sd = sd(:);
    rho = 0; cc = 1;
    if m.R == 1
        r = sd*z;
    else
        rho = tanh(theta(m.iCor)); cc = sqrt(1 - rho^2);
        r = [sd(1)*z(:, 1), sd(2)*(rho*z(:, 1) + cc*z(:, 2))];
    end
end
