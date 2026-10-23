function [lp, g] = priorLogPdf(x, pr)
    % Prior log density: [1 mean sd] normal, [2 mu kappa] von Mises
    if pr(1) == 1
        lp = -0.5*((x - pr(2))/pr(3))^2 - log(pr(3)*sqrt(2*pi));
        g = -(x - pr(2))/pr(3)^2;
    else
        lp = pr(3)*cos(x - pr(2)) - log(2*pi) - logBesselI0(pr(3));
        g = -pr(3)*sin(x - pr(2));
    end
end
