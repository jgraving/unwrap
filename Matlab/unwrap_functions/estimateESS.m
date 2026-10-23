function ess = estimateESS(x)
    % Effective sample size (as in MatlabOnly.m)
    [numSamples, ~, numChains] = size(x);
    lagMax = min(numSamples - 1, 1000);
    meanRho = zeros(lagMax, 1);
    for c = 1:numChains
        xc = x(:, 1, c);
        xc = xc - mean(xc);
        denom = sum(xc.^2);
        if denom == 0
            ess = numSamples*numChains;
            return
        end
        rho = zeros(lagMax, 1);
        for lag = 1:lagMax
            rho(lag) = sum(xc(1:end - lag).*xc(1 + lag:end))/denom;
        end
        meanRho = meanRho + rho;
    end
    meanRho = meanRho/numChains;
    tau = 1;
    for k = 1:2:(lagMax - 1)
        pairSum = meanRho(k) + meanRho(k + 1);
        if pairSum <= 0
            break
        end
        tau = tau + 2*pairSum;
    end
    ess = numSamples*numChains/tau;
end
