function m = indexModel(m)
    % Parameter indices
    m.N = numel(m.Y);
    m.Pmu = size(m.Xmu, 2); m.Pk = size(m.Xk, 2); m.Q = size(m.Zb, 2);
    m.nskm = size(m.skmPrior, 1); m.R = size(m.Zr, 2);
    m.J = (m.R > 0) * max(m.g);
    m.nCor = double(m.R == 2);
    off = 0;
    m.iMu = off + (1:m.Pmu)'; off = off + m.Pmu;
    m.iK = off + (1:m.Pk)'; off = off + m.Pk;
    m.iB = off + (1:m.Q)'; off = off + m.Q;
    m.iSkm = off + (1:m.nskm)'; off = off + m.nskm;
    m.iZ = off + (1:m.J*m.R)'; off = off + m.J*m.R;
    m.iSd = off + (1:m.R)'; off = off + m.R;
    m.iCor = off + (1:m.nCor)'; off = off + m.nCor;
    m.nPar = off;
end
