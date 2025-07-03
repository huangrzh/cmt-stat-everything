function [S,Ls,Gij,C] = GetEntropy(H, ids)

Gij = CorrelationMatrix(H);
C = Gij(ids, ids);
Ls = eig(C);
Ls = abs(real(Ls));
nLs = length(Ls);
S = 0.;
for iL = 1:nLs
    if not(or(abs(Ls(iL)-1)<1e-14, abs(Ls(iL))<1e-14))
        S = S - Ls(iL)*log(Ls(iL)) - (1-Ls(iL))*log(1-Ls(iL));
    end
end
end