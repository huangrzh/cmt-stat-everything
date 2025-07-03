
vL = 20;
nL = numel(vL);
vS = zeros(nL,1);
vE = cell(nL,1);
vG = vE;
vGc = vE;
vLs = vE;
miu = 1.0;
for id = 1:nL
    id
    nsite = vL(id);
    ids = 1:2:nsite;
    H = GetHam_1d_2site4(nsite,miu);
    %H = GetHam_1d(1, 0.5, 0, nsite);
    [vS(id),vLs{id},vG{id},vGc{id}] = GetEntropy(H,ids);
end