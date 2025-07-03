
clear;
miu = 0;
vL = 20:10:200;
nL = numel(vL);
zs = 1:0.2:2;
zs = 100;
nz = length(zs);
vS = zeros(nL,nz);
vE = cell(nL,nz);
vLs = vE;

t1 = 1i;
ta = 1;
tb = -1;
tab = 1;
ua = 0.1;
ub = -0.1;
for iz = 1:nz
    z = zs(iz);
    for id = 1:nL
        [id,iz]
        nsite = vL(id);
        %[H,Ek,k] = GetHam(nsite, z);
        %[H] = GetHam_1d(nsite, miu);
        [H, k1, Ek1, k2, Ek2] = GetHam_power_bose(nsite,z);
        %[H, k1, Ek1, k2, Ek2] = GetHam_cos(nsite);
        %[H] = GetHam_1d_2site(nsite, t1, ta, tb, tab, ua, ub);
        %[H] = GetHam_1d_2site2(nsite, ta, tb, ua, ub);
        CT = inv(H);
        M = nsite/2;
        C = CT(1:M, 1:M);
        Ls = eig(C);
        
        %ve = real(ve);
        
        Ls = abs(real(Ls));
        nLs = length(Ls);
        S = 0.;
        for iL = 1:nLs
            if not(or(abs(Ls(iL)-1)<1e-14, abs(Ls(iL))<1e-14))
                S = S + (0.5+Ls(iL))*log(0.5+Ls(iL)) - (Ls(iL)-0.5)*log(Ls(iL)-0.5);
            end
        end
        %E = sort([Ek1,Ek2]);
        vS(id,iz) = S;
        %vE{id,iz} = E;
        vLs{id,iz} = Ls;
    end
end