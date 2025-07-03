

miu = 1.0;
Ns = 20:20:400;
nN = length(Ns);
vS2 = zeros(nN,1);
t1 = 1;
t2 = 0.2;
miu = 0;

for iN = 1:nN
    [iN,nN]
    N = Ns(iN);
    ks = 2*pi*(1:N)/N - pi;
    ks = ks(:);
    Ek2 = 2*t1*cos(ks) + 2*t2*cos(2*ks) + miu; % t1, t2
    iks = find(Ek2<0);
    Vcc = 0;
    for ik = iks'
        k = ks(ik);
        rs = 0:N-1;
        Vcc = Vcc + exp(-1i*k*rs);
    end
    Vcc = Vcc/N;
    Gcc = zeros(N,N);
    for ir = 1:N
        Gcc(:,ir) = Vcc([N-ir+2:N,1:N-ir+1]);
    end
    Gp = Gcc(1:N/2,1:N/2);
    Ls = eig(Gp);
    Ls = abs(real(Ls));
    nLs = length(Ls);
    S = 0.;
    for iL = 1:nLs
        if not(or(abs(Ls(iL)-1)<1e-14, abs(Ls(iL))<1e-14))
            S = S - Ls(iL)*log(Ls(iL)) - (1-Ls(iL))*log(1-Ls(iL));
        end
    end
    vS2(iN) = S;
end
