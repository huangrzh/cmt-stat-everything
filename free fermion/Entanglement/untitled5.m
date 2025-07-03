

miu = 1.0;
Ns = 200;
nN = length(Ns);
vS2 = zeros(nN,1);
for iN = 1:nN
    [iN,nN]
    N = Ns(iN);
    Vcc = zeros(N,1);
    Ek2 = zeros(2,N);
    
    k = 2*pi*(1:N)/N - pi;
    k = k(:);
    for ik = 1:N
        Hk = zeros(2,2);
        Hk(1,1) = 0;
        Hk(2,2) = cos(2*k(ik));
        Hk(1,2) = 1./sqrt(2) - cos(k(ik));
        Hk(2,1) = Hk(1,2);
        [Uk,Dk] = eig(Hk);
        Dk = diag(Dk);
        Ek2(:,ik) = Dk;
        Gig = 0;
        for ig=1:2
            if (Dk(ig)<0)
                Gig = Gig + exp(1i*k(ik)*(0:N-1))*abs(Uk(1,ig))^2;
            end
        end
        Vcc = Vcc + Gig(:);
    end
    Vcc = Vcc/N;
    Gcc = zeros(N,N);
    for ir = 1:N
        Gcc(:,ir) = Vcc([N-ir+2:N,1:N-ir+1]);
    end
    Gp = Gcc(1:N/2,1:N/2);
    Gp = Gcc; % sublattice entropy
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
