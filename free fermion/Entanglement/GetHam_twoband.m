function [H] = GetHam_twoband(nsite)
H = zeros(nsite, nsite);
k = (2*pi)*(1:nsite)/nsite;

Ek1 = cos(k);
Ek2 = cos(2*k);

for i1 = 1:nsite
    for i2 = 1:nsite
        Coe1 = sum(Ek1.*exp(1i*k.*(i1-i2)))/nsite;
        Coe2 = sum(Ek2.*exp(1i*k.*(i1-i2)))/nsite;
        H(i1,i2) = Coe1 + Coe2;
    end
end
end