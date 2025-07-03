function [H, k1, Ek1, k2, Ek2] = GetHam_cos(nsite)
H = zeros(nsite, nsite);
k = (2*pi)*(1:nsite)/nsite;

k1 = k(k>pi);
k2 = setdiff(k,k1);
Ek1 = zeros(size(k1));
for ik1 = 1:length(k1)
    if k1(ik1)-pi*7/4<0
        Ek1(ik1) = cos(abs(k1(ik1)-pi*7/4))-1;
    else
        Ek1(ik1) = 1-cos(k1(ik1)-pi*7/4);
    end
end

Ek2 = zeros(size(k2));
for ik2 = 1:length(k2)
    if k2(ik2)-pi/4<0
        Ek2(ik2) = 1-cos(abs(k2(ik2)-pi/4));
    else
        Ek2(ik2) = cos((k2(ik2)-pi/4))-1;
    end
end

for i1 = 1:nsite
    for i2 = 1:nsite
        Coe1 = sum(Ek1.*exp(1i*k1.*(i1-i2)))/nsite;
        Coe2 = sum(Ek2.*exp(1i*k2.*(i1-i2)))/nsite;
        H(i1,i2) = Coe1 + Coe2;
    end
end
end