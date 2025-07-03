function [H, k1, Ek1, k2, Ek2] = GetHam_power(nsite,z)
H = zeros(nsite, nsite);
k = (2*pi)*(1:nsite)/nsite;

k1 = k(k>pi);
k2 = setdiff(k,k1);
Ek1 = zeros(size(k1));
for ik1 = 1:length(k1)
    if k1(ik1)-pi*3/2<0
        Ek1(ik1) = - abs(k1(ik1)-pi*3/2).^z;
    else
        Ek1(ik1) = (k1(ik1)-pi*3/2)^z;
    end
end

Ek2 = zeros(size(k2));
for ik2 = 1:length(k2)
    if k2(ik2)-pi/2<0
        Ek2(ik2) = abs(k2(ik2)-pi/2)^z;
    else
        Ek2(ik2) = -(k2(ik2)-pi/2)^z;
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