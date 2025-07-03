function [H, Ek, k] = GetHam(Num,z)
L = 2*Num;
H = zeros(L, L);
k = (2*pi)*(1:Num)/Num - pi;
Ek = abs(k).^z;
for i1 = 1:Num
    for i2 = 1:Num
        Coe = sum(Ek.*exp(1i*k.*(i1-i2)))/Num;
        a1 = 2*i1-1;
        a2 = 2*i2-1;
        b1 = 2*i1;
        b2 = 2*i2;
        H(a1,a2) = +0.48*Coe;
        H(b1,b2) = -0.48*Coe;
        H(a1,b2) = +0.64*Coe;
        H(b1,a2) = +0.64*Coe;
    end
end
end