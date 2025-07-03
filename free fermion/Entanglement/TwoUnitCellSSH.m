function [ks, E] = TwoUnitCellSSH(Num, v, w, u)
%% Hamlitonian used in Shensei's paper
E = zeros(Num+1, 2);
counter = 1;

for k = [-pi:2*pi/Num:pi]
    Hk = [1j*u,v+w*exp(-1j*k);v+w*exp(1j*k),-1j*u];
    [~,Ek] = eig(Hk);
    Ek = diag(Ek);
    E(counter, :) = Ek;
    counter = counter + 1;
end
ks = [-pi:2*pi/Num:pi];
end