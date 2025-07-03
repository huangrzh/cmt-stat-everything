function [w, deltaE] = findgap
%% v1=v2=v, w=v2, w^2+v^2-u^2=0
% v1 = 0.8;
% v2 = 0.8;
% v3 = 1;
% u = sqrt(2);
flux = 0;
w = 1;
Num = 1000;
v = 0.7;
w = 0.8;
u = 0.5;
for m = 1:numel(w)
    %[ks, Ek] = Hamiltoniank(Num, v1, v2, v3, w(m), u);
    [ks, Ek] = TwoUnitCellSSH(Num, v, w, u);
    figure(m)
    plot(ks, real(Ek), '. b ')
    hold on
    plot(ks, imag(Ek), '. r ')
end



end