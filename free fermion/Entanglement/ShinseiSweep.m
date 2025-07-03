Num = 30;
M = [2:2:30];
w1 = 0.8;
w2 = 0.8;
v1 = 0.7;
v2 = 0.7;
u = 0.5;
filling = 0.5;
S = zeros(numel(M),1);
%flux = 1.0e-8;
flux = 0;
for m = 1:numel(M)
    [S(m), ~] = ShinseiEntropy(v1, v2, w1, w2, u, Num, filling, M(m), flux);
    %[2*Num(m), gamma(n), S(m,n)]
end



Ceffa = zeros(1,1);
Ceffb = zeros(1,1);
R2 = zeros(1,1);

for l = 1
    p1 = polyfit(log(M'), S, 1);
    Ceffa(l,1) = p1(1);
    Ceffb(l,1) = p1(2);      
    R2(l,1) = norm(p1(1)*log(M')+p1(2)-mean(S))^2/norm(S...
        - mean(S))^2;
end
figure(1)
plot(M', real(S), 's - b', 'LineWidth', 2);
figure(2)
plot(M', imag(S), 'o - r', 'LineWidth', 2);
%semilogx(Num*2, S, 's - b', 'LineWidth', 2);