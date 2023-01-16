% Diagonalize the free fermion Hamiltonian and displace operator
% simultaneously

N= 100;

H = zeros(N,N);
T = H;
for i=1:N
    j = mod(i,N)+1;
    
    H(j,i) = -1.0;
    H(i,j) = -1.0;
    T(j,i) = 1.0;
end

[U,Ek,Ks] = simdiag(H,T);

f = @(x,y) norm(U*x*U'-y);
disp([f(Ek,H), f(Ks,T), norm(imag(Ek))]);

Ks = angle(diag(Ks))/pi;
[Ks,id] = sort(Ks);
Ek = real(diag(Ek));
Ek = Ek(id);

figure; hold on; box on;
plot(Ks, Ek, 'o-');
plot(Ks, zeros(N,1));
xlabel('k');
ylabel('E_k');
