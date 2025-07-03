function C = CorrelationMatrix(H)
% Correlation Matrix
% Half-cut
[U, E] = eig(H);
E = diag(E);
fillind = find(E<0);
Us = U(:, fillind);
Uinv = inv(U);
Uinvs = Uinv(fillind, :);
C = Us*Uinvs;
end