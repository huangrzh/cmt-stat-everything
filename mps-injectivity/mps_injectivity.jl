# 1. MPS = psi + psi 
# 2. MPS = psi + g*psi, g is symmetric action and psi is SSB.


using LinearAlgebra,TensorKit,KrylovKit
D = 4;
m = 3

A1 = randn(D,D,m)
A2 = randn(D,D,m)
#@load "As.jld2"
T = zeros(2*D,2*D,m)
T[1:D,1:D,:] = A1
T[D+1:2*D,D+1:2*D,:] = A1 # if we choose A1, no way to get normal MPS??


U = randn(2*D,2*D)
U,_ = svd(U)
@tensor T[-1 -2 -3] := U[-1 1]*T[1 2 -3]*conj(U[-2 2])



@tensor Tm[-1 -11; -2 -12] := T[-1 -2 1]*T[-11 -12 1]
Tm = reshape(Tm, 4*D*D, 4*D*D)
d,v0,info = eigsolve(Tm, 6, :LM) # there are 4 deg states, since D/2*D/2 = D^2/4
@show svdvals(reshape(v0[1],8,8))
@show svdvals(reshape(v0[2],8,8))

v = randn(2*D, 2*D)
v = v*v'
d,u = eigen(v)
p = sortperm(abs.(d))
d = d[p]
u = u[:,p]
d[1:D] .= 0.0
v = u*diagm(d)*u'

v = v[:]
n0 = norm(v)
tol = 1e-10
for i = 1:100
    global v = Tm*v
    n1 = norm(v)
    v = v/n1
    @show i,n1,abs(n1-n0)
    if abs(n1-n0) < tol
        break
    end
    global n0 = n1
end
v = reshape(v, 2*D, 2*D)
vh = sqrt(v)
ivh = inv(vh)
@tensor T1[-1 -2 -3] := ivh[-1 1]*T[1 2 -3]*vh[2 -2];
@tensor rho[-1 -11] := T1[-1 2 3]*T1[-11 2 3];
#=
v = reshape(v0[1], 2*D, 2*D)
v = v + v';
s,u = eigen(v)
p = sortperm(s,rev=true)
s = s[p]
u = u[:,p]

s1 = zeros(2*D)
s1[1:D] .= 1.0
s2 = zeros(2*D)
s2[D+1:end] .= 1.0

P1 = u*diagm(s1)*u'
P2 = u*diagm(s2)*u'


@tensor T1[-1 -2; -3] := P1[-1 1]*T[1 2 -3]*P1[2 -2]
@tensor T2[-1 -2; -3] := P2[-1 1]*T[1 2 -3]*P2[2 -2]

A1 = permutedims(A1, [2,3,1])
A2 = permutedims(A2, [2,3,1])
T1 = permutedims(T1, [2,3,1])
T2 = permutedims(T2, [2,3,1])

A1 = TensorMap(A1, ℂ^D*ℂ^m, ℂ^D)
A2 = TensorMap(A2, ℂ^D*ℂ^m, ℂ^D)
T1 = TensorMap(T1, ℂ^8*ℂ^m, ℂ^8)
T2 = TensorMap(T2, ℂ^8*ℂ^m, ℂ^8)

psi1 = InfiniteMPS([A1])
psi2 = InfiniteMPS([A2])
psi3 = InfiniteMPS([T1])
psi4 = InfiniteMPS([T2])
=#