# For a G-symmetric 1d local model (G is on-site) whose ground state \psi is a SPT (belong to a non-trivial class in H2(G,U(1))), if we use G-symmetric MPS to represent \psi, the MPS must be non-injective. Here the G-symmetric MPS is defined by using G-linear irreps on both the physical and virtual level. 





using MPSKit, MPSKitModels, TensorKit, KrylovKit, Printf, JLD2, MAT, LinearAlgebra, Revise

# H = -σx τz σx - τx σz τx - g * σx σx
function get_hamiltonian(g::Float64, nsite::Int)
    q(i,j) = Irrep[ℤ₂ × ℤ₂](i,j)
    vp = Vect[Z2Irrep⊠Z2Irrep]((0, 0) => 1, (0, 1) => 1, (1, 0) => 1, (1, 1) => 1)
    vi1 = Vect[Z2Irrep⊠Z2Irrep]((1, 0) => 1)
    vi2 = Vect[Z2Irrep⊠Z2Irrep]((0, 1) => 1)

    ZI = TensorMap(zeros, ComplexF64, vp, vp)
    blocks(ZI)[q(1,0)] .= -1
    blocks(ZI)[q(1,1)] .= -1
    blocks(ZI)[q(0,0)] .= +1
    blocks(ZI)[q(0,1)] .= +1

    IZ = TensorMap(zeros, ComplexF64, vp, vp)
    blocks(IZ)[q(1,0)] .= +1
    blocks(IZ)[q(1,1)] .= -1
    blocks(IZ)[q(0,0)] .= +1
    blocks(IZ)[q(0,1)] .= -1

    XIa = TensorMap(ones, ComplexF64, vp * vi1, vp)
    XIb = TensorMap(ones, ComplexF64, vp, vi1 * vp)
    @tensor σσ[-1 -2; -3 -4] := XIa[-1 1; -3] * XIb[-2; 1 -4]

    XZ = TensorMap(zeros, ComplexF64, vp * vi1, vp)
    blocks(XZ)[q(1,0)] .= +1
    blocks(XZ)[q(1,1)] .= -1
    blocks(XZ)[q(0,0)] .= +1
    blocks(XZ)[q(0,1)] .= -1
    @tensor στσ[-1 -2; -3 -4] := XZ[-1 1; -3] * XIb[-2; 1 -4]

    ZX = TensorMap(zeros, ComplexF64, vp, vi2 * vp)
    blocks(ZX)[q(1,0)] .= -1
    blocks(ZX)[q(1,1)] .= -1
    blocks(ZX)[q(0,0)] .= +1
    blocks(ZX)[q(0,1)] .= +1
    IXa = TensorMap(ones, ComplexF64, vp * vi2, vp)
    @tensor τστ[-1 -2; -3 -4] := IXa[-1 1; -3] * ZX[-2; 1 -4]

    H = MPOHamiltonian(-στσ - τστ - g * σσ)
    return repeat(H,nsite), vp
end

function o1o2(o1,o2)
    @tensor o12[-1 -2; -3 -4] := o1[-1 -3]*o2[-2 -4] 
    o12 = reshape(o12, 4, 4)
    return TensorMap(o12, ℂ^4, ℂ^4)
end

function o1o2_(o1,o2)
    @tensor o12[-1 -2; -3 -4] := o1[-1 -3]*o2[-2 -4] 
    return o12
end

function get_ham_nonsym(g,nsite)
    X = ComplexF64[0 1; 1 0]
    Z = ComplexF64[1 0; 0 -1]
    Oi = ComplexF64[1 0; 0 1]
    
    XI = o1o2(X,Oi)
    IX = o1o2(Oi,X)
    XZ = o1o2(X,Z)    
    ZX = o1o2(Z,X)

    σσ = o1o2_(XI, XI)
    στσ = o1o2_(XZ, XI)
    τστ = o1o2_(IX, ZX)

    H = MPOHamiltonian(-στσ - τστ - g * σσ)
    return repeat(H,nsite), ComplexSpace(4)
end


const g = 0.0
const nsite = 4

Ham, vp = get_hamiltonian(g, nsite)
#Ham, vp = get_ham_nonsym(g, nsite)

st0 = InfiniteMPS(fill(vp,nsite), fill(vp,nsite))

st1,_ = find_groundstate(st0, Ham, IDMRG2(trscheme=truncerr(1e-6)))

st2,_ = find_groundstate(st1, Ham, VUMPS(tol=1e-10))