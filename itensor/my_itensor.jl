using ITensors, KrylovKit, ITensors.HDF5

mutable struct DemoObserver <: AbstractObserver
    energy_tol::Float64
    min_sweep::Int64
    last_energy::Float64

    DemoObserver(energy_tol=0.0, min_sweep=6, last_energy=0.0) = new(energy_tol, min_sweep, last_energy)
end

function ITensors.checkdone!(o::DemoObserver; kwargs...)
    sw = kwargs[:sweep]
    energy = kwargs[:energy]
    dE = real(energy - o.last_energy) / abs(energy)
    #println("@", sw, " ", o.min_sweep, "  ", dE, " ", o.energy_tol)
    if -dE < o.energy_tol && sw > o.min_sweep
        println("E=$energy, dE=$dE,  sweep $sw, Stopping DMRG after sweep $sw\n\n")
        return true
    end
    o.last_energy = energy
    return false
end

function dmrg_eigs(H::MPO, psi0::MPS, b::Int64, neig; kwargs...)
    sw_time = @elapsed begin
        ITensors.check_hascommoninds(siteinds, H, psi0)
        ITensors.check_hascommoninds(siteinds, H, psi0')
        H = ITensors.permute(H, (linkind, siteinds, linkind))
        PH = ProjMPO(H)

        # eigsolve kwargs
        eigsolve_tol::Float64 = get(kwargs, :eigsolve_tol, 1e-14)
        eigsolve_krylovdim::Int = get(kwargs, :eigsolve_krylovdim, 2 * neig)
        eigsolve_maxiter::Int = get(kwargs, :eigsolve_maxiter, 1000)
        eigsolve_verbosity::Int = get(kwargs, :eigsolve_verbosity, 0)
        ishermitian::Bool = get(kwargs, :ishermitian, true)
        eigsolve_which_eigenvalue::Symbol = :SR

        psi = copy(psi0)
        N = length(psi)

        if !isortho(psi) || ITensors.orthocenter(psi) != b
            ITensors.orthogonalize!(psi, b)
        end
        @assert isortho(psi) && ITensors.orthocenter(psi) == b

        position!(PH, psi, b)

        phi = psi[b] * psi[b+1]

        vals, vecs = eigsolve(
            PH,
            phi,
            neig,
            eigsolve_which_eigenvalue;
            ishermitian=ishermitian,
            tol=eigsolve_tol,
            krylovdim=eigsolve_krylovdim,
            maxiter=eigsolve_maxiter
        )
    end
    println("dmrg_exits time = ", sw_time)
    return vals
end


function my_dmrg(Hs::Vector{MPO}, psi0::MPS, sweeps::Sweeps; kwargs...)::Tuple{Number,MPS}
    for H in Hs
        ITensors.check_hascommoninds(siteinds, H, psi0)
        ITensors.check_hascommoninds(siteinds, H, psi0')
    end
    Hs .= ITensors.permute.(Hs, Ref((linkind, siteinds, linkind)))
    PHS = ProjMPOSum(Hs)
    return my_dmrg(PHS, psi0, sweeps; kwargs...)
end

function my_dmrg_onesite(Hs::Vector{MPO}, psi0::MPS, sweeps::Sweeps; kwargs...)::Tuple{Number,MPS}
    for H in Hs
        ITensors.check_hascommoninds(siteinds, H, psi0)
        ITensors.check_hascommoninds(siteinds, H, psi0')
    end
    Hs .= ITensors.permute.(Hs, Ref((linkind, siteinds, linkind)))
    PHS = ProjMPOSum(Hs)
    ITensors.set_nsite!(PHS, 1)
    return my_dmrg_onesite(PHS, psi0, sweeps; kwargs...)
end


function my_dmrg(H::MPO, Ms::Vector{MPS}, psi0::MPS, sweeps::Sweeps; kwargs...)::Tuple{Number,MPS}
    ITensors.check_hascommoninds(siteinds, H, psi0)
    ITensors.check_hascommoninds(siteinds, H, psi0')
    for M in Ms
        ITensors.check_hascommoninds(siteinds, M, psi0)
    end
    H = ITensors.permute(H, (linkind, siteinds, linkind))
    Ms .= ITensors.permute.(Ms, Ref((linkind, siteinds, linkind)))
    weight = get(kwargs, :weight, 1.0)
    PMM = ProjMPO_MPS(H, Ms; weight=weight)
    return my_dmrg(PMM, psi0, sweeps; kwargs...)
end

function my_dmrg_onesite(H::MPO, Ms::Vector{MPS}, psi0::MPS, sweeps::Sweeps; kwargs...)::Tuple{Number,MPS}
    ITensors.check_hascommoninds(siteinds, H, psi0)
    ITensors.check_hascommoninds(siteinds, H, psi0')
    for M in Ms
        ITensors.check_hascommoninds(siteinds, M, psi0)
    end
    H = ITensors.permute(H, (linkind, siteinds, linkind))
    Ms .= ITensors.permute.(Ms, Ref((linkind, siteinds, linkind)))
    weight = get(kwargs, :weight, 1.0)
    PMM = ProjMPO_MPS(H, Ms; weight=weight)
    ITensors.set_nsite!(PMM, 1)
    return my_dmrg_onesite(PMM, psi0, sweeps; kwargs...)
end

function my_dmrg(H::MPO, psi0::MPS, sweeps::Sweeps; kwargs...)::Tuple{Number,MPS}
    ITensors.check_hascommoninds(siteinds, H, psi0)
    ITensors.check_hascommoninds(siteinds, H, psi0')
    # Permute the indices to have a better memory layout
    # and minimize permutations
    H = ITensors.permute(H, (linkind, siteinds, linkind))
    PH = ProjMPO(H)
    return my_dmrg(PH, psi0, sweeps; kwargs...)
end

function my_dmrg_onesite(H::MPO, psi0::MPS, sweeps::Sweeps; kwargs...)::Tuple{Number,MPS}
    ITensors.check_hascommoninds(siteinds, H, psi0)
    ITensors.check_hascommoninds(siteinds, H, psi0')
    # Permute the indices to have a better memory layout
    # and minimize permutations
    H = ITensors.permute(H, (linkind, siteinds, linkind))
    PH = ProjMPO(H)
    ITensors.set_nsite!(PH, 1)
    return my_dmrg_onesite(PH, psi0, sweeps; kwargs...)
end


function my_dmrg(PH, psi0::MPS, sweeps::Sweeps; kwargs...)
    ITensors.@debug_check begin
        checkflux(psi0)
        checkflux(PH)
    end

    which_decomp::Union{String,Nothing} = get(kwargs, :which_decomp, nothing)
    svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")
    obs = get(kwargs, :observer, NoObserver())
    outputlevel::Int = get(kwargs, :outputlevel, 1)

    write_when_maxdim_exceeds::Union{Int,Nothing} = get(
        kwargs, :write_when_maxdim_exceeds, nothing
    )

    # eigsolve kwargs
    eigsolve_tol::Float64 = get(kwargs, :eigsolve_tol, 1e-14)
    eigsolve_krylovdim::Int = get(kwargs, :eigsolve_krylovdim, 3)
    eigsolve_maxiter::Int = get(kwargs, :eigsolve_maxiter, 1)
    eigsolve_verbosity::Int = get(kwargs, :eigsolve_verbosity, 0)
    fname_mps::String = get(kwargs, :fname_mps, "")
    println("fname_mps: ", fname_mps)

    ishermitian::Bool = get(kwargs, :ishermitian, true)
    eigsolve_which_eigenvalue::Symbol = :SR
    if haskey(kwargs, :maxiter)
        error("""maxiter keyword has been replaced by eigsolve_krylovdim.
                 Note: compared to the C++ version of ITensor,
                 setting eigsolve_krylovdim 3 is the same as setting
                 a maxiter of 2.""")
    end

    if haskey(kwargs, :errgoal)
        error("errgoal keyword has been replaced by eigsolve_tol.")
    end

    if haskey(kwargs, :quiet)
        error("quiet keyword has been replaced by outputlevel")
    end

    psi = copy(psi0)
    N = length(psi)

    if !isortho(psi) || ITensors.orthocenter(psi) != 1
        ITensors.orthogonalize!(psi, 1)
    end
    @assert isortho(psi) && ITensors.orthocenter(psi) == 1

    position!(PH, psi, 1)
    energy = 0.0

    for sw in 1:nsweep(sweeps)
        sw_time = @elapsed begin
            maxtruncerr = 0.0

            if !isnothing(write_when_maxdim_exceeds) &&
               maxdim(sweeps, sw) > write_when_maxdim_exceeds
                if outputlevel >= 2
                    println(
                        "write_when_maxdim_exceeds = $write_when_maxdim_exceeds and maxdim(sweeps, sw) = $(maxdim(sweeps, sw)), writing environment tensors to disk",
                    )
                end
                PH = disk(PH)
            end

            for (b, ha) in sweepnext(N)
                b_time = @elapsed begin
                    ITensors.@debug_check begin
                        checkflux(psi)
                        checkflux(PH)
                    end

                    ITensors.@timeit_debug timer "dmrg: position!" begin
                        position!(PH, psi, b)
                    end

                    ITensors.@debug_check begin
                        checkflux(psi)
                        checkflux(PH)
                    end

                    ITensors.@timeit_debug timer "dmrg: psi[b]*psi[b+1]" begin
                        phi = psi[b] * psi[b+1]
                    end

                    ITensors.@timeit_debug timer "dmrg: eigsolve" begin
                        vals, vecs = eigsolve(
                            PH,
                            phi,
                            1,
                            eigsolve_which_eigenvalue;
                            ishermitian=ishermitian,
                            tol=eigsolve_tol,
                            krylovdim=eigsolve_krylovdim,
                            maxiter=eigsolve_maxiter
                        )
                    end
                    energy, phi = vals[1], vecs[1]

                    ortho = ha == 1 ? "left" : "right"

                    drho = nothing
                    if noise(sweeps, sw) > 0.0
                        ITensors.@timeit_debug timer "dmrg: noiseterm" begin
                            # Use noise term when determining new MPS basis
                            drho = noise(sweeps, sw) * noiseterm(PH, phi, ortho)
                        end
                    end

                    ITensors.@debug_check begin
                        checkflux(phi)
                    end

                    ITensors.@timeit_debug timer "dmrg: replacebond!" begin
                        spec = replacebond!(
                            psi,
                            b,
                            phi;
                            maxdim=maxdim(sweeps, sw),
                            mindim=mindim(sweeps, sw),
                            cutoff=cutoff(sweeps, sw),
                            eigen_perturbation=drho,
                            ortho=ortho,
                            normalize=true,
                            which_decomp=which_decomp,
                            svd_alg=svd_alg
                        )
                    end

                    phi = psi[b] * psi[b+1]
                    energy = real(dot(phi, product(PH, phi)))

                    maxtruncerr = max(maxtruncerr, spec.truncerr)

                    ITensors.@debug_check begin
                        checkflux(psi)
                        checkflux(PH)
                    end

                end

                SvN = get_EnS_from_spe(spec)
                if outputlevel >= 2
                    ITensors.@printf(
                        "Twosite: %d,  %d, (%d,%d) %.12f,  %.2E,   %d,   %.3f,  %.3f\n", sw, ha, b, b + 1, energy,
                        spec.truncerr, dim(linkind(psi, b)), SvN, b_time
                    )
                    flush(stdout)
                end
            end
        end
        if outputlevel >= 1
            ITensors.@printf(
                "After sweep %d energy=%.12f maxlinkdim=%d maxerr=%.2E time=%.3f\n",
                sw,
                energy,
                maxlinkdim(psi),
                maxtruncerr,
                sw_time
            )
            flush(stdout)
        end
        isdone = checkdone!(obs; energy=energy, psi=psi, sweep=sw, outputlevel=outputlevel)

        if length(fname_mps) > 0
            fmps = h5open(fname_mps, "w")
            write(fmps, "psi", psi)
            close(fmps)
        end

        isdone && break
    end
    return (energy, psi)
end

function my_dmrg_onesite(PH, psi0::MPS, sweeps::Sweeps; kwargs...)
    #ITensors.check_hascommoninds(siteinds, H, psi0)
    #ITensors.check_hascommoninds(siteinds, H, psi0')
    #H = ITensors.permute(H, (linkind, siteinds, linkind))
    #PH = ProjMPO(H)
    #PH.nsite = 1

    ITensors.@debug_check begin
        checkflux(psi0)
        checkflux(PH)
    end

    which_decomp::Union{String,Nothing} = get(kwargs, :which_decomp, nothing)
    svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")
    obs = get(kwargs, :observer, NoObserver())
    outputlevel::Int = get(kwargs, :outputlevel, 1)

    write_when_maxdim_exceeds::Union{Int,Nothing} = get(
        kwargs, :write_when_maxdim_exceeds, nothing
    )

    # eigsolve kwargs
    eigsolve_tol::Float64 = get(kwargs, :eigsolve_tol, 1e-14)
    eigsolve_krylovdim::Int = get(kwargs, :eigsolve_krylovdim, 3)
    eigsolve_maxiter::Int = get(kwargs, :eigsolve_maxiter, 1)
    eigsolve_verbosity::Int = get(kwargs, :eigsolve_verbosity, 0)
    fname_mps::String = get(kwargs, :fname_mps, "")
    println("fname_mps: ", fname_mps)

    ishermitian::Bool = get(kwargs, :ishermitian, true)
    eigsolve_which_eigenvalue::Symbol = :SR
    if haskey(kwargs, :maxiter)
        error("""maxiter keyword has been replaced by eigsolve_krylovdim.
                 Note: compared to the C++ version of ITensor,
                 setting eigsolve_krylovdim 3 is the same as setting
                 a maxiter of 2.""")
    end

    if haskey(kwargs, :errgoal)
        error("errgoal keyword has been replaced by eigsolve_tol.")
    end

    if haskey(kwargs, :quiet)
        error("quiet keyword has been replaced by outputlevel")
    end

    psi = copy(psi0)
    N = length(psi)

    if !isortho(psi) || ITensors.orthocenter(psi) != 2
        ITensors.orthogonalize!(psi, 2)
    end
    @assert isortho(psi) && ITensors.orthocenter(psi) == 2

    position!(PH, psi, 2)
    energy = 0.0

    for sw in 1:nsweep(sweeps)
        sw_time = @elapsed begin
            maxtruncerr = 0.0

            if !isnothing(write_when_maxdim_exceeds) &&
               maxdim(sweeps, sw) > write_when_maxdim_exceeds
                if outputlevel >= 2
                    println(
                        "write_when_maxdim_exceeds = $write_when_maxdim_exceeds and maxdim(sweeps, sw) = $(maxdim(sweeps, sw)), writing environment tensors to disk",
                    )
                end
                PH = disk(PH)
            end

            for (b, ha) in sweepnext(N)
                if b == 1
                    continue
                end
                ITensors.orthogonalize!(psi, b)

                b_time = @elapsed begin
                    ITensors.@debug_check begin
                        checkflux(psi)
                        checkflux(PH)
                    end

                    ITensors.@timeit_debug timer "dmrg: position!" begin
                        position!(PH, psi, b)
                    end

                    ITensors.@debug_check begin
                        checkflux(psi)
                        checkflux(PH)
                    end

                    ITensors.@timeit_debug timer "dmrg: psi[b]" begin
                        phi = psi[b]
                    end

                    ITensors.@timeit_debug timer "dmrg: eigsolve" begin
                        vals, vecs = eigsolve(
                            x -> product1(x, PH),
                            phi,
                            1,
                            eigsolve_which_eigenvalue;
                            ishermitian=ishermitian,
                            tol=eigsolve_tol,
                            krylovdim=eigsolve_krylovdim,
                            maxiter=eigsolve_maxiter
                        )
                    end
                    energy, phi = vals[1], vecs[1]

                    ortho = ha == 1 ? "left" : "right"

                    drho = nothing
                    if noise(sweeps, sw) > 0.0
                        ITensors.@timeit_debug timer "dmrg: noiseterm" begin
                            # Use noise term when determining new MPS basis
                            drho = noise(sweeps, sw) * noiseterm(PH, phi, ortho)
                        end
                    end

                    ITensors.@debug_check begin
                        checkflux(phi)
                    end
                    phi ./= norm(phi)
                    psi[b] = phi

                    ITensors.@debug_check begin
                        checkflux(psi)
                        checkflux(PH)
                    end
                end


                if outputlevel >= 2
                    SvN, Es = GetEntanglementEntropy(psi, b)
                    ITensors.@printf(
                        "Onesite: %d,  %d, @%d %.12f,  %d,   %.3f    %.3f\n", sw, ha, b, energy,
                        dim(linkind(psi, b)), SvN, b_time
                    )
                    flush(stdout)
                end
            end
        end
        if outputlevel >= 1
            ITensors.@printf(
                "After sweep %d energy=%.12f maxlinkdim=%d maxerr=%.2E time=%.3f\n",
                sw,
                energy,
                maxlinkdim(psi),
                maxtruncerr,
                sw_time
            )
            flush(stdout)
        end

        if typeof(PH) == ProjMPO
            energy = inner(psi, PH.H, psi)
        else
            energy = inner(psi, PH.PH.H, psi)
        end
        isdone = checkdone!(obs; energy=energy, psi=psi, sweep=sw, outputlevel=outputlevel)

        if length(fname_mps) > 0
            fmps = h5open(fname_mps, "w")
            write(fmps, "psi", psi)
            close(fmps)
        end

        isdone && break
    end
    return (energy, psi)
end

#PH = ProjMPO(0, length(H) + 1, 1, H, Vector{ITensor}(undef, length(H)))
function GetEntanglementEntropy(psi::MPS, b::Int64)
    orthogonalize!(psi, b)
    U, S, V = svd(psi[b], (linkind(psi, b - 1), siteind(psi, b)))
    D = dim(S, 1)
    EnSp = zeros(D)
    for n = 1:D
        EnSp[n] = S[n, n]^2
    end
    EnSp = EnSp / sum(EnSp)
    SvN = -sum(EnSp .* log.(EnSp))

    return SvN, EnSp
end

function product1(v::ITensor, P::ProjMPO)::ITensor
    Pv = ITensors.contract(P, v)
    if ITensors.order(Pv) != ITensors.order(v)
        error(
            string(
                "The order of the ProjMPO-ITensor product P*v is not equal to the order of the ITensor v, ",
                "this is probably due to an index mismatch.\nCommon reasons for this error: \n",
                "(1) You are trying to multiply the ProjMPO with the $(nsite(P))-site wave-function at the wrong position.\n",
                "(2) `orthogonalize!` was called, changing the MPS without updating the ProjMPO.\n\n",
                "P*v inds: $(inds(Pv)) \n\n",
                "v inds: $(inds(v))",
            ),
        )
    end
    return ITensors.noprime(Pv)
end



function product1(v::ITensor, P::ITensors.ProjMPO_MPS)::ITensor
    Pv = product1(v, P.PH)
    for p in P.pm
        Pv += P.weight * product1(v, p)
    end
    return Pv
end


function product1(v::ITensor, P::ITensors.ProjMPS)::ITensor
    #if nsite(P) != 2
    #    error("Only two-site ProjMPS currently supported")
    #end

    if ITensors.nsite(P) == 2
        Lpm = dag(prime(P.M[P.lpos+1], "Link"))
        !isnothing(lproj(P)) && (Lpm *= lproj(P))

        Rpm = dag(prime(P.M[P.rpos-1], "Link"))
        !isnothing(rproj(P)) && (Rpm *= rproj(P))

        pm = Lpm * Rpm

        pv = scalar(pm * v)

        Mv = pv * dag(pm)

        return noprime(Mv)
    else
        pm = dag(prime(P.M[P.lpos+1], "Link"))
        !isnothing(lproj(P)) && (pm *= lproj(P))

        #Rpm = dag(prime(P.M[P.rpos-1], "Link"))
        !isnothing(rproj(P)) && (pm *= rproj(P))

        pv = scalar(pm * v)

        Mv = pv * dag(pm)

        return noprime(Mv)
    end
end





#= function ITensors.truncate!(P::Vector{Float64}; kwargs...)::Tuple{Float64,Float64}
    # Keyword argument deprecations
    use_absolute_cutoff = false
    if haskey(kwargs, :absoluteCutoff)
        @warn "In truncate!, keyword argument absoluteCutoff is deprecated in favor of use_absolute_cutoff"
        use_absolute_cutoff = get(kwargs, :absoluteCutoff, use_absolute_cutoff)
    end
    use_relative_cutoff = true
    if haskey(kwargs, :doRelCutoff)
        @warn "In truncate!, keyword argument doRelCutoff is deprecated in favor of use_relative_cutoff"
        use_relative_cutoff = get(kwargs, :doRelCutoff, use_relative_cutoff)
    end

    maxdim::Int = min(get(kwargs, :maxdim, length(P)), length(P))
    mindim::Int = max(get(kwargs, :mindim, 1), 1)
    cutoff::Float64 = max(get(kwargs, :cutoff, 0.0), 0.0)
    use_absolute_cutoff::Bool = get(kwargs, :use_absolute_cutoff, use_absolute_cutoff)
    use_relative_cutoff::Bool = get(kwargs, :use_relative_cutoff, use_relative_cutoff)

    origm = length(P)
    docut = 0.0

    #if P[1] <= 0.0
    #  P[1] = 0.0
    #  resize!(P, 1)
    #  return 0.0, 0.0
    #end

    if origm == 1
        docut = abs(P[1]) / 2
        return 0.0, docut
    end

    s = sign(P[1])
    s < 0 && (P .*= s)

    #Zero out any negative weight
    for n in origm:-1:1
        (P[n] >= 0.0) && break
        P[n] = 0.0
    end

    n = origm
    truncerr = 0.0
    while n > maxdim
        truncerr += P[n]
        n -= 1
    end

    if use_absolute_cutoff
        #Test if individual prob. weights fall below cutoff
        #rather than using *sum* of discarded weights
        while P[n] <= cutoff && n > mindim
            truncerr += P[n]
            n -= 1
        end
    else
        scale = 1.0
        if use_relative_cutoff
            scale = sum(P)
            (scale == 0.0) && (scale = 1.0)
        end

        #Continue truncating until *sum* of discarded probability 
        #weight reaches cutoff reached (or m==mindim)
        while (truncerr + P[n] <= cutoff * scale) && (n > mindim)
            truncerr += P[n]
            n -= 1
        end

        truncerr /= scale
    end

    if n < 1
        n = 1
    end

    if n < origm
        docut = (P[n] + P[n+1]) / 2
    end

    s < 0 && (P .*= s)
    resize!(P, n)

    return truncerr, docut
end



function mycorrelation(psi::MPS, Op1::AbstractString, Op2::AbstractString, i, end_site; kwargs...)
    N = length(psi)
    ElT = ITensors.promote_itensor_eltype(psi)

    start_site = i

    psi = copy(psi)
    orthogonalize!(psi, start_site)
    norm2_psi = norm(psi[start_site])^2

    s = siteinds(psi)
    onsiteOp = "$Op1*$Op2"
    fermionic2 = has_fermion_string(Op2, s[1])
    if !using_auto_fermion() && fermionic2
        Op1 = "$Op1*F"
    end

    if start_site == 1
        L = ITensor(1.0)
    else
        lind = commonind(psi[start_site], psi[start_site-1])
        L = delta(dag(lind), lind')
    end

    Li = L * psi[i]

    # Get j > i correlations
    Li = Li * op(Op1, s, i) * dag(prime(psi[i]))
    for j in (i+1):end_site-1
        cj = j - start_site + 1
        lind = commonind(psi[j], Li)
        Li *= psi[j]

        if !using_auto_fermion() && fermionic2
            Li *= op("F", s, j) * dag(prime(psi[j]))
        else
            Li *= dag(prime(psi[j], "Link"))
        end
    end
    val = Li * op(Op2, s, end_site) * dag(prime(prime(psi[end_site], "Site"), lind))
    C = scalar(val) / norm2_psi

    return C
end



function my_correlation_matrix(psi::MPS, Op1::AbstractString, Op2::AbstractString; kwargs...)
    N = length(psi)
    ElT = ITensors.promote_itensor_eltype(psi)

    site_range::UnitRange{Int} = get(kwargs, :site_range, 1:N)
    start_site = first(site_range)
    end_site = last(site_range)

    psi = copy(psi)
    ITensors.orthogonalize!(psi, start_site)
    norm2_psi = norm(psi[start_site])^2

    s = siteinds(psi)
    onsiteOp = "$Op1*$Op2"
    fermionic2 = has_fermion_string(Op2, s[1])
    if !ITensors.using_auto_fermion() && fermionic2
        Op1 = "$Op1*F"
    end

    # Nb = size of block of correlation matrix
    Nb = end_site - start_site + 1

    C = zeros(ElT, Nb, Nb)

    if start_site == 1
        L = ITensor(1.0)
    else
        lind = commonind(psi[start_site], psi[start_site-1])
        L = delta(dag(lind), lind')
    end

    for i in start_site:(end_site-1)
        println("----------------------------\n i = ", i)
        ci = i - start_site + 1

        Li = L * psi[i]

        # Get j == i diagonal correlations
        rind = commonind(psi[i], psi[i+1])
        C[ci, ci] = scalar(Li * op(onsiteOp, s, i) * prime(dag(psi[i]), not(rind))) / norm2_psi

        # Get j > i correlations
        Li = Li * op(Op1, s, i) * dag(prime(psi[i]))
        for j in (i+1):end_site
            println("j = ", j)
            cj = j - start_site + 1
            lind = commonind(psi[j], Li)
            Li *= psi[j]

            val = Li * op(Op2, s, j) * dag(prime(prime(psi[j], "Site"), lind))
            C[ci, cj] = scalar(val) / norm2_psi
            C[cj, ci] = conj(C[ci, cj])

            if !ITensors.using_auto_fermion() && fermionic2
                Li *= op("F", s, j) * dag(prime(psi[j]))
            else
                Li *= dag(prime(psi[j], "Link"))
            end
        end
        L *= psi[i] * dag(prime(psi[i], "Link"))
    end

    # Get last diagonal element of C
    i = end_site
    lind = commonind(psi[i], psi[i-1])
    C[Nb, Nb] =
        scalar(L * psi[i] * op(onsiteOp, s, i) * prime(prime(dag(psi[i]), "Site"), lind)) /
        norm2_psi

    return C
end
 =#
#=
  function LinearAlgebra.svd(T::DenseTensor{ElT,2,IndsT}; kwargs...) where {ElT,IndsT}
    truncate = haskey(kwargs, :maxdim) || haskey(kwargs, :cutoff)

    #
    # Keyword argument deprecations
    #
    use_absolute_cutoff = false
    if haskey(kwargs, :absoluteCutoff)
      @warn "In svd, keyword argument absoluteCutoff is deprecated in favor of use_absolute_cutoff"
      use_absolute_cutoff = get(kwargs, :absoluteCutoff, use_absolute_cutoff)
    end

    use_relative_cutoff = true
    if haskey(kwargs, :doRelCutoff)
      @warn "In svd, keyword argument doRelCutoff is deprecated in favor of use_relative_cutoff"
      use_relative_cutoff = get(kwargs, :doRelCutoff, use_relative_cutoff)
    end

    if haskey(kwargs, :fastsvd) || haskey(kwargs, :fastSVD)
      error(
        "In svd, fastsvd/fastSVD keyword arguments are removed in favor of alg, see documentation for more details.",
      )
    end

    maxdim::Int = get(kwargs, :maxdim, minimum(dims(T)))
    mindim::Int = get(kwargs, :mindim, 1)
    cutoff::Float64 = get(kwargs, :cutoff, 0.0)
    use_absolute_cutoff::Bool = get(kwargs, :use_absolute_cutoff, use_absolute_cutoff)
    use_relative_cutoff::Bool = get(kwargs, :use_relative_cutoff, use_relative_cutoff)
    alg::String = get(kwargs, :alg, "divide_and_conquer")

    #@timeit_debug timer "dense svd" begin
    if alg == "divide_and_conquer"
      MUSV = svd_catch_error(matrix(T); alg=LinearAlgebra.DivideAndConquer())

      if isnothing(MUSV)
        println("Start using QRIteration method for SVD")
        MUSV = svd_catch_error(matrix(T); alg=LinearAlgebra.QRIteration())
      end

      if isnothing(MUSV)
        println("Start using svd_recursive method for SVD")
        MUSV = svd_recursive(matrix(T))
      end

    elseif alg == "qr_iteration"
      MUSV = svd_catch_error(matrix(T); alg=LinearAlgebra.QRIteration())
    elseif alg == "recursive"
      MUSV = svd_recursive(matrix(T))
    else
      error(
        "svd algorithm $alg is not currently supported. Please see the documentation for currently supported algorithms.",
      )
    end
    if isnothing(MUSV)
      if any(isnan, T)
        println("SVD failed, the matrix you were trying to SVD contains NaNs.")
      else
        println(lapack_svd_error_message(alg))
      end
      return nothing
    end
    MU, MS, MV = MUSV
    conj!(MV)
    #end # @timeit_debug

    P = MS .^ 2
    if truncate
      truncerr, _ = truncate!(
        P;
        mindim=mindim,
        maxdim=maxdim,
        cutoff=cutoff,
        use_absolute_cutoff=use_absolute_cutoff,
        use_relative_cutoff=use_relative_cutoff,
      )
    else
      truncerr = 0.0
    end
    spec = Spectrum(P, truncerr)
    dS = length(P)
    if dS < length(MS)
      MU = MU[:, 1:dS]
      resize!(MS, dS)
      MV = MV[:, 1:dS]
    end

    # Make the new indices to go onto U and V
    u = eltype(IndsT)(dS)
    v = eltype(IndsT)(dS)
    Uinds = IndsT((ind(T, 1), u))
    Sinds = IndsT((u, v))
    Vinds = IndsT((ind(T, 2), v))
    U = tensor(Dense(vec(MU)), Uinds)
    S = tensor(Diag(MS), Sinds)
    V = tensor(Dense(vec(MV)), Vinds)
    return U, S, V, spec
  end  
  =#