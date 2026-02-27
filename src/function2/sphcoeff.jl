mutable struct ConvolutionSky{I<:Int, MC<:Matrix{Complex{Float64}}}
    lmax::I
    alm::MC
    realization::I
end

function ConvolutionSky(;
        lmax::Int = 3*128 - 1,
        alm::Matrix{ComplexF64} = fill(1.0 + 1.0im, 2, 2),
        realization::Int = 1
    )
    return ConvolutionSky(
        lmax,
        alm,
        realization
    )
end

function alm_lrange(cs, cc)
    n = sum(2*l + 1 for l in cc.lstart:cc.lstop)
    alm_calc = Matrix{ComplexF64}(undef,  3, n)
    for l in cc.lstart:cc.lstop
        m = 0
        idx_in = alm_idx(l=l, m=m, lmax=cs.lmax)
        idx_out= lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cs.lmax)
        alm_calc[1,idx_out] = cs.alm[idx_in,1]
        alm_calc[2,idx_out] = -(cs.alm[idx_in,2] + 1im*cs.alm[idx_in,3])
        alm_calc[3,idx_out] = -(cs.alm[idx_in,2] - 1im*cs.alm[idx_in,3])
        for m in 1:l
            phase = isodd(m) ? -1.0 : 1.0
            idx_in = alm_idx(l=l, m=m, lmax=cs.lmax)
            idx_out_positive= lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cs.lmax)
            idx_out_negative= lmr_idx(l=l, m=-m, lstart=cc.lstart, mmax=cs.lmax)
            # m positive
            alm_calc[1,idx_out_positive] = cs.alm[idx_in,1]                           # spin0
            alm_calc[2,idx_out_positive] = -(cs.alm[idx_in,2] + 1im*cs.alm[idx_in,3]) # spin2 = -(E +iB)
            alm_calc[3,idx_out_positive] = -(cs.alm[idx_in,2] - 1im*cs.alm[idx_in,3]) # spin2 = -(E +iB)
            # m negative
            alm_calc[1,idx_out_negative] = conj(cs.alm[idx_in,1])*phase             # conj(spin0)
            alm_calc[2,idx_out_negative] = conj(alm_calc[3,idx_out_positive])*phase # conj(spin2)*(-1)^m
            alm_calc[3,idx_out_negative] = conj(alm_calc[2,idx_out_positive])*phase # conj(spin-2)*(-1)^m
        end
    end
    return alm_calc
end

mutable struct ConvolutionBeam{I<:Int, MC<:Matrix{Complex{Float64}}}
    lmax::I
    mmax::I
    blm::MC
end


function ConvolutionBeam(;
    lmax::Int = 3*128-1,
    mmax::Int = 2,
    blm::Matrix{ComplexF64} = fill(1.0 + 1.0im, 2, 2)
    )
    return ConvolutionBeam(lmax, mmax, blm)
end


function blm_lrange(cb, cc)
    n = sum(2*min(l, cc.mmax_calculate) + 1 for l in cc.lstart:cc.lstop)
    blm_calc = Matrix{ComplexF64}(undef, n, 3)
    for l in cc.lstart:cc.lstop
        m = 0
        idx_in = alm_idx(l=l, m=m, lmax=cb.lmax, mmax = cb.mmax)
        idx_out= lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.mmax_calculate)
        blm_calc[idx_out,1] = cb.blm[idx_in,1]
        blm_calc[idx_out,2] = -(cb.blm[idx_in,2] + 1im*cb.blm[idx_in,3]) # spin2 = -(E +iB)
        blm_calc[idx_out,3] = -(cb.blm[idx_in,2] - 1im*cb.blm[idx_in,3]) # spin-2 = -(E -iB)
        for m in 1:min(l, cc.mmax_calculate)
            phase = isodd(m) ? -1.0 : 1.0
            idx_in = alm_idx(l=l, m=m, lmax=cb.lmax, mmax = cb.mmax)
            idx_out_positive= lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.mmax_calculate)
            idx_out_negative= lmr_idx(l=l, m=-m, lstart=cc.lstart, mmax=cc.mmax_calculate)
            # m positive
            blm_calc[idx_out_positive,1] = cb.blm[idx_in,1]
            blm_calc[idx_out_positive,2] = -(cb.blm[idx_in,2] + 1im*cb.blm[idx_in,3]) # spin2 = -(E +iB)
            blm_calc[idx_out_positive,3] = -(cb.blm[idx_in,2] - 1im*cb.blm[idx_in,3]) # spin-2 = -(E -iB)
            # m negative
            blm_calc[idx_out_negative,1] = conj(cb.blm[idx_in,1])*phase # conj(spin0)*(-1)^m 
            blm_calc[idx_out_negative,2] = conj(blm_calc[idx_out_positive,3])*phase # conj(spin-2)*(-1)^m 
            blm_calc[idx_out_negative,3] = conj(blm_calc[idx_out_positive,2])*phase # conj(spin2)*(-1)^m 
        end
    end
    return conj.(blm_calc)
end


@inline function lmr_idx(; l::Int, m::Int, lstart::Int, mmax::Int)
    l ≥ lstart || throw(ArgumentError("l must be >= lstart"))

    mcap = min(l, mmax)
    (-mcap ≤ m ≤ mcap) || throw(ArgumentError("m must be in [-min(l,mmax), min(l,mmax)]"))

    offset = 0
    for k in lstart:l-1
        offset += 2*min(k, mmax) + 1
    end

    return offset + (m + mcap) + 1
end

@inline function alm_idx(; l::Integer, m::Integer, lmax::Integer, mmax::Integer=lmax)
    (0 ≤ m ≤ mmax) || throw(ArgumentError("m must be in [0, mmax]"))
    (m ≤ l ≤ lmax) || throw(ArgumentError("l must satisfy m ≤ l ≤ lmax"))

    # offset for all previous m-blocks
    offset = m * (lmax + 1) - (m * (m - 1)) ÷ 2
    return Int(offset + (l - m) + 1)  # 1-based
end