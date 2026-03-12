mutable struct ConvolutionSky{I<:Int, AC<:Array{ComplexF64,3}}
    numofsky::I
    lmax::I
    alm::AC   # (Nsky, 3, Nalm)
end

function ConvolutionSky(;
    numofsky::Int = 1,
    lmax::Int = 3*128 - 1,
    alm::Union{Nothing, Array{ComplexF64,3}} = nothing,
)
    ncoeff = alm_idx(lmax = lmax, mmax = lmax, l=lmax, m=lmax)
    alm_arr = isnothing(alm) ?
        fill(1.0 + 1.0im, numofsky, 3, ncoeff) :
        alm

    size(alm_arr, 1) == numofsky || throw(ArgumentError("alm first axis must be 3"))
    size(alm_arr, 2) == 3 || throw(ArgumentError("alm second axis must be Nalm=$ncoeff"))
    size(alm_arr, 3) == ncoeff || throw(ArgumentError("alm third axis must match numofsky=$numofsky"))

    return ConvolutionSky(numofsky, lmax, alm_arr)
end

function slice_spin_alm_by_l(cs, cc)
    #n = alm_idx(lmax = cs.lmax, mmax = cs.lmax, l=cs.lmax, m=cs.lmax)
    n = sum(2*l + 1 for l in cc.lstart:cc.lstop)
    alm_calc = Array{ComplexF64,3}(undef, cs.numofsky, 3, n)

    for i in 1:cs.numofsky
        for l in cc.lstart:cc.lstop
            m = 0
            idx_in = alm_idx(l=l, m=m, lmax=cs.lmax)
            idx_out = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cs.lmax)
            alm_calc[i, 1, idx_out] = cs.alm[i, 1, idx_in] # spin0
            alm_calc[i, 2, idx_out] = -(cs.alm[i, 2, idx_in] + 1im*cs.alm[i, 3, idx_in])
            alm_calc[i, 3, idx_out] = -(cs.alm[i, 2, idx_in] - 1im*cs.alm[i, 3, idx_in])

            for m in 1:l
                phase = isodd(m) ? -1.0 : 1.0
                idx_in = alm_idx(l=l, m=m, lmax=cs.lmax)
                idx_out_positive = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cs.lmax)
                idx_out_negative = lmr_idx(l=l, m=-m, lstart=cc.lstart, mmax=cs.lmax)
                # m positive
                alm_calc[i, 1, idx_out_positive] = cs.alm[i, 1, idx_in]                               # spin0
                alm_calc[i, 2, idx_out_positive] = -(cs.alm[i, 2, idx_in] + 1im*cs.alm[i, 3, idx_in]) # spin2 = -(E +iB)
                alm_calc[i, 3, idx_out_positive] = -(cs.alm[i, 2, idx_in] - 1im*cs.alm[i, 3, idx_in]) # spin2 = -(E +iB)
                # m negative
                alm_calc[i, 1, idx_out_negative] = conj(cs.alm[i, 1, idx_in])*phase              # conj(spin0)
                alm_calc[i, 2, idx_out_negative] = conj(alm_calc[i, 3, idx_out_positive])*phase   # conj(spin2)*(-1)^m
                alm_calc[i, 3, idx_out_negative] = conj(alm_calc[i, 2, idx_out_positive])*phase   # conj(spin-2)*(-1)^m
            end
        end
    end

    return alm_calc
end
#=
function slice_spin_alm_by_l(cs, cc,numofsky)
    n = alm_idx(lmax = cs.lmax, mmax = cs.lmax, l=cs.lmax, m=cs.lmax)
    alm_calc = Matrix{ComplexF64}(undef, n, 3, )

        for l in cc.lstart:cc.lstop
            m = 0
            idx_in = alm_idx(l=l, m=m, lmax=cs.lmax)
            idx_out = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cs.lmax)
            alm_calc[idx_out, 1] = cs.alm[1, idx_in, numofsky] # spin0
            alm_calc[idx_out, 2] = -(cs.alm[2, idx_in, numofsky] + 1im*cs.alm[3, idx_in, numofsky])
            alm_calc[idx_out, 3] = -(cs.alm[2, idx_in, numofsky] - 1im*cs.alm[3, idx_in, numofsky])

            for m in 1:l
                phase = isodd(m) ? -1.0 : 1.0
                idx_in = alm_idx(l=l, m=m, lmax=cs.lmax)
                idx_out_positive = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cs.lmax)
                idx_out_negative = lmr_idx(l=l, m=-m, lstart=cc.lstart, mmax=cs.lmax)
                # m positive
                alm_calc[idx_out_positive, 1] = cs.alm[1, idx_in, numofsky]                               # spin0
                alm_calc[idx_out_positive, 2] = -(cs.alm[2, idx_in, numofsky] + 1im*cs.alm[3, idx_in, numofsky]) # spin2 = -(E +iB)
                alm_calc[idx_out_positive, 3] = -(cs.alm[2, idx_in, numofsky] - 1im*cs.alm[3, idx_in, numofsky]) # spin2 = -(E +iB)
                # m negative
                alm_calc[idx_out_negative, 1] = conj(cs.alm[1, idx_in, numofsky])*phase              # conj(spin0)
                alm_calc[idx_out_negative, 2] = conj(alm_calc[idx_out_positive, 3])*phase   # conj(spin2)*(-1)^m
                alm_calc[idx_out_negative, 3] = conj(alm_calc[idx_out_positive, 2])*phase   # conj(spin-2)*(-1)^m
            end
        end
    

    return alm_calc
end
=#

mutable struct ConvolutionBeam{I<:Int, AC<:Array{ComplexF64,3}}
    numofbeams::I
    lmax::I
    mmax::I
    blm::AC
end

function ConvolutionBeam(;
    numofbeams::Int = 1,
    lmax::Int = 3*128 - 1,
    mmax::Int = 2,
    blm::Union{Nothing, Array{ComplexF64,3}} = nothing,
)
    mmax <= lmax || throw(ArgumentError("mmax must be <= lmax"))

    ncoeff = alm_idx(lmax = lmax, mmax = mmax, l=lmax, m=mmax)
    blm_arr = isnothing(blm) ?
        fill(1.0 + 1.0im,ncoeff, 3, numofbeams) :
        blm

    size(blm_arr, 1) == ncoeff     || throw(ArgumentError("blm first axis must be ncoeff=$ncoeff"))
    size(blm_arr, 2) == 3          || throw(ArgumentError("blm second axis must be 3"))
    size(blm_arr, 3) == numofbeams || throw(ArgumentError("blm third axis must be numofbeams=$numofbeams"))

    return ConvolutionBeam(numofbeams, lmax, mmax, blm_arr)
end

function slice_spin_blm_by_l(cb, cc)
    #n = alm_idx(lmax = cb.lmax, mmax = cb.mmax_calculate, l=cb.lmax, m=cc.mmax_calculate)
    n = sum(2*min(l, cb.mmax) + 1 for l in cc.lstart:cc.lstop)
    blm_calc = Array{ComplexF64,3}(undef, n, 3, cb.numofbeams)

    for i in 1:cb.numofbeams
        for l in cc.lstart:cc.lstop
            m = 0
            idx_in = alm_idx(l=l, m=m, lmax=cb.lmax, mmax = cb.mmax)
            idx_out = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.mmax_calculate)
            blm_calc[idx_out, 1, i] = cb.blm[idx_in, 1, i]
            blm_calc[idx_out, 2, i] = -(cb.blm[idx_in, 2, i] + 1im*cb.blm[idx_in, 3, i]) # spin2 = -(E +iB)
            blm_calc[idx_out, 3, i] = -(cb.blm[idx_in, 2, i] - 1im*cb.blm[idx_in, 3, i]) # spin-2 = -(E -iB)

            for m in 1:min(l, cc.mmax_calculate)
                phase = isodd(m) ? -1.0 : 1.0
                idx_in = alm_idx(l=l, m=m, lmax=cb.lmax, mmax = cb.mmax)
                idx_out_positive = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.mmax_calculate)
                idx_out_negative = lmr_idx(l=l, m=-m, lstart=cc.lstart, mmax=cc.mmax_calculate)
                # m positive
                blm_calc[idx_out_positive, 1, i] = cb.blm[idx_in, 1, i]
                blm_calc[idx_out_positive, 2, i] = -(cb.blm[idx_in, 2, i] + 1im*cb.blm[idx_in, 3, i]) # spin2 = -(E +iB)
                blm_calc[idx_out_positive, 3, i] = -(cb.blm[idx_in, 2, i] - 1im*cb.blm[idx_in, 3, i]) # spin-2 = -(E -iB)
                # m negative
                blm_calc[idx_out_negative, 1, i] = conj(cb.blm[idx_in, 1, i])*phase            # conj(spin0)*(-1)^m
                blm_calc[idx_out_negative, 2, i] = conj(blm_calc[idx_out_positive, 3, i])*phase # conj(spin-2)*(-1)^m
                blm_calc[idx_out_negative, 3, i] = conj(blm_calc[idx_out_positive, 2, i])*phase # conj(spin2)*(-1)^m
            end
        end
    end

    return conj.(blm_calc)
end

function slice_spin_blm_by_l(cb, cc, numofbeams)
    n = alm_idx(lmax = cb.lmax, mmax = cb.mmax, l=cb.lmax, m=cb.mmax)
    blm_calc = Matrix{ComplexF64}(undef, 3, n)

    for l in cc.lstart:cc.lstop
        m = 0
        idx_in = alm_idx(l=l, m=m, lmax=cb.lmax, mmax = cb.mmax)
        idx_out = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.mmax_calculate)
        blm_calc[idx_out, 1] = cb.blm[idx_in, 1, numofbeams]
        blm_calc[idx_out, 2] = -(cb.blm[idx_in, 2, numofbeams] + 1im*cb.blm[idx_in, 3, numofbeams]) # spin2 = -(E +iB)
        blm_calc[idx_out, 3] = -(cb.blm[idx_in, 2, numofbeams] - 1im*cb.blm[idx_in, 3, numofbeams]) # spin-2 = -(E -iB)

        for m in 1:min(l, cc.mmax_calculate)
            phase = isodd(m) ? -1.0 : 1.0
            idx_in = alm_idx(l=l, m=m, lmax=cb.lmax, mmax = cb.mmax)
            idx_out_positive = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.mmax_calculate)
            idx_out_negative = lmr_idx(l=l, m=-m, lstart=cc.lstart, mmax=cc.mmax_calculate)
            # m positive
            blm_calc[idx_out_positive, 1] = cb.blm[idx_in, 1, numofbeams]
            blm_calc[idx_out_positive, 2] = -(cb.blm[idx_in, 2, numofbeams] + 1im*cb.blm[idx_in, 3, numofbeams]) # spin2 = -(E +iB)
            blm_calc[idx_out_positive, 3] = -(cb.blm[idx_in, 2, numofbeams] - 1im*cb.blm[idx_in, 3, numofbeams]) # spin-2 = -(E -iB)
            # m negative
            blm_calc[idx_out_negative, 1] = conj(cb.blm[idx_in, 1, numofbeams])*phase            # conj(spin0)*(-1)^m
            blm_calc[idx_out_negative, 2] = conj(blm_calc[idx_out_positive, 3])*phase # conj(spin-2)*(-1)^m
            blm_calc[idx_out_negative, 3] = conj(blm_calc[idx_out_positive, 2])*phase # conj(spin2)*(-1)^m
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
