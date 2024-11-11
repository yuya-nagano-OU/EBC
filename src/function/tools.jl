function alm_idx(l, m::Integer, lmax::Integer)
    return Int(m * (2 * lmax + 1 - m) // 2 + l)+1
end

function for_healpy_order(l, m::Integer, lmax::Integer)
    a = 1
    if m != 0
        for i in 1:m
            a += lmax + 2 - i
        end
    end
    a += l - m
    return a
end

function unique_theta_num(num, cp)
    # This function calculates the unique theta number range for a given pixel number and nside parameter.
    #
    # # Arguments
    # - `num::Int`: The pixel number for which the unique theta number range is to be calculated.
    # - `nside::Int`: The nside parameter, which determines the resolution of the HEALPix grid.
    #
    # # Returns
    # - `Tuple{Int, Int}`: A tuple containing the start and stop indices of the unique theta number range.
    #
    # The function first calculates the total number of pixels (npix) for the given nside using the `nside2npix` function.
    # It then determines the start and stop indices based on the value of `num` relative to `nside` and `3nside`.
    # The calculations are divided into three cases:
    # 1. When `num` is less than `nside + 1`.
    # 2. When `num` is less than `3nside + 1`.
    # 3. When `num` is greater than or equal to `3nside + 1`.
    #
    # The function returns the start and stop indices as integers.
    npix = nside2npix(cp.nside)
    n = num-1
    a1 = 4
    d = 4
    if num < cp.nside + 1
        start = 1/2*n*(2*a1+(n-1)*d)+1
        stop  = 1/2*num*(2*a1+(num-1)*d)
    elseif num < 3cp.nside+1
        n_2 = num - cp.nside
        n = cp.nside - 1
        start = 1/2*cp.nside*(2*a1+(cp.nside-1)*d) + 4*(n_2-1)*cp.nside + 1
        stop = 1/2*cp.nside*(2*a1+(cp.nside-1)*d) + 4*(n_2)*cp.nside 
    else
        n_2 = 4cp.nside-1 - num
        start = npix - 1/2*(n_2+1)*(2*a1+(n_2)*d) +1
        stop = npix - 1/2*(n_2)*(2*a1+(n_2-1)*d)
    end
    return Int(start), Int(stop)
end

function unique_theta_val(cp)
    npix = nside2npix(cp.nside)
    res = Resolution(cp.nside)
    θ= zeros(npix)
    for i in 1:npix
        ang = pix2angRing(res, i)
        θ[i] = ang[1]
    end
    return unique(θ)
end