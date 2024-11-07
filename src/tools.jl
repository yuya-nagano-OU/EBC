function unique_theta_num(num, NSIDE)
    # This function calculates the unique theta number range for a given pixel number and NSIDE parameter.
    #
    # # Arguments
    # - `num::Int`: The pixel number for which the unique theta number range is to be calculated.
    # - `NSIDE::Int`: The NSIDE parameter, which determines the resolution of the HEALPix grid.
    #
    # # Returns
    # - `Tuple{Int, Int}`: A tuple containing the start and stop indices of the unique theta number range.
    #
    # The function first calculates the total number of pixels (NPIX) for the given NSIDE using the `nside2npix` function.
    # It then determines the start and stop indices based on the value of `num` relative to `NSIDE` and `3NSIDE`.
    # The calculations are divided into three cases:
    # 1. When `num` is less than `NSIDE + 1`.
    # 2. When `num` is less than `3NSIDE + 1`.
    # 3. When `num` is greater than or equal to `3NSIDE + 1`.
    #
    # The function returns the start and stop indices as integers.
    NPIX = nside2npix(NSIDE)
    n = num-1
    a1 = 4
    d = 4
    if num < NSIDE + 1
        start = 1/2*n*(2*a1+(n-1)*d)+1
        stop  = 1/2*num*(2*a1+(num-1)*d)
    elseif num < 3NSIDE+1
        n_2 = num - NSIDE
        n = NSIDE - 1
        start = 1/2*NSIDE*(2*a1+(NSIDE-1)*d) + 4*(n_2-1)*NSIDE + 1
        stop = 1/2*NSIDE*(2*a1+(NSIDE-1)*d) + 4*(n_2)*NSIDE 
    else
        n_2 = 4NSIDE-1 - num
        start = NPIX - 1/2*(n_2+1)*(2*a1+(n_2)*d) +1
        stop = NPIX - 1/2*(n_2)*(2*a1+(n_2-1)*d)
    end
    return Int(start), Int(stop)
end

function　unique_theta_val(nside, res)
    npix = nside2npix(nside)
    θ= zeros(npix)
    for i in 1:npix
        ang = pix2angRing(res, i)
        θ[i] = ang[1]
    end
    return unique(θ)
end