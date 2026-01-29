function unique_theta_num(num, nside)
    #=
    # This function calculates the unique theta number range for a given pixel number and nside parameter.
    #
    # # Arguments
    # - `num::Int`: The pixel number for which the unique theta number range is to be calculated.
    # - `nside::Int`: The nside parameter, which determines the resolution of the HEALPix grid.
    #
    # # Returns
    # - `Tuple{Int, Int}`: A tuple containing the start and stop indices of the unique theta number range.
    =#
    npix = nside2npix(nside)
    if num < 1 || num > npix
        throw(ArgumentError("num must be in 1:$npix for nside=$nside, got $num"))
    end
    n = num-1
    a1 = 4
    d = 4
    if num < nside + 1
        start = 1/2*n*(2*a1+(n-1)*d)+1
        stop  = 1/2*num*(2*a1+(num-1)*d)
    elseif num < 3nside+1
        n_2 = num - nside
        n = nside - 1
        start = 1/2*nside*(2*a1+(nside-1)*d) + 4*(n_2-1)*nside + 1
        stop = 1/2*nside*(2*a1+(nside-1)*d) + 4*(n_2)*nside
    else
        n_2 = 4nside-1 - num
        start = npix - 1/2*(n_2+1)*(2*a1+(n_2)*d) +1
        stop = npix - 1/2*(n_2)*(2*a1+(n_2-1)*d)
    end
    return Int(start), Int(stop)
end

function unique_theta_val(nside)
    npix = nside2npix(nside)
    res = Resolution(nside)
    θ= zeros(npix)
    for i in 1:npix
        ang = pix2angRing(res, i)
        θ[i] = ang[1]
    end
    return unique(θ)
end
