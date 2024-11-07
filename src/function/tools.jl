function unique_theta_num(num, NSIDE)
# This function calculates the start and stop indices for unique theta values
# based on the given pixel number (num) and NSIDE parameter.
#
# Arguments:
# - num: The pixel number for which the unique theta indices are to be calculated.
# - NSIDE: The NSIDE parameter that determines the resolution of the map.
#
# Returns:
# - A tuple containing the start and stop indices as integers.
#
# The function works by first calculating the total number of pixels (NPIX) using
# the nside2npix function. It then determines the start and stop indices based on
# the value of num in relation to NSIDE. The calculations are divided into three
# cases:
# 1. When num is less than NSIDE + 1.
# 2. When num is less than 3 * NSIDE + 1.
# 3. When num is greater than or equal to 3 * NSIDE + 1.
#
# Each case uses a different formula to calculate the start and stop indices.
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
    #map_TQU[1,:] .= beam_map
    for i in 1:npix
        ang = pix2angRing(res, i)
        θ[i] = ang[1]
    end
    return unique(θ)
end
