function unique_theta_num(num, NSIDE)
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