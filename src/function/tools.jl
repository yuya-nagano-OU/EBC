"""
Return the inclusive pixel index range (in RING ordering) for a given ring.

Arguments
- ring::Int : HEALPix ring index, 1 <= ring <= 4nside-1
- nside::Int

Returns
- (start, stop)::Tuple{Int,Int}
"""
function ring_pixel_range(ring::Int, nside::Int)
    nrings = 4nside - 1
    if ring < 1 || ring > nrings
        throw(ArgumentError("ring must be in 1:$nrings for nside=$nside, got $ring"))
    end

    nphi = pixels_in_ring(ring, nside)
    start = first_pixel_in_ring(ring, nside)
    stop  = start + nphi - 1
    return start, stop
end

"""
Number of pixels in a given ring (RING ordering).
"""
function pixels_in_ring(ring::Int, nside::Int)
    if ring <= nside
        return 4ring
    elseif ring <= 3nside
        return 4nside
    else
        return 4(4nside - ring)
    end
end

"""
First pixel index (1-based, RING ordering) of a given ring.
"""
function first_pixel_in_ring(ring::Int, nside::Int)
    if ring <= nside
        # sum_{k=1}^{ring-1} 4k + 1
        return 2ring*(ring - 1) + 1
    elseif ring <= 3nside
        # north cap pixels + equatorial offset + 1
        northcap = 2nside*(nside - 1)
        return northcap + (ring - nside - 1) * 4nside + 1
    else
        # use south-side symmetry
        nrings = 4nside - 1
        r′ = nrings - ring + 1   # mirrored ring counted from south pole
        npix = 12nside^2
        return npix - 2r′*(r′ + 1) + 1
    end
end


function unique_theta_val(nside::Int)
    nrings = 4nside - 1
    θ = Vector{Float64}(undef, nrings)

    for r in 1:nrings
        z = if r <= nside
            1 - r^2 / (3nside^2)
        elseif r <= 3nside
            (4/3) - (2r)/(3nside)
        else
            rr = 4nside - r
            -1 + rr^2 / (3nside^2)
        end
        θ[r] = acos(z)
    end

    return θ
end

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
#=
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
=#
