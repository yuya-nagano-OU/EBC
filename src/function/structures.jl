mutable struct ConvolutionParams{I<:Int, MC<:Matrix{Complex{Float64}}, VI<:Vector{Int64}}
    nside::I
    lmax::I
    alm::MC
    blm::MC
    l_range::VI
end

function gen_ConvolutionParams(;
        nside=12*128^2,
        lmax=3*nside-1,
        alm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        blm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        l_range = [0,lmax]
    )
    return ConvolutionParams(
        nside,
        lmax,
        alm,
        blm,
        l_range
    )
end