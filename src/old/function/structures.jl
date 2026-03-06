mutable struct ConvolutionParams{I<:Int, MC<:Matrix{Complex{Float64}}, VI<:Vector{Int64},Bl<:Bool}
    nside::I
    lmax::I
    alm::MC
    blm::MC
    beam_mmax::I
    l_range::VI
    HWP::Bl
    realization::I
end

function gen_ConvolutionParams(;
        nside=12*128^2,
        lmax=3*nside-1,
        alm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        blm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        beam_mmax = 2,
        l_range = [0,lmax],
        HWP = false,
        realization = 1
    )
    return ConvolutionParams(
        nside,
        lmax,
        alm,
        blm,
        beam_mmax,
        l_range,
        HWP,
        realization
    )
end

mutable struct ConvolutionParams_pc{I<:Int, MC<:Matrix{Complex{Float64}}, VI<:Vector{Int64}, AC<:Array{ComplexF64,2}}
    nside::I
    lmax::I
    alm::MC
    blm::MC
    beam_mmax::I
    l_range::VI
    ini_wignerd::AC
end

function gen_ConvolutionParams_pc(;
        nside=12*128^2,
        lmax=3*nside-1,
        alm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        blm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        beam_mmax = 2,
        l_range = [0,lmax],
        ini_wignerd = zeros(ComplexF64, 3, 2*lmax+1)
    )
    return ConvolutionParams_pc(
        nside,
        lmax,
        alm,
        blm,
        beam_mmax,
        l_range, 
        ini_wignerd
    )
end

mutable struct ConvolutionParams_{I<:Int, MC<:Matrix{Complex{Float64}}, VI<:Vector{Int64}, AC<:Array{Complex{Float64}, 3}}
    nside::I
    lmax::I
    alm::MC
    blm::MC
    l_range::VI
    ini_wignerd::AC
    calc_wignerD::AC
end

function gen_ConvolutionParams_(;
        nside=12*128^2,
        lmax=3*nside-1,
        alm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        blm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        l_range = [0,lmax],
        ini_wignerd = zeros(ComplexF64,l_range[2]-l_range[1]+1, 2l_range[2]+1,2l_range[2]+1),
        calc_wignerD = zeros(ComplexF64,l_range[2]-l_range[1]+1, 2l_range[2]+1,2l_range[2]+1)
    )
    return ConvolutionParams_(
        nside,
        lmax,
        alm,
        blm,
        l_range,
        ini_wignerd,
        calc_wignerD
    )
end