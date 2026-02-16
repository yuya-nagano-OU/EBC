mutable struct ConvolutionCalculate{I<:Int, Bl<:Bool}
    nside_output::I
    lstart::I
    lstop::I
    mmax_calculate::I
    HWP::Bl
end

function ConvolutionCalculate(;
    nside_output::Int = 128,
    lstart::Int = 0,
    lstop::Int = 3*nside_output - 1,
    mmax_calculate::Int = 2,
    HWP::Bool = false
)
    lstart <= lstop || throw(ArgumentError("lstart must be <= lstop (got lstart=$lstart, lstop=$lstop)"))
    mmax_calculate <= lstop || throw(ArgumentError("mmax_calculate must be <= lstop (got mmax_calculate=$mmax_calculate, lstop=$lstop)"))
    return ConvolutionCalculate(
        nside_output,
        lstart,
        lstop,
        mmax_calculate,
        HWP
    )
end


