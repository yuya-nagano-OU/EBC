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

function convolver_1pixel(cs, cb, cc, alm_slice, blm_slice, globalD, localD)
    if cc.HWP == false
        result = zeros(ComplexF64, cs.numofsky, cb.numofbeams, 3)
        for i in 1:cs.numofsky
            for j in 1:cb.numofbeams
                result[i,j,1] = transpose(alm_slice[1,1,:])*globalD*localD[1]*blm_slice[:,1,1] + 1/2*transpose(alm_slice[1,2,:])*globalD*localD[1]*blm_slice[:,2,1] + 1/2* transpose(alm_slice[1,3,:])*globalD*localD[1]*blm_slice[:,3,1]
                result[i,j,2] = conj(transpose(alm_slice[1,1,:])*globalD*localD[2]*blm_slice[:,1,1]) + 1/2*transpose(alm_slice[1,2,:])*globalD*localD[1]*blm_slice[:,2,1] + 1/2*conj(transpose(alm_slice[1,3,:])*globalD*localD[3]*blm_slice[:,3,1])
                result[i,j,3] = conj(result[i,j,2])
            end
        end
    end
    return result
end