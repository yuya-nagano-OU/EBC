mutable struct ConvolutionCalculate{I<:Int, Bl<:Bool}
    nside_output::I
    resol
    lstart::I
    lstop::I
    mmax_calculate::I
    HWP::Bl
end

function ConvolutionCalculate(;
    nside_output::Int = 128,
    resol = Resolution(nside_output),
    lstart::Int = 0,
    lstop::Int = 3*nside_output - 1,
    mmax_calculate::Int = 2,
    HWP::Bool = false
)
    lstart <= lstop || throw(ArgumentError("lstart must be <= lstop (got lstart=$lstart, lstop=$lstop)"))
    mmax_calculate <= lstop || throw(ArgumentError("mmax_calculate must be <= lstop (got mmax_calculate=$mmax_calculate, lstop=$lstop)"))
    return ConvolutionCalculate(
        nside_output,
        resol,
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
                result[i,j,1] = transpose(alm_slice[i,1,:])*globalD*localD[1]*blm_slice[:,1,j] + 1/2*transpose(alm_slice[i,2,:])*globalD*localD[1]*blm_slice[:,2,j] + 1/2* transpose(alm_slice[i,3,:])*globalD*localD[1]*blm_slice[:,3,j]
                result[i,j,2] = (transpose(alm_slice[i,1,:])*globalD*localD[2]*blm_slice[:,1,j]) + 1/2*transpose(alm_slice[i,2,:])*globalD*localD[2]*blm_slice[:,2,j] + 1/2*(transpose(alm_slice[i,3,:])*globalD*localD[2]*blm_slice[:,3,j])
                result[i,j,3] = (transpose(alm_slice[i,1,:])*globalD*localD[3]*blm_slice[:,1,j]) + (1/2*transpose(alm_slice[i,2,:])*globalD*localD[3]*blm_slice[:,2,j]) + 1/2*(transpose(alm_slice[i,3,:])*globalD*localD[3]*blm_slice[:,3,j])
                #result[i,j,3] = conj(result[i,j,2])
            end
        end
    end
    return result
end

function compute_pixel_convolution(pix_idx, globalD, φ, θ, ψ; τ=5)
    α_local, β_local, γ_local = calc_local_euiler_angles(cc.resol, pix_idx, φ, θ, ψ)
    localD = local_effective_wignerD_conj_reduced_formapmake(cb, cc, α_local, β_local, γ_local, ψ, τ=τ)
    result = convolver_1pixel(cs, cb, cc, alm_slice, blm_slice, globalD, localD)
    return result
end