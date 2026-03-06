function eff_convolver(alm_new, blm_new, eff_wignerD)
    result_1 = 0.0
    result_2 = 0.0
    result_3 = 0.0
    for l in 0:cp.lmax
        result_1 += transpose(alm_new[1,l+1,lmax+1-l:lmax+1+l])*conj(eff_wignerD[l+1][:,:]*blm_new[1,l+1,lmax+1-l:lmax+1+l])
        result_2 += transpose(alm_new[2,l+1,lmax+1-l:lmax+1+l])*conj(eff_wignerD[l+1][:,:]*blm_new[2,l+1,lmax+1-l:lmax+1+l])
        result_3 += transpose(alm_new[3,l+1,lmax+1-l:lmax+1+l])*conj(eff_wignerD[l+1][:,:]*blm_new[3,l+1,lmax+1-l:lmax+1+l])
    end
    return result_1 + (result_2 + result_3)/2
end

function eff_convolver_optimized(alm_new, blm_new, eff_wignerD, mmax)
    result_1 = 0.0 + 0.0im 
    result_2 = 0.0 + 0.0im
    result_3 = 0.0 + 0.0im
    for l in 0:mmax
        result_1 += transpose(alm_new[1][l+1][:])*conj(eff_wignerD[l+1][:,:]*blm_new[1][l+1][mmax+1-l:mmax+1+l])
        result_2 += transpose(alm_new[2][l+1][:])*conj(eff_wignerD[l+1][:,:]*blm_new[2][l+1][mmax+1-l:mmax+1+l])
        result_3 += transpose(alm_new[3][l+1][:])*conj(eff_wignerD[l+1][:,:]*blm_new[3][l+1][mmax+1-l:mmax+1+l])
    end
    for l in mmax+1:cp.lmax
        result_1 += transpose(alm_new[1][l+1][:])*conj(eff_wignerD[l+1][:,l+1-mmax:l+1+mmax]*blm_new[1][l+1][:])
        result_2 += transpose(alm_new[2][l+1][:])*conj(eff_wignerD[l+1][:,l+1-mmax:l+1+mmax]*blm_new[2][l+1][:])
        result_3 += transpose(alm_new[3][l+1][:])*conj(eff_wignerD[l+1][:,l+1-mmax:l+1+mmax]*blm_new[3][l+1][:])
    end
    return result_1 + (result_2 + result_3)/2
end

function eff_convolver_optimized(alm_new, blm_new, eff_wignerD)
    result_1 = 0.0 + 0.0im 
    result_2 = 0.0 + 0.0im
    result_3 = 0.0 + 0.0im
    for l in 0:cp.lmax
        #print(l)
        result_1 += transpose(alm_new[1][l+1][:])*conj(eff_wignerD[l+1][:,:]*blm_new[1][l+1][:])
        result_2 += transpose(alm_new[2][l+1][:])*conj(eff_wignerD[l+1][:,:]*blm_new[2][l+1][:])
        result_3 += transpose(alm_new[3][l+1][:])*conj(eff_wignerD[l+1][:,:]*blm_new[3][l+1][:])
    end
    return result_1 + (result_2 + result_3)/2
end

function eff_convolver_optimized_2(alm_new, blm_new, eff_wignerD)
    result_1 = 0.0 + 0.0im 
    result_2 = 0.0 + 0.0im
    result_3 = 0.0 + 0.0im
    @views for l in 0:cp.lmax
        #@show l, length(blm_new[1][l+1][:]), length(eff_wignerD[:,:])
        result_1 += transpose(alm_new[1][l+1])*conj.(eff_wignerD[l+1]*blm_new[1][l+1])
        result_2 += transpose(alm_new[2][l+1])*conj.(eff_wignerD[l+1]*blm_new[2][l+1])
        result_3 += transpose(alm_new[3][l+1])*conj.(eff_wignerD[l+1]*blm_new[3][l+1])
    end
    return result_1 .+ (result_2 .+ result_3)./2
end

function eff_convolver_optimized_2(cp,alm_new, blm_new, eff_wignerD)
    result_1 = 0.0 + 0.0im 
    result_2 = 0.0 + 0.0im
    result_3 = 0.0 + 0.0im
    @inbounds for l in 0:cp.lmax
        result_1 += transpose(alm_new[1][l+1])*conj.(eff_wignerD[l+1]*blm_new[1][l+1])
        result_2 += transpose(alm_new[2][l+1])*conj.(eff_wignerD[l+1]*blm_new[2][l+1])
        result_3 += transpose(alm_new[3][l+1])*conj.(eff_wignerD[l+1]*blm_new[3][l+1])
    end
    return result_1 .+ (result_2 .+ result_3)./2
end

function eff_convolver_optimized_22(cp,alm_new, blm_new, eff_wignerD)
    result_1 = 0.0 + 0.0im
    temp_1 = 0.0 + 0.0im
    result_2 = 0.0 + 0.0im
    result_3 = 0.0 + 0.0im
    @inbounds for l in 0:cp.lmax
        @show l
        @time temp_1 = (eff_wignerD[l+1]*blm_new[1][l+1])
        @time result_1 += dot(transpose(alm_new[1][l+1]),conj.(temp_1))
        @time temp_2 = conj.(eff_wignerD[l+1]*blm_new[2][l+1])
        @time result_2 += transpose(alm_new[2][l+1])*temp_2
        @time temp_3 = conj.(eff_wignerD[l+1]*blm_new[3][l+1])
        @time result_3 += transpose(alm_new[3][l+1])*temp_3
    end
    return result_1 .+ (result_2 .+ result_3)./2
end

function eff_convolver_optimized_222(cp,alm_new, blm_new, eff_wignerD)
    result_1 = 0.0 + 0.0im
    temp_1 = 0.0 + 0.0im
    result_2 = 0.0 + 0.0im
    result_3 = 0.0 + 0.0im
    @inbounds for l in 0:cp.lmax
        temp_1 = transpose(alm_new[1][l+1])*conj.(eff_wignerD[l+1])
        result_1 += temp_1*conj.(blm_new[1][l+1])
        temp_2 = transpose(alm_new[2][l+1])*conj.(eff_wignerD[l+1])
        result_2 += temp_2*conj.(blm_new[2][l+1])
        temp_3 = transpose(alm_new[3][l+1])*conj.(eff_wignerD[l+1])
        result_3 += temp_3*conj.(blm_new[3][l+1])
    end
    return result_1 .+ (result_2 .+ result_3)./2
end

function eff_convolver_optimized_fastest(cp,alm_new, blm_new, eff_wignerD)
    result_1 = 0.0 + 0.0im
    temp_1 = 0.0 + 0.0im
    result_2 = 0.0 + 0.0im
    result_3 = 0.0 + 0.0im
    @inbounds for l in 0:cp.lmax
        temp_1 = transpose(alm_new[1][l+1])*conj.(eff_wignerD[l+1])
        result_1 += temp_1*conj.(blm_new[1][l+1])
        temp_2 = transpose(alm_new[2][l+1])*conj.(eff_wignerD[l+1])
        result_2 += temp_2*conj.(blm_new[2][l+1])
        temp_3 = transpose(alm_new[3][l+1])*conj.(eff_wignerD[l+1])
        result_3 += temp_3*conj.(blm_new[3][l+1])
    end
    return result_1 .+ (result_2 .+ result_3)./2
end


function eff_convolver_optimized_threading(alm_new, blm_new, eff_wignerD)
    result_1 = zeros(Complex{Float64}, nthreads())
    result_2 = zeros(Complex{Float64}, nthreads())
    result_3 = zeros(Complex{Float64}, nthreads())
    
    @views Threads.@threads for l in 0:cp.lmax
        tid = threadid()
        result_1[tid] += transpose(alm_new[1][l+1]) * conj.(eff_wignerD[l+1] * blm_new[1][l+1])
        result_2[tid] += transpose(alm_new[2][l+1]) * conj.(eff_wignerD[l+1] * blm_new[2][l+1])
        result_3[tid] += transpose(alm_new[3][l+1]) * conj.(eff_wignerD[l+1] * blm_new[3][l+1])
    end
    
    final_result_1 = sum(result_1)
    final_result_2 = sum(result_2)
    final_result_3 = sum(result_3)
    
    return final_result_1 .+ (final_result_2 .+ final_result_3) ./ 2
end

#=
function eff_convolver_optimized_lv(alm_new, blm_new, eff_wignerD)
    result_1 = 0.0 + 0.0im 
    result_2 = 0.0 + 0.0im
    result_3 = 0.0 + 0.0im
    @views for l in 0:cp.lmax
        #@show l, length(blm_new[1][l+1][:]), length(eff_wignerD[:,:])
        @turbo result_1 += transpose(alm_new[1][l+1])*conj.(eff_wignerD[l+1]*blm_new[1][l+1])
        @turbo result_2 += transpose(alm_new[2][l+1])*conj.(eff_wignerD[l+1]*blm_new[2][l+1])
        @turbo result_3 += transpose(alm_new[3][l+1])*conj.(eff_wignerD[l+1]*blm_new[3][l+1])
    end
    return result_1 .+ (result_2 .+ result_3)./2
end
=#

function eff_convolver_optimized_transpose(alm_new, blm_new, eff_wignerD)
    result_1 = 0.0 + 0.0im 
    result_2 = 0.0 + 0.0im
    result_3 = 0.0 + 0.0im
    @views for l in 0:cp.lmax
        #@show l, length(blm_new[1][l+1][:]), length(eff_wignerD[:,:])
        result_1 += dot(alm_new[1][l+1],conj.(eff_wignerD[l+1]*blm_new[1][l+1]))
        result_2 += dot(alm_new[2][l+1],conj.(eff_wignerD[l+1]*blm_new[2][l+1]))
        result_3 += dot(alm_new[3][l+1],conj.(eff_wignerD[l+1]*blm_new[3][l+1]))
    end
    return result_1 .+ (result_2 .+ result_3)./2
end