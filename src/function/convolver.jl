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
        #print(l)
        result_1 += transpose(alm_new[1][l+1][:])*conj(eff_wignerD[l+1][:,:]*blm_new[1][l+1][mmax+1-l:mmax+1+l])
        result_2 += transpose(alm_new[2][l+1][:])*conj(eff_wignerD[l+1][:,:]*blm_new[2][l+1][mmax+1-l:mmax+1+l])
        result_3 += transpose(alm_new[3][l+1][:])*conj(eff_wignerD[l+1][:,:]*blm_new[3][l+1][mmax+1-l:mmax+1+l])
    end
    for l in mmax+1:cp.lmax
        #print(l)
        result_1 += transpose(alm_new[1][l+1][:])*conj(eff_wignerD[l+1][:,l+1-mmax:l+1+mmax]*blm_new[1][l+1][:])
        result_2 += transpose(alm_new[2][l+1][:])*conj(eff_wignerD[l+1][:,l+1-mmax:l+1+mmax]*blm_new[2][l+1][:])
        result_3 += transpose(alm_new[3][l+1][:])*conj(eff_wignerD[l+1][:,l+1-mmax:l+1+mmax]*blm_new[3][l+1][:])
    end
    return result_1 + (result_2 + result_3)/2
end