function simple_rotater(lmax, θ, φ, ψ, Blm)
    blm_temp = zeros(ComplexF64, length(Blm))
    for l in 0:lmax
        Blm_l = make_order_alm_2(Blm, lmax, l, l)
        for m in 0:l
            index_result =  for_healpy_order(l, m, lmax)
            for n in -l:l
                W = WignerD.wignerDjmn(l, m, n, φ, θ, ψ)
                blm_temp[index_result] += Blm_l[n.+l.+1]*W
            end
        end
    end
    return blm_temp
end

function effective_wignerD_onz(cp, ψs, initial_wignerd_onz)
    for m in -cp.lmax:cp.lmax
        initial_wignerd_onz[1,cp.lmax+1+m] = mean(exp.(-1im*ψs*m))
        initial_wignerd_onz[2,cp.lmax+1+m] = mean(exp.(-1im*(ψs*(m+2))))
        initial_wignerd_onz[3,cp.lmax+1+m] = mean(exp.(-1im*(ψs*(m-2))))
    end
    return initial_wignerd_onz
end

function effective_beam_on_zaxis(cp, ψs, initial_wignerd_onz)
    for m in -cp.beam_mmax:cp.beam_mmax
        initial_wignerd_onz[1,cp.beam_mmax+1+m] = mean(exp.(-1im*ψs*m))
        initial_wignerd_onz[2,cp.beam_mmax+1+m] = mean(exp.(-1im*(ψs*(m+2))))
        initial_wignerd_onz[3,cp.beam_mmax+1+m] = mean(exp.(-1im*(ψs*(m-2))))
    end
    return initial_wignerd_onz
end

function pc_effective_wignerD_(cp, ψs, initial_wignerd_pc)
    for m in -cp.lmax:cp.lmax
        initial_wignerd_pc[1,cp.lmax+1+m] = mean(exp.(-1im*ψs*m))
        initial_wignerd_pc[2,cp.lmax+1+m] = mean(exp.(-1im*(ψs*(m+2))))
        initial_wignerd_pc[3,cp.lmax+1+m] = mean(exp.(-1im*(ψs*(m-2))))
    end
    return initial_wignerd_pc
end

function initialwignerds_dict(cp, theta, path, initial_wignerd)
    for l in cp.l_range[1]:cp.l_range[2]
        initial_wignerd[l] = swignerd_calc(l, theta, path)
    end
   return initial_wignerd
end

function initialwignerds_array(cp, theta, path, initial_wignerd)
    count = 1
    for l in cp.l_range[1]:cp.l_range[2]
        initial_wignerd[count] = swignerd_calc(l, theta, path)
        count += 1
    end
   return initial_wignerd
end

function get_pc_total_effective_wignerD(cp, phi,  wd_onz, initial_wignerd)
    for l in 0:cp.lmax
        for m in -l:l
            initial_wignerd[l+1][l+1+m,:] .*= exp(-1im*m*(phi))
            initial_wignerd[l+1][l+1+m,:] .*= (wd_onz[cp.lmax+1-l:cp.lmax+1+l])
        end
    end
    return initial_wignerd
end

function get_pc_total_effective_wignerD(cp, phi,  wd_onz, initial_wignerd, convolve_wignerd)
    for l in 0:cp.lmax
        for m in -l:l
            convolve_wignerd[l+1][l+1+m,:] .= exp(-1im*m*(phi)).*initial_wignerd[l+1][l+1+m,:]
            convolve_wignerd[l+1][l+1+m,:] .= (wd_onz[cp.lmax+1-l:cp.lmax+1+l]).*convolve_wignerd[l+1][l+1+m,:]
        end
    end
    return convolve_wignerd
end

function get_pc_total_effective_wignerD_mmaxver(cp, phi,  wd_onz, initial_wignerd, convolve_wignerd)
    for l in 0:cp.beam_mmax
        for m in -l:l
            convolve_wignerd[l+1][l+1+m,:] .= exp(-1im*m*(phi)).*initial_wignerd[l+1][l+1+m,:]
            convolve_wignerd[l+1][l+1+m,:] .= (wd_onz[cp.beam_mmax+1-l:cp.beam_mmax+1+l]).*convolve_wignerd[l+1][l+1+m,:]
        end
    end
    for l in cp.beam_mmax+1:cp.lmax
        for m in -l:l
            convolve_wignerd[l+1][l+1+m,:] .= exp(-1im*m*(phi)).*initial_wignerd[l+1][l+1+m,:]
            convolve_wignerd[l+1][l+1+m,:] .= (wd_onz[:]).*convolve_wignerd[l+1][l+1+m,:]
        end
    end
    return convolve_wignerd
end

function swignerd_calc(ell, theta, path)
    result = zeros(ComplexF64, 2*ell+1, 2*ell+1)
    #d = WignerD.wignerd(ell,pi/2)
    d = npzread(path*"dmatrices=ell$ell"*".npy")
    d_left = zeros(ComplexF64, 2*ell+1, 2*ell+1)
    d_right = zeros(ComplexF64, 2*ell+1, 2*ell+1)
    d_left += transpose(d)
    d_right += d
    zero_position = ell+1
    for m in -ell:ell
        d_left[zero_position+m,:] .*= exp(-im*m*(-pi/2))
        d_left[:,zero_position+m] .*= exp(-im*m*(theta))
        d_right[:,zero_position+m] .*= exp(-im*m*(pi/2))
    end
    result = d_left*d_right
   return result 
    d=d_left=d_right=0
end

function swignerd_calc(ell, theta)
    result = zeros(ComplexF64, 2*ell+1, 2*ell+1)
    d = WignerD.wignerd(ell,pi/2)
    #d = npzread(path*"dmatrices=ell$ell"*".npy")
    d_left = zeros(ComplexF64, 2*ell+1, 2*ell+1)
    d_right = zeros(ComplexF64, 2*ell+1, 2*ell+1)
    d_left += transpose(d)
    d_right += d
    zero_position = ell+1
    for m in -ell:ell
        d_left[zero_position+m,:] .*= exp(-im*m*(-pi/2))
        d_left[:,zero_position+m] .*= exp(-im*m*(theta))
        d_right[:,zero_position+m] .*= exp(-im*m*(pi/2))
    end
    result = d_left*d_right
   return result 
    d=d_left=d_right=0
end

function swignerd_calc(cp, ell, theta, path, d_left, d_right)
    #d = WignerD.wignerd(ell,pi/2)
    zero_position = cp.lmax+1
    zero_position_d = ell+1
    zero_position_r = min(cp.beam_mmax,ell)+1
    t = min(ell, cp.beam_mmax)
    d = npzread(path*"dmatrices=ell$ell"*".npy")
    d_left[zero_position-ell:zero_position+ell,zero_position-ell:zero_position+ell] = transpose(d)
    d_right[zero_position-ell:zero_position+ell,zero_position-t:zero_position+t] = d[:,zero_position_d-t:zero_position_d+t]
    for m in t+1:ell
        d_left[zero_position+m,zero_position-ell:zero_position+ell] .*= exp(-im*m*(-pi/2))
        d_left[zero_position-ell:zero_position+ell,zero_position+m] .*= exp(-im*m*(theta))
        d_left[zero_position-m,zero_position-ell:zero_position+ell] .*= exp(im*m*(-pi/2))
        d_left[zero_position-ell:zero_position+ell,zero_position-m] .*= exp(im*m*(theta))
    end
    for m in -t:t
        #@show m
        d_left[zero_position+m,zero_position-ell:zero_position+ell] .*= exp(-im*m*(-pi/2))
        d_left[zero_position-ell:zero_position+ell,zero_position+m] .*= exp(-im*m*(theta))
        d_right[zero_position-ell:zero_position+ell,zero_position+m] .*= exp(-im*m*(pi/2))
    end
    result = d_left[zero_position-ell:zero_position+ell,zero_position-ell:zero_position+ell]*d_right[zero_position-ell:zero_position+ell,zero_position-t:zero_position+t]
   return result 
    #d_left.*=0
    #d_right.*=0
end

function swignerd_calc_h5(cp, ell, theta, fp, d_left, d_right)
    #d = WignerD.wignerd(ell,pi/2)
    zero_position = cp.lmax+1
    zero_position_d = ell+1
    zero_position_r = min(cp.beam_mmax,ell)+1
    t = min(ell, cp.beam_mmax)
    d = read(fp,"ell$ell")
    d_left[zero_position-ell:zero_position+ell,zero_position-ell:zero_position+ell] = transpose(d)
    d_right[zero_position-ell:zero_position+ell,zero_position-t:zero_position+t] = d[:,zero_position_d-t:zero_position_d+t]
    for m in t+1:ell
        d_left[zero_position+m,zero_position-ell:zero_position+ell] .*= exp(-im*m*(-pi/2))
        d_left[zero_position-ell:zero_position+ell,zero_position+m] .*= exp(-im*m*(theta))
        d_left[zero_position-m,zero_position-ell:zero_position+ell] .*= exp(im*m*(-pi/2))
        d_left[zero_position-ell:zero_position+ell,zero_position-m] .*= exp(im*m*(theta))
    end
    for m in -t:t
        #@show m
        d_left[zero_position+m,zero_position-ell:zero_position+ell] .*= exp(-im*m*(-pi/2))
        d_left[zero_position-ell:zero_position+ell,zero_position+m] .*= exp(-im*m*(theta))
        d_right[zero_position-ell:zero_position+ell,zero_position+m] .*= exp(-im*m*(pi/2))
    end
    result = d_left[zero_position-ell:zero_position+ell,zero_position-ell:zero_position+ell]*d_right[zero_position-ell:zero_position+ell,zero_position-t:zero_position+t]
   return result 
    #d_left.*=0
    #d_right.*=0
end

function wignerd2pix_array_h5(cp, theta, fp, wignerd2pix)
    count = 1
    d_left = zeros(ComplexF64, 2*cp.lmax+1, 2*cp.lmax+1)
    d_right = zeros(ComplexF64,2*cp.lmax+1, 2*cp.lmax+1)
    for l in cp.l_range[1]:cp.l_range[2]
        wignerd2pix[count] = swignerd_calc_h5(cp,l, theta, fp, d_left, d_right)
        count += 1
    end
   #return wignerd2pix
    d_left = d_right = 0
end

function effective_beam_on_zaxis(cp, ψs, initial_wignerd_onz)
    for m in -cp.beam_mmax:cp.beam_mmax
        initial_wignerd_onz[1,cp.beam_mmax+1+m] = mean(exp.(-1im*ψs*m))
        initial_wignerd_onz[2,cp.beam_mmax+1+m] = mean(exp.(-1im*(ψs*(m+2))))
        initial_wignerd_onz[3,cp.beam_mmax+1+m] = mean(exp.(-1im*(ψs*(m-2))))
    end
    return initial_wignerd_onz
end

function get_pc_total_effective_wignerD_mmaxver(cp, phi,  wd_onz, initial_wignerd, convolve_wignerd)
    for l in 0:cp.beam_mmax
        for m in -l:l
            convolve_wignerd[l+1][l+1+m,:] .= exp(-1im*m*(phi)).*initial_wignerd[l+1][l+1+m,:]
            convolve_wignerd[l+1][l+1+m,:] .= (wd_onz[cp.beam_mmax+1-l:cp.beam_mmax+1+l]).*convolve_wignerd[l+1][l+1+m,:]
        end
    end
    for l in cp.beam_mmax+1:cp.lmax
        for m in -l:l
            convolve_wignerd[l+1][l+1+m,:] .= exp(-1im*m*(phi)).*initial_wignerd[l+1][l+1+m,:]
            convolve_wignerd[l+1][l+1+m,:] .= (wd_onz[:]).*convolve_wignerd[l+1][l+1+m,:]
        end
    end
    return convolve_wignerd
end