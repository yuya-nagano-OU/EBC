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