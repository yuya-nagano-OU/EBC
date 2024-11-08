function make_order_alm_2(alm, lmax::Integer, l::Integer, mmax::Integer) #almをl固定 m: -l~l までの順番にする
    new_alm = zeros(ComplexF64, 2mmax+1)
    for m in -mmax:mmax
        if m < 0
            idx = for_healpy_order(l, -m, lmax)
            new_alm[m + mmax + 1] = conj(alm[idx])*(-1)^m
        else
            idx = for_healpy_order(l, m, lmax)
            new_alm[m+1 + mmax] = alm[idx]
        end
    end
    return new_alm
end