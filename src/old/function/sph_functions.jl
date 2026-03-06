#=
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
=#

function get_reorder_alm(alm, lmax)
    new_alm = zeros(ComplexF64, 3, lmax+1, 2lmax+1)
    for l in 0:lmax
        for m in -l:1:-1
            idx = for_healpy_order(l, -m, lmax)
            new_alm[1,l+1,lmax+1+m] = conj(alm[1,idx])*(-1)^m
            new_alm[2,l+1,lmax+1+m] = -(conj(alm[2,idx])+1im*conj(alm[3,idx]))*(-1)^m
            new_alm[3,l+1,lmax+1+m] = -(conj(alm[2,idx])-1im*conj(alm[3,idx]))*(-1)^m
        end
        for m in 0:l
            idx = for_healpy_order(l,m, lmax)
            new_alm[1,l+1,lmax+1+m] = alm[1,idx]
            new_alm[2,l+1,lmax+1+m] = -(alm[2,idx]+1im*alm[3,idx])
            new_alm[3,l+1,lmax+1+m] = -(alm[2,idx]-1im*alm[3,idx])
        end
    end
    return new_alm
end

function get_reorder_alm_optimized(alm, lmax)
    new_alm = [[zeros(ComplexF64, 2*i+1) for i in 0:cp.lmax ],[zeros(ComplexF64,2*i+1) for i in 0:cp.lmax ],[zeros(ComplexF64,2*i+1) for i in 0:cp.lmax ] ]
    for l in 0:lmax
        for m in -l:1:-1
            idx = for_healpy_order(l, -m, lmax)
            new_alm[1][l+1][l+1+m] = conj(alm[1,idx])*(-1)^m
            new_alm[2][l+1][l+1+m] = -(conj(alm[2,idx])+1im*conj(alm[3,idx]))*(-1)^m
            new_alm[3][l+1][l+1+m] = -(conj(alm[2,idx])-1im*conj(alm[3,idx]))*(-1)^m
        end
        for m in 0:l
            idx = for_healpy_order(l,m, lmax)
            new_alm[1][l+1][l+1+m] = alm[1,idx]
            new_alm[2][l+1][l+1+m] = -(alm[2,idx]+1im*alm[3,idx])
            new_alm[3][l+1][l+1+m] = -(alm[2,idx]-1im*alm[3,idx])
        end
    end
    return new_alm
end

function get_reorder_blm_optimized(blm, lmax, mmax)
    new_blm = [[zeros(ComplexF64, 2*mmax+1) for i in 0:cp.lmax ],[zeros(ComplexF64,2*mmax+1) for i in 0:cp.lmax ],[zeros(ComplexF64,2*mmax+1) for i in 0:cp.lmax ] ]
    for l in 0:mmax
        for m in -l:1:-1
            idx = for_healpy_order(l, -m, lmax)
            new_blm[1][l+1][mmax+1+m] = conj(blm[1,idx])*(-1)^m
            new_blm[2][l+1][mmax+1+m] = -(conj(blm[2,idx])+1im*conj(blm[3,idx]))*(-1)^m
            new_blm[3][l+1][mmax+1+m] = -(conj(blm[2,idx])-1im*conj(blm[3,idx]))*(-1)^m
        end
        for m in 0:l
            idx = for_healpy_order(l,m, lmax)
            new_blm[1][l+1][mmax+1+m] = blm[1,idx]
            new_blm[2][l+1][mmax+1+m] = -(blm[2,idx]+1im*blm[3,idx])
            new_blm[3][l+1][mmax+1+m] = -(blm[2,idx]-1im*blm[3,idx])
        end
    end
    
    for l in mmax+1:lmax
        for m in -mmax:1:-1
            idx = for_healpy_order(l, -m, lmax)
            new_blm[1][l+1][mmax+1+m] = conj(blm[1,idx])*(-1)^m
            new_blm[2][l+1][mmax+1+m] = -(conj(blm[2,idx])+1im*conj(blm[3,idx]))*(-1)^m
            new_blm[3][l+1][mmax+1+m] = -(conj(blm[2,idx])-1im*conj(blm[3,idx]))*(-1)^m
        end
        for m in 0:mmax
            idx = for_healpy_order(l,m, lmax)
            new_blm[1][l+1][mmax+1+m] = blm[1,idx]
            new_blm[2][l+1][mmax+1+m] = -(blm[2,idx]+1im*blm[3,idx])
            new_blm[3][l+1][mmax+1+m] = -(blm[2,idx]-1im*blm[3,idx])
        end
    end
    return new_blm
end

function get_reorder_blm_optimized_2(blm, lmax, mmax)
    new_blm = [[zeros(ComplexF64, min(2*i+1, 2*mmax+1)) for i in 0:cp.lmax ],[zeros(ComplexF64, min(2*i+1, 2*mmax+1)) for i in 0:cp.lmax ],[zeros(ComplexF64, min(2*i+1, 2*mmax+1)) for i in 0:cp.lmax ] ]
    for l in 0:mmax
        for m in -l:1:-1
            idx = for_healpy_order(l, -m, lmax)
            new_blm[1][l+1][l+1+m] = conj(blm[1,idx])*(-1)^m
            new_blm[2][l+1][l+1+m] = -(conj(blm[2,idx])+1im*conj(blm[3,idx]))*(-1)^m
            new_blm[3][l+1][l+1+m] = -(conj(blm[2,idx])-1im*conj(blm[3,idx]))*(-1)^m
        end
        for m in 0:l
            idx = for_healpy_order(l,m, lmax)
            new_blm[1][l+1][l+1+m] = blm[1,idx]
            new_blm[2][l+1][l+1+m] = -(blm[2,idx]+1im*blm[3,idx])
            new_blm[3][l+1][l+1+m] = -(blm[2,idx]-1im*blm[3,idx])
        end
    end
    
    for l in mmax+1:lmax
        for m in -mmax:1:-1
            idx = for_healpy_order(l, -m, lmax)
            new_blm[1][l+1][mmax+1+m] = conj(blm[1,idx])*(-1)^m
            new_blm[2][l+1][mmax+1+m] = -(conj(blm[2,idx])+1im*conj(blm[3,idx]))*(-1)^m
            new_blm[3][l+1][mmax+1+m] = -(conj(blm[2,idx])-1im*conj(blm[3,idx]))*(-1)^m
        end
        for m in 0:mmax
            idx = for_healpy_order(l,m, lmax)
            new_blm[1][l+1][mmax+1+m] = blm[1,idx]
            new_blm[2][l+1][mmax+1+m] = -(blm[2,idx]+1im*blm[3,idx])
            new_blm[3][l+1][mmax+1+m] = -(blm[2,idx]-1im*blm[3,idx])
        end
    end
    return new_blm
end