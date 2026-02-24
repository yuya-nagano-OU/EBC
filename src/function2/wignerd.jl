function wignerd_mdown_from_l(l::Int, n::Int, beta::Real; eps::Real=1e-12)
    (l >= 0) || throw(ArgumentError("l must be non-negative"))
    (abs(n) <= l) || throw(ArgumentError("n must satisfy |n| <= l"))

    dvals = zeros(ComplexF64, l + 1)

    if abs(sin(beta)) <= eps
        for m in 0:l
            dvals[m + 1] = WignerD.wignerDjmn(l, m, n, 0.0, beta, 0.0)
        end
        return dvals
    end

    d_l = WignerD.wignerDjmn(l, l, n, 0.0, beta, 0.0)
    dvals[l + 1] = d_l

    l == 0 && return dvals

    cosβ = cos(beta)
    sinβ = sin(beta)

    d_m_plus1 = 0.0 + 0.0im
    d_m = d_l

    for m in l:-1:1
        am1 = 2 / sqrt((l - (m - 1)) * (l + (m - 1) + 1))
        λm = (n - m * cosβ) / sinβ

        if m == l
            d_m_minus1 = λm * am1 * d_m
        else
            am = 2 / sqrt((l - m) * (l + m + 1))
            d_m_minus1 = λm * am1 * d_m - (am1 / am) * d_m_plus1
        end

        dvals[m] = d_m_minus1
        d_m_plus1 = d_m
        d_m = d_m_minus1
    end

    return dvals
end
