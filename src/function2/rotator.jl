# --- basic rotation matrices (active, right-handed) ---
Rz(α) = @SMatrix [cos(α) -sin(α) 0.0;
                  sin(α)  cos(α) 0.0;
                  0.0     0.0    1.0]

Ry(β) = @SMatrix [ cos(β) 0.0 sin(β);
                   0.0    1.0 0.0;
                  -sin(β) 0.0 cos(β)]


# ZYZ Euler rotation: R(α,β,γ) = Rz(α) * Ry(β) * Rz(γ)
Rzyz(α, β, γ) = Rz(α) * Ry(β) * Rz(γ)

# angle wrap to (-π, π]
wrap_pm_pi(x) = mod(x + π, 2π) - π

function euler_zyz_from_R(R::AbstractMatrix{<:Real}; eps::Real=1e-12)
    #=
    Extract ZYZ Euler angles (α,β,γ) from a 3×3 rotation matrix R
    Input:
        R   : 3x3 rotation matrix (real-valued, right-handed, active)
        eps : threshold for detecting the singular cases β≈0 or β≈π
    Return:
        (α, β, γ) : ZYZ Euler angles in radians, with α and γ wrapped to (-π, π]
    =#
    # Clamp for numerical safety
    cβ = clamp(R[3,3], -1.0, 1.0)
    β = acos(cβ)

    if abs(sin(β)) > eps
        # Standard case
        α = atan(R[2,3], R[1,3])          # atan2(y, x) in Julia is atan(y, x)
        γ = atan(R[3,2], -R[3,1])
        return (wrap_pm_pi(α), β, wrap_pm_pi(γ))
    else
        # Singular: β ≈ 0 or π
        # If β≈0: R ≈ Rz(α+γ)
        # If β≈π: R ≈ Rz(α-γ) * Ry(π)  (ambiguity remains)
        # We fix α=0 and solve for γ from the top-left 2×2 block.
        if cβ > 0  # β≈0
            α = 0.0
            γ = atan(R[2,1], R[1,1])      # angle of Rz(γ)
            return (α, 0.0, wrap_pm_pi(γ))
        else       # β≈π
            α = 0.0
            # For β=π, top-left 2×2 block ~ [-cosγ -sinγ; -sinγ cosγ] up to convention
            # A stable choice:
            γ = atan(-R[2,1], -R[1,1])
            return (α, π, wrap_pm_pi(γ))
        end
    end
end


function second_euler_for_split(phi, theta, dphi, dtheta, psi; eps=1e-12)
    #=
    Given target Euler angles (phi+dphi, theta+dtheta, psi), split as
    R_tgt = R(phi,theta,0) * R2, and return Euler angles of R2: (α2, β2, γ2).
    Returns:
        (α2, β2, γ2, R2).
    =#
    R1   = Rzyz(phi, theta, 0.0)
    Rtgt = Rzyz(phi + dphi, theta + dtheta, psi)

    # For rotation matrices: inverse = transpose
    R2 = transpose(R1) * Rtgt

    α2, β2, γ2 = euler_zyz_from_R(Matrix(R2); eps=eps)
    return (α2, β2, γ2, Matrix(R2))
end

# --- optional: quick verification helper ---
function check_split(phi, theta, dphi, dtheta, psi; eps=1e-12)
    α2, β2, γ2, R2 = second_euler_for_split(phi, theta, dphi, dtheta, psi; eps=eps)
    R1   = Matrix(Rzyz(phi, theta, 0.0))
    Rtgt = Matrix(Rzyz(phi + dphi, theta + dtheta, psi))
    Rrec = R1 * Matrix(Rzyz(α2, β2, γ2))
    return maximum(abs.(Rtgt .- Rrec)), (α2, β2, γ2)
end

function second_euler(phi, theta, dphi, dtheta, psi; eps=1e-12)
    R1   = Rzyz(phi, theta, 0.0)
    Rtgt = Rzyz(phi + dphi, theta + dtheta, psi)

    # For rotation matrices: inverse = transpose
    R2 = transpose(R1) * Rtgt

    α2, β2, γ2 = euler_zyz_from_R(Matrix(R2); eps=eps)
    return (α2, β2, γ2)
end



function WignerD_calculator!(result::AbstractMatrix{ComplexF64},
                             d, ell::Int, phi, theta, psi;
                             tmpB=nothing, w=nothing, p=nothing, q=nothing)

    L = 2*ell + 1
    @assert size(d,1) == L && size(d,2) == L
    @assert size(result,1) == L && size(result,2) == L

    # 作業配列（外から渡せば allocations ほぼゼロ）
    tmpB === nothing && (tmpB = Matrix{ComplexF64}(undef, L, L))
    w    === nothing && (w    = Vector{ComplexF64}(undef, L))
    p    === nothing && (p    = Vector{ComplexF64}(undef, L))
    q    === nothing && (q    = Vector{ComplexF64}(undef, L))

    zero = ell + 1
    ϕ = phi - pi/2
    ψ = psi + pi/2

    # w[M] = exp(-i*M*theta) を exp ではなく cis で（実数角なら速い）
    @inbounds for i in 1:L
        M = i - zero
        w[i] = cis(-M * theta)
    end

    # tmpB = Diag(w) * d  （＝各行 i を w[i] 倍）
    @inbounds for j in 1:L, i in 1:L
        tmpB[i,j] = d[i,j] * w[i]
    end

    # result = transpose(d) * tmpB
    # 注意: あなたの式は共役を取っていないので adjoint(d) ではなく transpose(d)
    mul!(result, transpose(d), tmpB)

    # p[m], q[n]
    @inbounds for i in 1:L
        m = i - zero
        p[i] = cis(-m * ϕ)
        q[i] = cis(-m * ψ)
    end

    # result[m,n] *= p[m] * q[n]
    @inbounds for j in 1:L, i in 1:L
        result[i,j] *= p[i] * q[j]
    end

    return result
end

function _WignerD_calculator_band!(result::AbstractMatrix{ComplexF64},
                                   d, ell::Int, phi, theta, psi, nmax::Int;
                                   tmpB=nothing, w=nothing, p=nothing, q=nothing)
    L = 2*ell + 1
    @assert size(d,1) == L && size(d,2) == L
    @assert size(result,1) == L && size(result,2) == L
    @assert 0 <= nmax <= ell

    # Work arrays
    tmpB === nothing && (tmpB = Matrix{ComplexF64}(undef, L, L))
    w    === nothing && (w    = Vector{ComplexF64}(undef, L))
    p    === nothing && (p    = Vector{ComplexF64}(undef, L))
    q    === nothing && (q    = Vector{ComplexF64}(undef, L))

    zero = ell + 1
    ϕ = phi - pi/2
    ψ = psi + pi/2

    @inbounds for i in 1:L
        M = i - zero
        w[i] = cis(-M * theta)
    end
    @inbounds for j in 1:L, i in 1:L
        tmpB[i,j] = d[i,j] * w[i]
    end
    @inbounds for i in 1:L
        m = i - zero
        p[i] = cis(-m * ϕ)
        q[i] = cis(-m * ψ)
    end

    fill!(result, 0.0 + 0.0im)
    i0 = zero - nmax
    i1 = zero + nmax
    @inbounds for j in i0:i1, i in i0:i1
        acc = 0.0 + 0.0im
        for k in 1:L
            acc += d[k,i] * tmpB[k,j]
        end
        result[i,j] = acc * p[i] * q[j]
    end

    return result
end

# 複数角度サンプルの単純平均:
#   D_eff = (1/N) * Σ_i D(phi[i], theta[i], psi[i])
function WignerD_calculator!(result::AbstractMatrix{ComplexF64},
                             d, ell::Int,
                             phi::AbstractVector, theta::AbstractVector, psi::AbstractVector;
                             tmpD=nothing,
                             tmpB=nothing, w=nothing, p=nothing, q=nothing)
    N = length(phi)
    @assert length(theta) == N && length(psi) == N

    L = 2*ell + 1
    @assert size(result,1) == L && size(result,2) == L
    tmpD === nothing && (tmpD = Matrix{ComplexF64}(undef, L, L))

    fill!(result, 0.0 + 0.0im)
    @inbounds for i in eachindex(phi)
        WignerD_calculator!(tmpD, d, ell, phi[i], theta[i], psi[i];
                            tmpB=tmpB, w=w, p=p, q=q)
        result .+= tmpD
    end
    result ./= N
    return result
end

# 2step 平均（effective_wignerD 相当）:
#   D_eff = (1/N) * Σ_i D(phi_pix, theta_pix, 0) * D(α_i, β_i, γ_i)
#   where (α_i, β_i, γ_i) = check_split(phi_pix, theta_pix, dphi_i, dtheta_i, psi_i)
function WignerD_calculator!(result::AbstractMatrix{ComplexF64},
                             d, ell::Int,
                             phi::AbstractVector, theta::AbstractVector, psi::AbstractVector,
                             phi_pix, theta_pix;
                             eps=1e-12,
                             nmax=ell,
                             D1=nothing, D2=nothing, tmpMul=nothing,
                             tmpB1=nothing, w1=nothing, p1=nothing, q1=nothing,
                             tmpB2=nothing, w2=nothing, p2=nothing, q2=nothing)
    N = length(phi)
    @assert length(theta) == N && length(psi) == N

    L = 2*ell + 1
    @assert size(result,1) == L && size(result,2) == L
    D1     === nothing && (D1     = Matrix{ComplexF64}(undef, L, L))
    D2     === nothing && (D2     = Matrix{ComplexF64}(undef, L, L))
    tmpMul === nothing && (tmpMul = Matrix{ComplexF64}(undef, L, L))

    # 1st step: pixel center fixed rotation
    WignerD_calculator!(D1, d, ell, phi_pix, theta_pix, 0.0;
                        tmpB=tmpB1, w=w1, p=p1, q=q1)

    fill!(result, 0.0 + 0.0im)
    @inbounds for i in eachindex(phi)
        dphi = phi[i] - phi_pix
        dtheta = theta[i] - theta_pix
        _, (α, β, γ) = check_split(phi_pix, theta_pix, dphi, dtheta, psi[i]; eps=eps)

        _WignerD_calculator_band!(D2, d, ell, α, β, γ, nmax;
                                  tmpB=tmpB2, w=w2, p=p2, q=q2)
        mul!(tmpMul, D1, D2)
        result .+= tmpMul
    end

    result ./= N
    return result
end

# 互換：元と同じ返り値（ただし内部で result を確保）
function WignerD_calculator_fast(d, ell, phi, theta, psi)
    L = 2*ell + 1
    result = Matrix{ComplexF64}(undef, L, L)
    WignerD_calculator!(result, d, ell, phi, theta, psi)
end


function local_effective_wignerD(cb, cc, α, β, γ)
    n_beam = sum(2*min(l, cb.mmax) + 1 for l in cc.lstart:cc.lstop)
    n_sky = sum(2*l + 1 for l in cc.lstart:cc.lstop)
    D_beam = spzeros(ComplexF64, n_sky, n_beam)
    @inbounds for l in cc.lstart:cc.lstop
        n_ = min(l, cb.mmax)
        @inbounds for i in eachindex(α)
            @inbounds for m in -l:l
                m_idx = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.lstop)
                @inbounds for n in -n_:n_
                    n_idx = lmr_idx(l=l, m=n, lstart=cc.lstart, mmax=cb.mmax)
                    D_beam[m_idx, n_idx] += WignerD.wignerDjmn(l, m, n, α[i], β[i], γ[i])
                end
            end
        end
    end
    return D_beam./length(α)
end

function global_wignerD(cc, φ, θ, ψ)
    n_sky = sum(2*l + 1 for l in cc.lstart:cc.lstop)
    D_beam = spzeros(ComplexF64, n_sky, n_sky)
    for l in cc.lstart:cc.lstop
        m_ = l
        @inbounds for m in -m_:m_
            m_idx = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.lstop)
            @inbounds for n in -m_:m_
                n_idx = lmr_idx(l=l, m=n, lstart=cc.lstart, mmax=cc.lstop)
                D_beam[m_idx, n_idx] = WignerD.wignerDjmn(l, m, n, φ, θ, ψ)
            end
        end
    end
    return D_beam
end

# d_beam の l-ブロック( m,n in [-min(l,mmax), min(l,mmax)] )だけを埋める
function fill_wignerD_band_block!(d_beam::AbstractMatrix{ComplexF64},
                                  l::Int, phi, theta, psi;
                                  lstart::Int, mmax::Int)
    mcap = min(l, mmax)
    ms = collect(-mcap:mcap)
    idxs = [lmr_idx(l=l, m=m, lstart=lstart, mmax=mmax) for m in ms]

    @inbounds for (ii, m) in enumerate(ms)
        irow = idxs[ii]
        for (jj, n) in enumerate(ms)
            icol = idxs[jj]
            d_beam[irow, icol] = WignerD.wignerDjmn(l, m, n, phi, theta, psi)
        end
    end
    return d_beam
end


# 同じ l に対して複数角度サンプルの平均を d_beam の l-ブロックに書き込む
function fill_wignerD_band_block_avg!(d_beam::AbstractMatrix{ComplexF64},
                                      l::Int,
                                      phis::AbstractVector, thetas::AbstractVector, psis::AbstractVector;
                                      lstart::Int, mmax::Int)
    N = length(phis)
    @assert length(thetas) == N && length(psis) == N

    mcap = min(l, mmax)
    ms = collect(-mcap:mcap)
    idxs = [lmr_idx(l=l, m=m, lstart=lstart, mmax=mmax) for m in ms]

    # 対象ブロックをゼロ初期化
    @inbounds for i in idxs, j in idxs
        d_beam[i, j] = 0.0 + 0.0im
    end

    @inbounds for k in eachindex(phis)
        ϕ = phis[k]
        θ = thetas[k]
        ψ = psis[k]
        for (ii, m) in enumerate(ms)
            irow = idxs[ii]
            for (jj, n) in enumerate(ms)
                icol = idxs[jj]
                d_beam[irow, icol] += WignerD.wignerDjmn(l, m, n, ϕ, θ, ψ)
            end
        end
    end

    invN = 1.0 / N
    @inbounds for i in idxs, j in idxs
        d_beam[i, j] *= invN
    end

    return d_beam
end


@inline pm(m::Int, n::Int) = isodd(m - n) ? -1.0 : 1.0

function local_effective_wignerD_conj(cb, cc, α, β, γ)
    n_beam = sum(2*min(l, cb.mmax) + 1 for l in cc.lstart:cc.lstop)
    n_sky = sum(2*l + 1 for l in cc.lstart:cc.lstop)
    D_beam = spzeros(ComplexF64, n_sky, n_beam)
    @inbounds for l in cc.lstart:cc.lstop
        n_ = min(l, cb.mmax)
        @inbounds for i in eachindex(α)
            @inbounds for m in -l:l
                m_idx = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.lstop)
                @inbounds for n in -n_:n_
                    sgn = pm(m, n)
                    n_idx = lmr_idx(l=l, m=n, lstart=cc.lstart, mmax=cb.mmax)
                    D_beam[m_idx, n_idx] += sgn * WignerD.wignerDjmn(l, -m, -n, α[i], β[i], γ[i])
                end
            end
        end
    end
    return D_beam./length(α)
end

function local_effective_wignerD_conj_reduced(cb, cc, α, β, γ; τ::Int=5)
    n_beam = sum(2*min(l, cb.mmax) + 1 for l in cc.lstart:cc.lstop)
    n_sky = sum(2*l + 1 for l in cc.lstart:cc.lstop)
    D_beam = spzeros(ComplexF64, n_sky, n_beam)
    @inbounds for l in cc.lstart:cc.lstop
        n_ = min(l, cb.mmax)
        @inbounds for i in eachindex(α)
            @inbounds for m in -l:l
                m_idx = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.lstop)
                # ★帯制限：|m-n| <= τ だけ回す
                n_lo = max(-n_, m - τ)
                n_hi = min( n_, m + τ)
                #@show n_lo, n_hi, n_
                @inbounds for n in n_lo:n_hi
                    sgn = pm(m, n)
                    n_idx = lmr_idx(l=l, m=n, lstart=cc.lstart, mmax=cb.mmax)
                    D_beam[m_idx, n_idx] += sgn * WignerD.wignerDjmn(l, -m, -n, α[i], β[i], γ[i])
                end
            end
        end
    end
    return D_beam./length(α)
end


function global_wignerD_conj(cc, φ, θ, ψ)
    n_sky = sum(2*l + 1 for l in cc.lstart:cc.lstop)
    D_beam = spzeros(ComplexF64, n_sky, n_sky)
    for l in cc.lstart:cc.lstop
        m_ = l
        @inbounds for m in -m_:m_
            m_idx = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.lstop)
            @inbounds for n in -m_:m_
                sgn = pm(m, n)
                n_idx = lmr_idx(l=l, m=n, lstart=cc.lstart, mmax=cc.lstop)
                D_beam[m_idx, n_idx] = sgn * WignerD.wignerDjmn(l, -m, -n, φ, θ, ψ)
            end
        end
    end
    return D_beam
end


function calc_local_euiler_angles(res, pix_idx, φ, θ, ψ)
    alphas = zeros(size(θ))
    betas = zeros(size(θ))
    gammas = zeros(size(θ))
    θ_pix,φ_pix = Healpix.pix2angRing(res, pix_idx)
    dθ = θ .- θ_pix
    dφ = φ .- φ_pix
    for i in eachindex(θ)
        err, (alphas[i], betas[i], gammas[i]) = check_split(φ_pix, θ_pix, dφ[i], dθ[i], ψ[i])
    end
    return alphas, betas, gammas
end