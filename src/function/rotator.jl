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

function local_effective_wignerD_conj_reduced_formapmake(cb, cc, α, β, γ, ψ; τ::Int=5)
    n_beam = sum(2*min(l, cb.mmax) + 1 for l in cc.lstart:cc.lstop)
    n_sky = sum(2*l + 1 for l in cc.lstart:cc.lstop)
    D_temp = 0+0im
    D_beam = spzeros(ComplexF64, n_sky, n_beam)
    D_beam2 = spzeros(ComplexF64, n_sky, n_beam)
    D_beam3 = spzeros(ComplexF64, n_sky, n_beam)
    @inbounds for l in cc.lstart:cc.lstop
        n_ = min(l, cb.mmax)
        @inbounds for i in eachindex(α)
            phase = exp(2im*ψ[i])
            phase2 = exp(4im*ψ[i])
            @inbounds for m in -l:l
                m_idx = lmr_idx(l=l, m=m, lstart=cc.lstart, mmax=cc.lstop)
                # ★帯制限：|m-n| <= τ だけ回す
                n_lo = max(-n_, m - τ)
                n_hi = min( n_, m + τ)
                #@show n_lo, n_hi, n_
                @inbounds for n in n_lo:n_hi
                    sgn = pm(m, n)
                    n_idx = lmr_idx(l=l, m=n, lstart=cc.lstart, mmax=cb.mmax)
                    D_temp = sgn * WignerD.wignerDjmn(l, -m, -n, α[i], β[i], γ[i])
                    D_beam[m_idx, n_idx] += D_temp
                    D_beam2[m_idx, n_idx] += D_temp * phase
                    D_beam3[m_idx, n_idx] += D_temp * conj(phase)
                end
            end
        end
    end
    return D_beam./length(α), D_beam2./length(α), D_beam3./length(α)
end
