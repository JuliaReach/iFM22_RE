include("heat3d_model.jl")

# computes out <- exp(A * NSTEPS * dt) * b
function _expv(A, b, NSTEPS, dt; hermitian=false, m=min(30, size(A, 1)), tol=1e-7)

    # initialization of the krylov subspace
    TA, Tb = eltype(A), eltype(b)
    T = promote_type(TA, Tb)
    Ks = KrylovSubspace{T, real(T)}(length(b), m)
    arnoldi!(Ks, A, b; m=m, ishermitian=hermitian, tol=tol)

    out = similar(b)
    expv!(out, NSTEPS*dt, Ks)

    return out
end

function _phiv(A, b, NSTEPS, dt; hermitian=false, m=min(30, size(A, 1)), tol=1e-7)

    # initialization of the krylov subspace
    TA, Tb = eltype(A), eltype(b)
    T = promote_type(TA, Tb)
    Ks = KrylovSubspace{T, real(T)}(length(b), m)
    arnoldi!(Ks, A, b; m=m, ishermitian=hermitian, tol=tol)

    out = Matrix{Float64}(undef, size(A, 1), 3)
    phiv!(out, NSTEPS*dt, Ks, 2)

    return view(out, :, 3) .* dt^2
end

function _support_function_Phi_dir(ℓ, A, X0, δ; m=30, tol=1e-7, hermitian=false)
    At = transpose(A) 
    vec = _expv(At, ℓ, 1, δ, m=m, tol=tol, hermitian=hermitian)
    ρ(vec, X0)
end

function _support_function_Eplus(ℓ::AbstractVector{N}, A::AbstractMatrix{N}, X0::Hyperrectangle, δ::Float64; cutoff=eps(N), m=min(30, size(A, 1)), tol=1e-7, hermitian=false) where {N}
    c = center(X0)
    r = radius_hyperrectangle(X0)
    A2 = A^2
    rin = abs.(A2 * c) + abs.(A2) * r
    rin = Vector(rin)
    Aabs = abs.(A)
    Pv = _phiv(Aabs, rin, 1, δ; m=m, tol=tol, hermitian=hermitian)
    return dot(abs.(ℓ), Pv)
end

# Non-Krylov version
function forward_non_krylov(instance=heat01; δ=0.02)
    A, Aᵀδ, X0, ℓ = instance(δ=δ)
    prob = @ivp(x' = A*x, x(0) ∈ X0)
    D = discretize(normalize(prob), δ, Forward())
    Ω0 = initial_state(D)
    res = ρ(ℓ, Ω0)
end

# Krylov version of Forward:
# ρ(ℓ, CH(X0, Φ*X0 ⊕ E⨥)) = ρ(ℓ, CH(X0, Φ*X0 ⊕ E⨥)) = max(ρ(ℓ, X0), ρ(Φ^T * ℓ, X0) + ρ(ℓ, E⨥))
function forward_krylov(instance=heat01; δ=0.02, m=30, tol=1e-7, hermitian=false)
    A, Aᵀδ, X0, ℓ = instance(δ=δ)
    aux1 = ρ(ℓ, X0)
    aux2 = _support_function_Phi_dir(ℓ, A, X0, δ, m=m, tol=tol, hermitian=hermitian)
    aux3 = _support_function_Eplus(ℓ, A, X0, δ, m=m, tol=tol, hermitian=hermitian)
    res = max(aux1, aux2 + aux3)
end

# Let Ω0 = CH(X0, Φ*X0 ⊕ E⨥), then we are interested in ρ(ℓ, Φ^k Ω0) for k = 0, 1, ..., NSTEPS-1
function forward_krylov_propagation(instance=heat01; NSTEPS=1, δ=0.02, m=30, tol=1e-7, hermitian=false)
    A, Aᵀδ, X0, ℓ = instance(δ=δ)
    At = copy(transpose(A))

    # initialization of the krylov subspace
    TA, Tb = eltype(At), eltype(ℓ)
    T = promote_type(TA, Tb)
    Ks = KrylovSubspace{T, real(T)}(length(ℓ), m)
    arnoldi!(Ks, At, ℓ; m=m, ishermitian=hermitian, tol=tol)

    ℓ′ = copy(ℓ)
    out = Vector{Float64}()
    for k in 0:NSTEPS-1
        if k > 0
            # ℓ′ = _expv(At, ℓ, k, δ; hermitian=false, m=min(30, size(A, 1)), tol=1e-7)

            # more efficient: use always the same precomputed Krylov subspace
            expv!(ℓ′, k*δ, Ks)
        end
        aux1 = ρ(ℓ′, X0)
        aux2 = _support_function_Phi_dir(ℓ′, A, X0, δ, m=m, tol=tol, hermitian=hermitian)
        aux3 = _support_function_Eplus(ℓ′, A, X0, δ, m=m, tol=tol, hermitian=hermitian)
        res = max(aux1, aux2 + aux3)
        push!(out, res)
    end
    out
end

function run_benchmarks()
    filename = "experiment_krylov_table_1.dat"
    open(filename, "w") do f
        forward_non_krylov(heat01) # warm-up runs
        forward_krylov(heat01)

        # ------------------------------------------------------------------------------
        write(f, "Case Heat01: 5x5x5 grid\n")
        t = @elapsed begin
                res = forward_non_krylov(heat01)
            end
        t = round(t, digits=3)
        write(f, "    Forward non-krylov: runtime: $t, result: $res\n")

        t = @elapsed begin
                res = forward_krylov(heat01, m=30, tol=1e-10, hermitian=true)
            end
        t = round(t, digits=3)
        write(f, "    Forward krylov: runtime: $t, result: $res\n")

        # ------------------------------------------------------------------------------
        write(f, "Case Heat02: 10x10x10 grid\n")
        t = @elapsed begin
                res = forward_non_krylov(heat02)
            end
        t = round(t, digits=3)
        write(f, "    Forward non-krylov: runtime: $t, result: $res\n")

        t = @elapsed begin
                res = forward_krylov(heat02, m=100, tol=1e-10, hermitian=true)
            end
        t = round(t, digits=3)
        write(f,"    Forward krylov: runtime: $t, result: $res\n")

        # ------------------------------------------------------------------------------
        write(f, "Case Heat03: 20x20x20 grid\n")
        write(f, "    Forward non-krylov: runtime: N/A, result: N/A\n")

        t = @elapsed begin
                res = forward_krylov(heat03, m=100, tol=1e-10, hermitian=true)
            end
        t = round(t, digits=3)
        write(f, "    Forward krylov: runtime: $t, result: $res\n")

        # ------------------------------------------------------------------------------
        write(f, "Case Heat04: 50x50x50 grid\n")
        write(f, "    Forward non-krylov: runtime: N/A, result: N/A\n")

        t = @elapsed begin
                res = forward_krylov(heat04, m=100, tol=1e-10, hermitian=true)
            end
        t = round(t, digits=3)
        write(f, "    Forward krylov: runtime: $t, result: $res\n")
    end
end

println("Running Heat3D model ...")
run_benchmarks()

nothing