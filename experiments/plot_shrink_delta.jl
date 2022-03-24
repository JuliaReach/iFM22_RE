using ReachabilityAnalysis
using Plots
using LaTeXStrings

# IVP parameters
A = [0.0 1; -4π 0]
X0 = BallInf([0.0, 10.0], 0.1)
X = Universe(2)
P = @ivp(x' = A*x, x ∈ X, x(0) ∈ X0)

# discretization parameters
δ = 0.1
alg = Forward()
k = 10

# result for original time step
Pd = discretize(P, δ, alg)
Ω0 = Pd.x0

# result for shrunk time step
δ_small = δ / k
Pd_small = discretize(P, δ_small, alg)
Ωs = ConvexHullArray()
push!(array(Ωs), Pd_small.x0)
Φ_small_pow = Pd_small.s.A
for i in 1:k-1
    push!(array(Ωs), Φ_small_pow * Pd_small.x0)
    global Φ_small_pow *= Pd_small.s.A
end
Ω0_small = concretize(Ωs)

# plot sets
plot(leg=:bottomleft, legendfontsize=20, tickfontsize=20)
plot!(Ω0, alpha=0, color=:blue, linecolor=:blue, linealpha=1, lw=2, lab=L"\Omega_0\ \delta = %$δ")
plot!(Ω0_small, alpha=0.3, color=:green, linecolor=:green, linealpha=1, lw=2, lab=L"\Omega_0\ \gamma = %$δ_small")
for i in eachindex(array(Ωs))
    plot!(Ωs[i], color=:gray, alpha=0.2, lab="")
end
plot!(X0, alpha=1, color="yellow", lab=L"X_0")
x0s = vcat(sample(X0, 0, include_vertices=true),
           sample(X0, 30, sampler=LazySets.FaceSampler(1)),
           sample(X0, 30))
for x0 in x0s
    local P = @ivp(x' = A*x, x ∈ X, x(0) ∈ x0)
    sol = solve(P, δ, ORBIT(δ=δ/10))
    plot!(sol, vars=(1, 2), c=:magenta, seriestype=:path, marker=:none)
end

savefig("shrink_delta.pdf")
