using ReachabilityAnalysis
using Plots
using LaTeXStrings

# IVP parameters
A = [0.0 1; -4π 0]
b = [1.0, -1.0]
B = hcat(b)
X0 = BallInf([0.0, 10.0], 10.0)
X = Universe(2)
P = @ivp(x' = A*x, x ∈ X, x(0) ∈ X0)

# discretization parameters
δ = 0.03

# reachable states at time t
R(t) = linear_map(exp(A * t), X0)

# methods
alg = FirstOrderddt(oa=false)

# results
Pd = discretize(P, δ, alg)
Ω0 = Pd.x0

# plot sets
plot(leg=:left)
[plot!(R(δi), color=:gray, alpha=0.03, lab="") for δi in range(0.0; stop=δ, length=20)]
plot!(X0, alpha=0, color=:blue, linecolor=:blue, linealpha=1, linewidth=2, lab=L"X_0")
plot!(R(δ), alpha=0, color=:green, linecolor=:green, linealpha=1, linewidth=2, lab=L"R_\delta")
plot!(Ω0, alpha=0, color=:red, linecolor=:red, linealpha=1, linewidth=2, lab=L"\Omega_0")
x0s = vcat(sample(X0, 0, include_vertices=true),
           sample(X0, 200, sampler=LazySets.FaceSampler(1)),
           sample(X0, 400))
for x0 in x0s
    local P = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ x0)
    sol = solve(P, δ, ORBIT(δ=δ/10))
    plot!(sol, vars=(1, 2), c=:magenta, seriestype=:path, marker=:none)
end

savefig("underapproximation.pdf")
