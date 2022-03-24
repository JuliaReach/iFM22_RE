using ReachabilityAnalysis
using ReachabilityAnalysis: translate, solve
using Plots
using LaTeXStrings

### deterministic case

# IVP parameters
A = [0.0 1; -4π 0]
b = [1.0, -1.0]
B = hcat(b)
X0 = BallInf([0.0, 10.0], 0.1)
X = Universe(2)
u = -3.0
U_raw = Singleton([u])
U = ConstantInput(linear_map(B, U_raw))
X_homog = Universe(3)

# discretization parameters
δ = 0.1
alg = Forward()

# homogenization parameters
αs = [1.0]

# result for heterogeneous system
P_heterog = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ X0)
Pd_heterog = discretize(P_heterog, δ, alg)
Ω0_heterog = Pd_heterog.x0

# result for homogenized system
Ω0_homog = []
for α in αs
    A_homog = [A b*u/α; zeros(1, 3)]
    X0_homog = cartesian_product(X0, Singleton([α]))
    P_homog = @ivp(x' = A_homog*x, x ∈ X_homog, x(0) ∈ X0_homog)

    Pd_homog = discretize(P_homog, δ, alg)
    local Ω0 = Projection(Pd_homog.x0, [1, 2])
#     Ω0 = overapproximate(Ω0, PolarDirections(30))  # workaround because σ is not available
    push!(Ω0_homog, Ω0)  # project auxiliary variable
end

# plot sets
plot(leg=:bottomleft, legendfontsize=20, tickfontsize=20)
plot!(Ω0_heterog, alpha=0, color=:blue, linecolor=:blue, linealpha=1, lw=2, lab=L"\Omega_0~\textrm{heterog}")
for i in eachindex(αs)
    plot!(Ω0_homog[i], alpha=0.3, color=:green, linecolor=:green, linealpha=1, lw=2, lab=L"\Omega_0~\textrm{homog}")
end
plot!(X0, alpha=1, color="yellow", lab=L"X_0")
# sampling initial states at time t from vertices and edges of X0
# and then plotting the reachable states for some random input signal
x0s = vcat(sample(X0, 0, include_vertices=true),
           sample(X0, 100, sampler=LazySets.FaceSampler(1)))
for x0 in x0s
    local P = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ x0)
    sol = solve(P, δ, ORBIT(δ=δ/10))
    plot!(sol, vars=(1, 2), c=:magenta, seriestype=:path, marker=:none)
end

savefig("homogenize_deterministic.pdf")


### nondeterministic case

U_raw = BallInf([u], 0.2)
U = ConstantInput(linear_map(B, U_raw))

# result for heterogeneous system
P_heterog = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ X0)
Pd_heterog = discretize(P_heterog, δ, alg)
Ω0_heterog = Pd_heterog.x0

# result for homogenized system
Ω0_homog = []
for α in αs
    A_homog = [A b*u/α; zeros(1, 3)]
    X0_homog = cartesian_product(X0, Singleton([α]))
    U_raw_homog = cartesian_product(translate(U_raw, -[u]), Singleton([0.0]))
    B_homog = [B zeros(size(B, 1)); zeros(size(B, 2)) 0]
    U_homog = ConstantInput(linear_map(B_homog, U_raw_homog))
    P_homog = @ivp(x' = A_homog*x + u, x ∈ X_homog, u ∈ U_homog, x(0) ∈ X0_homog)

    Pd_homog = discretize(P_homog, δ, alg)
    push!(Ω0_homog, Projection(Pd_homog.x0, [1, 2]))  # project auxiliary variable
end

# plot sets
plot(leg=:bottomleft, legendfontsize=20, tickfontsize=20)
plot!(Ω0_heterog, alpha=0, color=:blue, linecolor=:blue, linealpha=1, lw=2, lab=L"\Omega_0~\textrm{heterog}")
for i in eachindex(αs)
    plot!(Ω0_homog[i], alpha=0.3, color=:green, linecolor=:green, linealpha=1, lw=2, lab=L"\Omega_0~\textrm{homog}")
end
plot!(X0, alpha=1, color="yellow", lab=L"X_0")
# sampling initial states at time t from vertices and edges of X0
# and then plotting the reachable states for some random input signal
x0s = vcat(sample(X0, 0, include_vertices=true),
           sample(X0, 100, sampler=LazySets.FaceSampler(1)))
for x0 in x0s
    local P = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ x0)
    sol = solve(P, δ, ORBIT(δ=δ/10))
    plot!(sol, vars=(1, 2), c=:magenta, seriestype=:path, marker=:none)
end

savefig("homogenize_nondeterministic.pdf")
