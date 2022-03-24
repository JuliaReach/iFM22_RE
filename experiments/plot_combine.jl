using ReachabilityAnalysis
using ReachabilityAnalysis: solve
using JuMP, Ipopt
using Plots
using LaTeXStrings

# IVP parameters
A = [0.0 1; -4π 0]
b = [1.0, -1.0]
B = hcat(b)
X0 = BallInf([0.0, 10.0], 10.3)
X = Universe(2)
u = -3.0
U_raw = BallInf([u], 10.0)
U = ConstantInput(linear_map(B, U_raw))
P = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ X0)

# discretization parameters
δ = 0.05

# methods
algs = [
        (FirstOrderZonotope(), "Zonotope"),
        (CorrectionHull(order=4), "Correction hull"),
        (FirstOrder(), "First order"),
        (ForwardBackward(solver=Ipopt.Optimizer()), "Forward/backward"),
        (Forward(), "Forward")
       ]

# results
Ω0s = LazySet{Float64}[]
for (alg, name) in algs
    local Pd = discretize(P, δ, alg)
    local Ω0 = Pd.x0
    if name == "Forward/backward"
        Ω0 = overapproximate(Ω0, PolarDirections(30))  # workaround because σ is not available
    end
    push!(Ω0s, Ω0)
end

# intersect results
Ω0_intersect = concretize(IntersectionArray(Ω0s))

# plot sets
plot(leg=:bottomleft)
for (Ω0, (alg, name)) in zip(Ω0s, algs)
    plot!(Ω0, lab="Ω₀ $(name)")
end
plot!(Ω0_intersect, lab="⋂ Ω₀")
plot!(X0, alpha=1, color="yellow", lab=L"X_0")

# sampling initial states at time t from vertices and edges of X0
# and then plotting the reachable states for some random input signal
x0s = vcat(sample(X0, 0, include_vertices=true),
           sample(X0, 200, sampler=LazySets.FaceSampler(1)),
           sample(X0, 400))
for x0 in x0s
    local P = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ x0)
    sol = solve(P, δ, ORBIT(δ=δ/10))
    plot!(sol, vars=(1, 2), c=:magenta, seriestype=:path, marker=:none)
end

savefig("intersection.pdf")
