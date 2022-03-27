using ReachabilityAnalysis
using Plots

#
# w = 2pi/T
#
A = [0.0 1; -4π 0]
X0 = BallInf([0.0, 10.0], 0.01)
X = Universe(2)
B = hcat([1.0, -1.0])
U = ConstantInput(linear_map(B, Singleton([-3.0])))

P = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ X0)

# flowpipe (tight approximation)
δ = 0.00001
T = 0.05

alg_fw = GLGM06(δ=δ, approx_model=Forward())
sol_fw = solve(P; T=T, alg=alg_fw);

# discretizations
Ω0_chull = ConvexHullArray(map(set, sol_fw))
Ω0_fw = discretize(P, T, Forward()).x0
Ω0_ch = discretize(P, T, CorrectionHull(order=4)).x0

# random trajectories
x0s = vcat(sample(X0, 0, include_vertices=true),
           sample(X0, 50, sampler=LazySets.FaceSampler(1)))

P_orbit = [@ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ x0) for x0 in x0s]
sol_orbit = [solve(Pi, T=T, alg=ORBIT(δ=T/20.)) for Pi in P_orbit]

fig = plot(legend=:bottomright, legendfontsize=15, xticks=[], yticks=[])

plot!(fig, Ω0_ch, alpha=0, color=:blue, linecolor=:blue, lw=2, linealpha=1)
plot!(fig, Ω0_chull, 1e-4, alpha=0.3, color=:purple, linecolor=:purple, lw=2, linealpha=1)

plot!(fig, X0, c=:yellow, alpha=1., lab="")
plot!(fig, sol_fw[1501:2000], vars=(1, 2), lw=0.1, c=:orange, alpha=0.8, lab="")
plot!(fig, sol_fw[end], vars=(1, 2), c=:green, alpha=0.8, lab="")
[plot!(fig, s, vars=(1, 2), c=:magenta, seriestype=:path, marker=:none) for s in sol_orbit]

savefig(fig, "intro_approximate.pdf")
