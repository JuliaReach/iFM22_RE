using ReachabilityAnalysis
using ReachabilityAnalysis: center
using Plots
using LaTeXStrings

# w = 2pi/T
A = [0.0 1; -4π 0]
X0 = BallInf([0.0, 10.0], 0.01)
X = Universe(2)
b = [1.0, -1.0]
B = hcat(b)
u = -3.0
U = ConstantInput(linear_map(B, Singleton([u])))

P = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ X0)

# homogenize system
α = 1.0
A_homog = [A b*u/α; zeros(1, 3)]
X0_homog = cartesian_product(X0, Singleton([α]))
X_homog = Universe(3)
P_homog = @ivp(x' = A_homog*x, x ∈ X_homog, x(0) ∈ X0_homog)
P = P_homog

# flowpipe (tight approximation)
δ = 0.00001
T = 0.05
T_full = 0.15

alg_fw = GLGM06(δ=δ, approx_model=Forward())
sol_fw = solve(P; T=T, alg=alg_fw);

# discretizations
Ω0_chull = ConvexHullArray(map(set, sol_fw))
Ω0_ch = discretize(P, T, CorrectionHull(order=4)).x0

# random trajectories
x0s = vcat(sample(X0, 0, include_vertices=true),
           sample(X0, 50, sampler=LazySets.FaceSampler(1)))

P_orbit = [@ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ x0) for x0 in x0s]
sol_orbit = [solve(Pi, T=T_full, alg=ORBIT(δ=T/20.)) for Pi in P_orbit]

fig = plot(legend=:bottomleft, legendfontsize=15, xticks=[], yticks=[], size=(900, 400))

Φ = exp(A_homog * T)

sol_fw_ch = [Projection(Φ^i * Ω0_ch, [1, 2]) for i in 0:2]
sol_fw_chull = [ConvexHullArray([Projection(Φ^i * X, [1, 2]) for X in array(Ω0_chull)]) for i in 0:2]

[plot!(fig, s, vars=(1, 2), alpha=0.2, c=:magenta, seriestype=:path, marker=:none) for s in sol_orbit]
plot!(fig, sol_fw_ch[2:end], alpha=0, color=:blue, linecolor=:blue, lw=2, linealpha=1)
plot!(fig, sol_fw_chull[2:end], color=:purple, linecolor=:purple, lw=2, linealpha=1)
plot!(fig, sol_fw_ch[1], alpha=0, color=:blue, linecolor=:blue, lw=2, linealpha=1)
plot!(fig, sol_fw_chull[1], 1e-4, alpha=0.3, color=:purple, linecolor=:purple, lw=2, linealpha=1)

plot!(fig, X0, c=:yellow, alpha=1., lab="")
plot!(fig, sol_fw[1501:2000], vars=(1, 2), lw=0.1, c=:orange, alpha=0.8, lab="")
plot!(fig, sol_fw[end], vars=(1, 2), c=:green, alpha=0.8, lab="")

p = center(box_approximation(sol_fw_ch[1]))
plot!(fig, color=:black, annotations=(p[1], p[2] - .12, text(L"\Omega_0", 22)))
p = center(box_approximation(sol_fw_ch[2]))
plot!(fig, color=:black, annotations=(p[1] - .025, p[2] - .1, text(L"\Omega_1", 22)))
p = center(box_approximation(sol_fw_ch[3]))
plot!(fig, color=:black, annotations=(p[1] - .05, p[2] - .07, text(L"\Omega_2", 22)))

savefig(fig, "intro_recurrence.pdf")
