using ReachabilityAnalysis
using ReachabilityAnalysis: solve
using JuMP, Ipopt
using Plots
using LaTeXStrings

# IVP parameters
A = [0.0 1; -4π 0]
X0 = BallInf([0.0, 10.0], 0.1)
X = Universe(2)

# discretization parameters
δ = 0.01
b = [1.0, -1.0]
B = hcat(b)
u1 = -3.0

# varying input domains
Us_raw = [
          Singleton([100 * u1]),
          Singleton([u1]),
          BallInf([0.0], abs(u1)),
          BallInf([u1], 1.0)
         ]

Us = [
      (ConstantInput(linear_map(B, Us_raw[1])), "{$(Us_raw[1].element[1])}"),
      (ConstantInput(linear_map(B, Us_raw[2])), "{$u1}"),
      (ConstantInput(linear_map(B, Us_raw[3])), "□(0, $(abs(u1)))"),
      (ConstantInput(linear_map(B, Us_raw[4])), "□($u1, 0.1)")
     ]

# methods
algs = [
        (FirstOrderZonotope(), "Zonotope"),
        (CorrectionHull(order=4), "Correction hull"),
        (FirstOrder(), "First order"),
        (ForwardBackward(solver=Ipopt.Optimizer()), "Forward/backward"),
        (Forward(), "Forward")
       ]

# results
U2Ω0s = []
for (U, name_U) in Us
    local Ω0s = []
    local P = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ X0)
    for (alg, name_alg) in algs
        local Pd = discretize(P, δ, alg)
        local Ω0 = Pd.x0
        push!(Ω0s, Ω0)
    end
    push!(U2Ω0s, Ω0s)
end

# plot sets
colors = [:orange, :green, :blue, :purple, :red, :black]
plots = []
for (i, (U_raw, (U, name_U))) in enumerate(zip(Us_raw, Us))
    local fig = plot(leg=:bottomleft, xlab="U = B·$name_U")
    for (j, (alg, name_alg)) in enumerate(algs)
        local Ω0 = U2Ω0s[i][j]
        if name_alg == "Forward/backward"
            Ω0 = overapproximate(Ω0, PolarDirections(30))  # workaround because σ is not available
        end
        # j+1 because d/dt is missing
        plot!(Ω0, alpha=0, color=colors[j+1], linecolor=colors[j+1], linealpha=1, lw=2, lab="Ω₀ $(name_alg)")
    end

    local P = @ivp(x' = A*x + u, x ∈ X, u ∈ U, x(0) ∈ X0)
    sol = solve(P, δ, alg=LGG09(δ=δ/10, dirs=collect(PolarDirections(40))))
    alpha = i == 1 ? 0.2 : 0.05  # first plot should use a darker gray
    plot!(sol, vars=(1, 2), color=:gray, alpha=alpha, linealpha=0, lab="")

    savefig("experiment_U_$i.pdf")
    push!(plots, fig)
end

return plots
