using ReachabilityAnalysis
using ReachabilityAnalysis: solve
using JuMP, Ipopt
using Plots
using LaTeXStrings

# IVP parameters
A = [0.0 1; -4π 0]
X = Universe(2)

# discretization parameter
δ = 0.005

# varying initial states
X0s = [
       (Singleton([0.0, 10.0]), "□(x₀, 0)"),
       (BallInf([0.0, 10.0], 0.01), "□(x₀, 0.01)"),
       (BallInf([0.0, 10.0], 0.1), "□(x₀, 0.1)"),
       (BallInf([0.0, 10.0], 1.0), "□(x₀, 1)")
      ]

# methods
algs = [
        (FirstOrderddt(), "d/dt"),
        (FirstOrderZonotope(), "Zonotope"),
        (CorrectionHull(order=4), "Correction hull"),
        (FirstOrder(), "First order"),
        (ForwardBackward(solver=Ipopt.Optimizer()), "Forward/backward"),
        (Forward(), "Forward")
       ]

# results
X02Ω0s = []
for (X0, name_X0) in X0s
    local Ω0s = []
    local P = @ivp(x' = A*x, x ∈ X, x(0) ∈ X0)
    for (alg, name_alg) in algs
        local Pd = discretize(P, δ, alg)
        local Ω0 = Pd.x0
        push!(Ω0s, Ω0)
    end
    push!(X02Ω0s, Ω0s)
end

# plot sets
colors = [:orange, :green, :blue, :purple, :red, :black]
plots = []
for (i, (X0, name_X0)) in enumerate(X0s)
    local fig = plot(leg=:bottomleft, xlab="X₀ = $name_X0")
    for (j, (alg, name_alg)) in enumerate(algs)
        local Ω0 = X02Ω0s[i][j]
        if name_alg == "Forward/backward"
            Ω0 = overapproximate(Ω0, PolarDirections(30))  # workaround because σ is not available
        end
        plot!(Ω0, alpha=0, color=colors[j], linecolor=colors[j], linealpha=1, lw=2, lab="Ω₀ $(name_alg)")
    end

    # reachable states at time t
    R(t) = linear_map(exp(A * t), X0)

    [plot!(R(δi), color=:gray, alpha=0.05, linealpha=0, lab="") for δi in range(0.0; stop=δ, length=20)]

    savefig("experiment_X0_$i.pdf")
    push!(plots, fig)
end

return plots
