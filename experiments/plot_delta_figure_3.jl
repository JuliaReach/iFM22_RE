using ReachabilityAnalysis
using JuMP, Ipopt
using Plots
using LaTeXStrings

# IVP parameters
A = [0.0 1; -4π 0]
X0 = BallInf([0.0, 10.0], 0.1)
X = Universe(2)
P = @ivp(x' = A*x, x ∈ X, x(0) ∈ X0)

# varying time steps
δs = [0.02, 0.01, 0.006, 0.003]

# methods
algs = [
        (FirstOrderddt(), "d/dt"),
        (FirstOrderZonotope(), "Zonotope"),
        (CorrectionHull(order=4), "Correction hull"),
        (FirstOrder(), "First order"),
        (ForwardBackward(solver=Ipopt.Optimizer()), "Forward/backward"),
        (Forward(), "Forward")
       ]

# reachable states at time t
R(t) = linear_map(exp(A * t), X0)

# results
δ2Ω0s = []
for δ in δs
    local Ω0s = []
    for (alg, name) in algs
        local Pd = discretize(P, δ, alg)
        local Ω0 = Pd.x0
        push!(Ω0s, Ω0)
    end
    push!(δ2Ω0s, Ω0s)
end

# plot sets
colors = [:orange, :green, :blue, :purple, :red, :black]
plots = []
for (i, δ) in enumerate(δs)
    local fig = plot(leg=:bottomleft, xlab=L"\delta = %$δ")
    for (j, (alg, name)) in enumerate(algs)
        local Ω0 = δ2Ω0s[i][j]
        if name == "Forward/backward"
            Ω0 = overapproximate(Ω0, PolarDirections(30))  # workaround because σ is not available
        end
        plot!(Ω0, alpha=0, color=colors[j], linecolor=colors[j], linealpha=1, lw=2, lab="Ω₀ $(name)")
    end

    [plot!(R(δi), color=:gray, alpha=0.05, linealpha=0, lab="") for δi in range(0.0; stop=δ, length=20)]

    savefig("experiment_delta_$i.pdf")
    push!(plots, fig)
end

return plots
