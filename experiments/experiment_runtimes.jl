using ReachabilityAnalysis
import DataStructures
using DataStructures: OrderedDict
using JuMP, Ipopt

include("models.jl")

# scaling factor for time unit (base is in seconds)
unit_scale = 1e3  # milliseconds

# methods
algs_homog = [
        (FirstOrderddt(), "d-dt"),
        (FirstOrderZonotope(), "Zonotope"),
        (CorrectionHull(order=4), "Correction-hull"),
        (FirstOrder(), "First-order"),
        (ForwardBackward(solver=Ipopt.Optimizer()), "Forward-backward"),
        (Forward(), "Forward")
       ]
# excludes the "d/dt" method
algs_heterog = [
        (FirstOrderZonotope(), "Zonotope"),
        (CorrectionHull(order=4), "Correction-hull"),
        (FirstOrder(), "First-order"),
        (ForwardBackward(solver=Ipopt.Optimizer()), "Forward-backward"),
        (Forward(), "Forward")
       ]

function get_time(collection, m)
    if m != 1
        # m is the number of runs per experiment
        sort!(collection)
        part = ceil(Int, m / 20)
        collection = collection[part:m-part]
    end
    time = sum(collection) / length(collection) * unit_scale
end

# run experiment
function experiment(P, δs, model_name, algs; m::Int=1)
    d = ones(dim(P.x0))  # direction for comparison
    δ2name2results = OrderedDict{Float64, OrderedDict}()
    for δ in δs
        local name2results = OrderedDict{String, Tuple}()
        for (alg, name) in algs
            collection = Vector{Float64}(undef, m)
            local Ω0
            @inbounds for i in 1:m
                res = @timed discretize(P, δ, alg)
                Ω0 = res.value.x0
                collection[i] = res.time
            end
            local time1 = get_time(collection, m)

            collection = Vector{Float64}(undef, m)
            local support_value
            @inbounds for i in 1:m
                res = @timed ρ(d, Ω0)
                support_value = res.value
                collection[i] = res.time
            end
            local time2 = get_time(collection, m)
            local tuple = (Ω0, support_value, time1 + time2)
            name2results[name] = tuple
        end
        δ2name2results[δ] = name2results
        println("δ = $δ")
    end

    # write results to table
    for (idx, filename) in [(2, "delta_support_$(model_name).dat"),
                            (3, "delta_time_$(model_name).dat")]
        open(filename, "w") do f
            # header
            write(f, "delta")
            for (_, name) in algs
                write(f, "\t$name")
            end
            write(f, "\n")

            # data
            for (δ, name2results) in δ2name2results
                write(f, "$δ")
                for (i, (name, tuple)) in enumerate(name2results)
                    @assert name == algs[i][2] "$name != $(algs[i][2])"
                    data = tuple[idx]
                    write(f, "\t$data")
                end
                write(f, "\n")
            end
        end
    end

    return δ2name2results
end

# delta steps
δs_oscillator = vcat(
          range(0.01, stop=0.06, step=0.01)
         )
δs_freedom = vcat(
          range(0.000005, stop=0.0001, step=0.000005),
          range(0.0002, stop=0.001, step=0.0001),
          range(0.002, stop=0.01, step=0.001)
         )
δs_iss = vcat(
          [0.00002],  # precompilation
          range(0.00002, stop=0.0001, step=0.00001),
          range(0.0002, stop=0.001, step=0.0001),
          range(0.002, stop=0.02, step=0.001)
         )

# run experiments
experiment(oscillator(), δs_oscillator, "oscillator", algs_homog; m=40)
experiment(freedom(), δs_freedom, "freedom", algs_homog; m=40)
experiment(iss(), δs_iss, "iss", algs_heterog)
