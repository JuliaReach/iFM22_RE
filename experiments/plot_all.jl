# script to create all plots

for file in ("plot_intro_exact.jl",
             "plot_intro_approximate.jl",
             "plot_underapproximation.jl",
             "plot_homogenize.jl",
             "plot_shrink_delta.jl",
             "plot_combine.jl",
             "experiment_delta.jl",
             "experiment_X0.jl",
             "experiment_U.jl")
    println("Running $file ...")
    include(file)
end

nothing
