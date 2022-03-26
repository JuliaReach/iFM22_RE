# script to create all plots

for file in ("plot_intro_exact_figure_1_a.jl",
             "plot_intro_approximate_figure_1_b.jl",
             "plot_intro_recurrence_figure_1_c.jl",
             "plot_homogenize_figure_2a_and_figure_8.jl",
             "plot_shrink_delta_figure_2b.jl",
             "plot_delta_figure_3.jl",
             "plot_X0_figure_4.jl",
             "plot_U_figure_5.jl",
             # figure 6 is generated from the experiment runtime data
             "plot_underapproximation_figure_7.jl")
    println("Running $file ...")
    include(file)
end

nothing
