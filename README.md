# iFM22_RE

Repeatibility evaluation for "Conservative Time Discretization: A Comparative Study" (iFM'21) by Marcelo Forets and Christian Schilling.

The article was accepted at the [17th International Conference on integrated Formal Methods](https://ifm22.si.usi.ch/), which took place in Lugano, Switzerland between June 7-10, 2022.

An arXiv version is available [here](https://arxiv.org/abs/2111.01454).

## 💾 Installation

All you need is

1. A Julia compiler ([available here](https://julialang.org/downloads/))
   * At the time of writing, we used version v1.7.2.
   * We assume you have integrated `julia` in your `PATH` variable (otherwise replace the `julia` command below accordingly).
3. The content of this repository.

## ➡️ Start

Move to the main folder of this repository.
To run the repeatability evaluation, first create a new `output` directory to store the generated results and start the Julia REPL from that folder:

```shell
$ mkdir output
$ cd output
$ julia --project=..
```

Once you are in the Julia REPL, first install all required dependencies:

```julia
julia> using Pkg; Pkg.instantiate()
```

Next run the experiments as described in the next section.

To exit the Julia REPL, use `exit()`.

For reference, this repository includes the `results` directory containing the files we obtained in the experiments for the paper.

## ☑️ Experiments

Below we assume that Julia was run from the `output` directory.

To run the whole evaluation, run the script `experiments/run_all.jl`.

```julia
julia> include("../experiments/run_all.jl")
```

This will create the plots in most of the figures as `.pdf` files and store the experimental results in `.dat` files.
The latter can be used to create Table 2 and Figure 6 in the paper as described below.
The data used to create Table 1 is contained in the file `experiment_krylov_table_1.dat`.

### ⚙️ Creating Table 2 and Figure 6

We provide the LaTeX script to create Table 2 and Figure 6 in the `latex` directory.
The script (`table_2_and_figure_6.tex`) requires that the experimental `.dat` files have been generated in the `output` directory.

## 📔 Further documentation

For documentation specific to the set representation library used in this repeatability evaluation, see [LazySets.jl](https://github.com/JuliaReach/LazySets.jl#lazysetsjl). For documentation specific to the systems modeling language, see [MathematicalSystems.jl](https://github.com/JuliaReach/MathematicalSystems.jl#mathematicalsystemsjl). Finally, the actual implementation of the conservative time discretization methods can be found in the library [ReachabilityAnalysis.jl](https://github.com/JuliaReach/ReachabilityAnalysis.jl#reachabilityanalysisjl).
