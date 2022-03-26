# iFM22_RE

Repeatibility evaluation for "Conservative Time Discretization: A Comparative Study" (iFM'21) by Marcelo Forets and Christian Schilling.

The article was accepted at the [17th International Conference on integrated Formal Methods](https://ifm22.si.usi.ch/), which took place in Lugano, Switzerland between June 7-10, 2022.

An arXiv version is available [here](https://arxiv.org/abs/2111.01454).

## 💾 Installation

All you need is a Julia compiler ([available here](https://julialang.org/downloads/)).
At the time of writing, we used version v1.7.2.
Newer versions may or may not work, and older versions are [available here](https://julialang.org/downloads/oldreleases/).

Download the content of this repository and, from the main folder, start the Julia REPL and activate the project:

```shell
$ julia --project=.
```

Install all recorded dependencies:

```julia
julia> using Pkg; Pkg.instantiate()
```
To exit the Julia REPL, use `exit()`.

To run this repeatability evaluation, it is recommended that you create a new `output` directory to store the generated results and start Julia from that folder:

```shell
$ mkdir output
$ cd output
$ julia --project=..
```

For reference, this repository includes the folder `results/` containing the files we obtained in the experiments for the paper.

## ☑️ Complete evaluation

To run the whole evaluation, run the script `experiments/run_all.jl`.

```julia
julia> include("../experiments/run_all.jl")
```

### 🖼️ Visual evaluation

To reproduce only the plots from the paper, run the script `experiments/plot_all.jl`:

```julia
julia> include("../experiments/plot_all.jl")
```

### ⚙️ Quantitative evaluation

To reproduce only the quantitative experiments from the paper, run the script `experiments/experiment_runtimes.jl`:

```julia
julia> include("../experiments/experiment_runtimes.jl")
```

The results will be written to `.dat` files. These can be used to create the tables and plots in the paper.

TODO describe

### :notebook_with_decorative_cover: Further documentation

For documentation specific to the set representation library used in this repeatability evaluation, see [LazySets.jl](https://github.com/JuliaReach/LazySets.jl#lazysetsjl). For documentation specific to the systems modeling language, see [MathematicalSystems.jl](https://github.com/JuliaReach/MathematicalSystems.jl#mathematicalsystemsjl). Finally, the actual implementation of the conservative time discretization methods used can be found in the library [ReachabilityAnalysis.jl](https://github.com/JuliaReach/ReachabilityAnalysis.jl#reachabilityanalysisjl).
