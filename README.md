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

To reproduce the plots from the paper, run the script `experiments/plot_all.jl`:

```julia
julia> include("experiments/plot_all.jl")
```
