# NiGaussianMixture

Differentiable (reversible) Gaussian Mixture model.

The motivation is to beat the benchmark in this [paper](https://arxiv.org/abs/1807.10129), for the glory of Julia community!

[![Build Status](https://travis-ci.com/JuliaReverse/NiGaussianMixture.jl.svg?branch=master)](https://travis-ci.com/JuliaReverse/NiGaussianMixture.jl)
[![Codecov](https://codecov.io/gh/JuliaReverse/NiGaussianMixture.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaReverse/NiGaussianMixture.jl)

## Get started

Open a Julia REPL and type `]` to enter `pkg` mode and then type
```julia pkg
pkg> dev git@github.com:JuliaReverse/NiGaussianMixture.jl.git
pkg> add ForwardDiff BenchmarkTools
```

Then in a bash shell, type the following commands to open the benchmark file in Atom.
```bash
$ julia ~/.julia/dev/NiGaussianMixture/benchmarks/benchmark.jl
```

You will see results like:
```julia repl
loading data: /home/leo/.julia/dev/NiGaussianMixture/data/gmm/10k/gmm_d10_K5.txt
Normal Objective
BenchmarkTools.Trial:
  memory estimate:  4.43 MiB
  allocs estimate:  69121
  --------------
  minimum time:     8.937 ms (0.00% GC)
  median time:      9.088 ms (0.00% GC)
  mean time:        9.683 ms (5.64% GC)
  maximum time:     15.601 ms (33.67% GC)
  --------------
  samples:          516
  evals/sample:     1
Reversible Objective
BenchmarkTools.Trial:
  memory estimate:  19.44 MiB
  allocs estimate:  160434
  --------------
  minimum time:     49.781 ms (0.00% GC)
  median time:      52.293 ms (0.00% GC)
  mean time:        53.114 ms (3.52% GC)
  maximum time:     65.018 ms (8.92% GC)
  --------------
  samples:          95
  evals/sample:     1

NiLang Gradient
BenchmarkTools.Trial:
  memory estimate:  51.74 MiB
  allocs estimate:  320860
  --------------
  minimum time:     136.134 ms (2.61% GC)
  median time:      138.738 ms (2.92% GC)
  mean time:        148.102 ms (3.16% GC)
  maximum time:     238.323 ms (2.57% GC)
  --------------
  samples:          34
  evals/sample:     1

ForwardDiff Gradient
nparams = 330
BenchmarkTools.Trial:
  memory estimate:  659.90 MiB
  allocs estimate:  1935756
  --------------
  minimum time:     1.291 s (3.62% GC)
  median time:      1.307 s (3.78% GC)
  mean time:        1.306 s (3.75% GC)
  maximum time:     1.319 s (3.76% GC)
  --------------
  samples:          4
  evals/sample:     1
```

Note: this result is not yet fully optimized.
