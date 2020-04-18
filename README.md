# NiGaussianMixture

Differentiable (reversible) Gaussian Mixture model.

The motivation is to beat the benchmark in this [paper](https://arxiv.org/abs/1807.10129), for the glory of Julia community!

[![Build Status](https://travis-ci.com/JuliaReverse/NiGaussianMixture.jl.svg?branch=master)](https://travis-ci.com/JuliaReverse/NiGaussianMixture.jl)
[![Codecov](https://codecov.io/gh/JuliaReverse/NiGaussianMixture.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaReverse/NiGaussianMixture.jl)

## Get started!

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
  minimum time:     5.799 ms (0.00% GC)
  median time:      5.849 ms (0.00% GC)
  mean time:        6.063 ms (2.97% GC)
  maximum time:     12.997 ms (50.32% GC)
  --------------
  samples:          825
  evals/sample:     1
Reversible Objective
BenchmarkTools.Trial: 
  memory estimate:  6.91 MiB
  allocs estimate:  80334
  --------------
  minimum time:     29.817 ms (0.00% GC)
  median time:      30.005 ms (0.00% GC)
  mean time:        30.473 ms (1.29% GC)
  maximum time:     33.439 ms (8.34% GC)
  --------------
  samples:          165
  evals/sample:     1
NiLang Gradient
BenchmarkTools.Trial: 
  memory estimate:  19.04 MiB
  allocs estimate:  160662
  --------------
  minimum time:     79.517 ms (0.00% GC)
  median time:      80.665 ms (0.00% GC)
  mean time:        81.452 ms (1.28% GC)
  maximum time:     91.425 ms (9.19% GC)
  --------------
  samples:          62
  evals/sample:     1
ForwardDiff Gradient
nparams = 330
BenchmarkTools.Trial: 
  memory estimate:  659.90 MiB
  allocs estimate:  1935756
  --------------
  minimum time:     1.086 s (2.56% GC)
  median time:      1.087 s (2.57% GC)
  mean time:        1.091 s (2.57% GC)
  maximum time:     1.110 s (2.62% GC)
  --------------
  samples:          5
  evals/sample:     1
```

Note: CPU: Intel(R) Xeon(R) Gold 6230 CPU @ 2.10GHz.

It corresponds to the second column of ADBench paper
![ADBench](benchmarks/adbench.png)
