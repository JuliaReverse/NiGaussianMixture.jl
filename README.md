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
Normal Objective
BenchmarkTools.Trial: 
  memory estimate:  18.69 MiB
  allocs estimate:  159061
  --------------
  minimum time:     11.712 ms (0.00% GC)
  median time:      12.576 ms (0.00% GC)
  mean time:        13.073 ms (4.32% GC)
  maximum time:     19.791 ms (24.39% GC)
  --------------
  samples:          383
  evals/sample:     1
Reversible Objective
BenchmarkTools.Trial: 
  memory estimate:  967.94 KiB
  allocs estimate:  20119
  --------------
  minimum time:     14.487 ms (0.00% GC)
  median time:      14.834 ms (0.00% GC)
  mean time:        15.059 ms (0.37% GC)
  maximum time:     23.149 ms (0.00% GC)
  --------------
  samples:          332
  evals/sample:     1
NiLang Gradient
BenchmarkTools.Trial: 
  memory estimate:  3.74 MiB
  allocs estimate:  40271
  --------------
  minimum time:     37.693 ms (0.00% GC)
  median time:      38.862 ms (0.00% GC)
  mean time:        40.315 ms (0.71% GC)
  maximum time:     50.516 ms (18.68% GC)
  --------------
  samples:          124
  evals/sample:     1
ForwardDiff Gradient
nparams = 330
BenchmarkTools.Trial: 
  memory estimate:  3.40 GiB
  allocs estimate:  4454076
  --------------
  minimum time:     1.319 s (6.45% GC)
  median time:      1.435 s (11.68% GC)
  mean time:        1.459 s (12.34% GC)
  maximum time:     1.647 s (18.67% GC)
  --------------
  samples:          4
  evals/sample:     1
```

Note: CPU: Intel(R) Xeon(R) Gold 6230 CPU @ 2.10GHz.

It corresponds to the second column of ADBench paper
![ADBench](benchmarks/adbench.png)

## Make NiLang even faster!
We call for an efficient implementation of linear algebra!
The following routine cost `2/3` of NiLang's computing time, optimizing it will speed up the program by a factor of 2!

```julia
julia> A, b, v = randn(64, 64), randn(64), randn(64);

julia> @benchmark BLAS.gemv!('N', 1.0, $A, $b, 1.0, $v)
BenchmarkTools.Trial:
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     499.487 ns (0.00% GC)
  median time:      585.426 ns (0.00% GC)
  mean time:        589.599 ns (0.00% GC)
  maximum time:     4.237 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     195

julia> @benchmark NiGaussianMixture.igemv!($v, $A, $b)
BenchmarkTools.Trial:
  memory estimate:  32 bytes
  allocs estimate:  1
  --------------
  minimum time:     1.834 μs (0.00% GC)
  median time:      1.863 μs (0.00% GC)
  mean time:        1.955 μs (0.00% GC)
  maximum time:     7.121 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     10
```
