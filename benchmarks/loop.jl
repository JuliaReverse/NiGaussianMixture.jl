using NiGaussianMixture
using StaticArrays
using BenchmarkTools

main_term = zeros(1, 5)
Qs = [MMatrix{10,10}(randn(10,10)) for i=1:5]
xs = [MVector{10}(randn(10)) for i=1:10000]
ms = [MVector{10}(randn(10)) for i=1:5]
alphas = randn(1, 5)
sum_qs = randn(1, 5)
@benchmark loop!(0.0, $main_term, $Qs, $xs, $ms, $alphas, $sum_qs, 10000)
