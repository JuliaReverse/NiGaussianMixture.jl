using BenchmarkTools
using NiGaussianMixture

function load(d, k)
    dir_in = joinpath(dirname(dirname(@__FILE__)), "data", "gmm")
    dir_out = dirname(@__FILE__)
    fn = joinpath("10k", "gmm_d$(d)_K$k")

    fn_in = joinpath(dir_in, fn)
    fn_out = joinpath(dir_out, fn)
    NiGaussianMixture.read_gmm_instance(string(fn_in,".txt"), false)
end

alphas, means, icf, x, wishart = load(10, 5)
do_benchmark(alphas, means, icf, x, wishart)

# Objective
# Call once in case of precompilation etc
@benchmark gmm_objective($alphas,$means,$icf,$x,$wishart)
