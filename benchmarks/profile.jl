using BenchmarkTools
using NiGaussianMixture

function load(d, k)
    dir_in = joinpath(dirname(dirname(@__FILE__)), "data", "gmm")
    dir_out = dirname(@__FILE__)
    fn = joinpath("10k", "gmm_d$(d)_K$k")

    fn_in = joinpath(dir_in, fn)
    fn_out = joinpath(dir_out, fn)
    println("loading data: $(fn_in).txt")
    NiGaussianMixture.read_gmm_instance(string(fn_in,".txt"), false)
end

alphas, means, icf, x, wishart = load(64, 10)
using Profile
Profile.clear()
@profile gmm_objective(0.0, alphas, means, icf, x, wishart)
Profile.print(mincount=10, format=:flat)
