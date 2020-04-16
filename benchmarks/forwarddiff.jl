using BenchmarkTools
using ForwardDiff
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

println("Objective")
display(@benchmark gmm_objective($alphas,$means,$icf,$x,$wishart))

println("Gradients")
function benchmark_forwarddiff(alphas, means, icf, x, wishart)
    d = size(x, 1)
    k = size(alphas, 2)
    function wrapper_gmm_objective(packed)
        alphas,means,icf = NiGaussianMixture.unpack(d,k,packed)
        gmm_objective(alphas,means,icf,x,wishart)
    end
    params = NiGaussianMixture.pack(alphas,means,icf)
    println("nparams = $(length(params))")
    display(@benchmark ForwardDiff.gradient($wrapper_gmm_objective, $(params)))
end
benchmark_forwarddiff(alphas, means, icf, x, wishart)
println()
