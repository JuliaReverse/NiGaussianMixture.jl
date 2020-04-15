using BenchmarkTools
using NiGaussianMixture

function do_benchmark(alphas, means, icf, x, wishart)
  d = size(x, 1)
  k = size(alphas, 2)
  function wrapper_gmm_objective(packed)
    alphas,means,icf = unpack(d,k,packed)
    gmm_objective(alphas,means,icf,x,wishart)
  end
  params = pack(alphas,means,icf)
  println("nparams = $(length(params))")
  display(@benchmark ForwardDiff.gradient($wrapper_gmm_objective, $(params)))
end

alphas, means, icf, x, wishart = load(10, 5)
do_benchmark(alphas, means, icf, x, wishart)

alphas, means, icf, x, wishart = load(10, 5)
do_benchmark(alphas, means, icf, x, wishart)

# Objective
# Call once in case of precompilation etc
@benchmark gmm_objective($alphas,$means,$icf,$x,$wishart)
