using BenchmarkTools
using NiGaussianMixture

alphas, means, icf, x, wishart = load(10, 5)
@benchmark gmm_objective(0.0, $alphas,$means,$icf,$x,$wishart)
@benchmark gmm_objective($alphas,$means,$icf,$x,$wishart)
@benchmark gmm_objective'(Val(1), 0.0, $alphas,$means,$icf,$x,$wishart)
