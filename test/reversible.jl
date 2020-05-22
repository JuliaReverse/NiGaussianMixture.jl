using Test
using NiGaussianMixture: logsumexp, ltri_unpack, ascending!, log_gamma_distrib
using NiLang
using SpecialFunctions
using ForwardDiff, NiLang.AD

function load(d, k)
    dir_in = joinpath(dirname(dirname(@__FILE__)), "data", "gmm")
    dir_out = dirname(@__FILE__)
    fn = joinpath("10k", "gmm_d$(d)_K$k")

    fn_in = joinpath(dir_in, fn)
    fn_out = joinpath(dir_out, fn)
    NiGaussianMixture.read_gmm_instance(string(fn_in,".txt"), false)
end

@testset "logsumexp" begin
	function logsumexp2(x)
	  	mx = maximum(x)
	  	log.(sum(exp.(x .- mx))) .+ mx
	end

	x = randn(100)
	@test ascending!(Float64[], Int[], x)[1][end] == maximum(x)
	@test logsumexp(0.0, 0.0, Float64[], Int[], x)[1] ≈ logsumexp2(x)
	@test check_inv(logsumexp, (0.0, 0.0, Float64[], Int[], x))
end

@testset "get_Q" begin
	function ltri_unpack2(D, LT::AbstractArray{T}) where T
	  	d=length(D)
	  	res = zeros(T, d, d)
	  	for r=1:d
			row_start = div((r-1)*(r-2),2)
			for i=1:r-1
	    		@inbounds res[r,i] = LT[row_start + i]
			end
			@inbounds res[r, r] = D[r]
	  	end
	  	res
	end
	function get_Q2(d,icf)
	  	ltri_unpack2(exp.(icf[1:d]),icf[d+1:end])
	end

	icf = randn(100)
	@test NiGaussianMixture.get_Q(zeros(5, 5), icf)[1] ≈ get_Q2(5, icf)
	@test check_inv(NiGaussianMixture.get_Q, (zeros(5, 5), icf))
end

@testset "log gamma distri" begin
	function log_gamma_distrib2(a, p)
  		out = 0.25 * p * (p - 1) * 1.1447298858494002 #convert(Float64, log(pi))
		for j in 1:p
    		out += (logabsgamma(a + 0.5*(1 - j)))[1]
  		end
		out
	end
	@test log_gamma_distrib(0.0, 0.7, 2)[1] ≈ log_gamma_distrib2(0.7, 2)
end

@testset "gmm objective" begin
	empty!(NiLang.GLOBAL_STACK)
	alphas, means, icf, x, wishart = load(10, 5)
	@test isapprox(gmm_objective(0.0,alphas,means,icf,x,wishart)[1], -317785.8722759027; atol=0.1)
	@test check_inv(gmm_objective, (0.0,alphas,means,icf,x,wishart); verbose=true)
	empty!(NiLang.GLOBAL_STACK)
end

@testset "forward diff" begin
	alphas, means, icf, x, wishart = load(10, 5)
	jf = NiGaussianMixture.jac_forwarddiff(alphas, means, icf, x, wishart)
	_, _, galphas,gmeans,gicf,_,_ = Grad(gmm_objective)(Val(1), 0.0, alphas,means,icf,x,wishart)
  	gparams = NiGaussianMixture.pack(grad(galphas), grad(gmeans), grad(gicf))
	@test isapprox(gparams, jf, rtol=1e-2)
end

@testset "gmm objective, irreversible" begin
    alphas, means, icf, x, wishart = load(10, 5)
    @test isapprox(gmm_objective(alphas,means,icf,x,wishart), -317178.3340912648; atol=0.1)
end
