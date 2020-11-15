function log_gamma_distrib(a, p)
    out = 0.25 * p * (p - 1) * 1.1447298858494002 #convert(Float64, log(pi))
	for j in 1:p
        out += (logabsgamma(a + 0.5*(1 - j)))[1]
    end
	out
end

function ltri_unpack(D, LT::AbstractArray{T}) where T
    d=length(D)
	res = zeros(T, d, d)
	ltri_unpack!(res, D, LT::AbstractArray{T}) where T
end

function ltri_unpack!(res, D, LT::AbstractArray{T}) where T
    d=length(D)
    for r=1:d
	    row_start = div((r-1)*(r-2),2)
	    for i=1:r-1
    	    @inbounds res[r,i] = LT[row_start + i]
	    end
	    @inbounds res[r, r] = D[r]
    end
    res
end

# input should be 1 dimensional
function logsumexp(x)
    mx = maximum(x)
    log(sum(xi->exp(xi - mx), x)) + mx
end

function log_wishart_prior(wishart::Wishart, sum_qs, Qs, icf, d, k)
    p = size(Qs,1)
    n = p + wishart.m + 1
    C = n*p*(log(wishart.gamma) - 0.5*log(2)) - log_gamma_distrib(0.5*n, p)

    frobenius = 0.
    @inbounds for iq=1:size(Qs, 3)
        frobenius += sum(abs2,diag(view(Qs,:,:,iq)))
    end
    frobenius += sum(abs2,view(icf,d+1:size(icf, 1),:))
	0.5*wishart.gamma^2 * frobenius - wishart.m*sum(sum_qs) - k*C
end

function get_Q!(Q,d,icf)
  	ltri_unpack!(Q,exp.(view(icf,1:d)),view(icf,d+1:length(icf)))
end

function gmm_objective(alphas,means::AbstractArray{T},icf,x,wishart::Wishart) where T
  	d = size(x,1)
  	n = size(x,2)
  	m = size(means,1)
  	k = size(means,2)
  	CONSTANT = -n*d*0.5*log(2 * pi)

  	sum_qs = sum(view(icf,1:d,:),dims=1)
    Qs = zeros(T,d,d,k)
    @inbounds for ik = 1:k
        get_Q!(view(Qs,:,:,ik),d,view(icf,:,ik))
    end
	slse = loop!(Qs, x, means, alphas, sum_qs, n)
  	CONSTANT + slse - n*logsumexp(alphas) + log_wishart_prior(wishart, sum_qs, Qs, icf, d, k)
end

function loop!(Qs, x, means, alphas::AbstractArray{T}, sum_qs, n) where T
  	k = size(means, 2)
    d = size(means, 1)
  	main_term = zeros(T,1,k)
  	slse = 0.0
	work = zeros(T, d)
	workB = zeros(T, d)
  	@inbounds for ix=1:n
    	for ik=1:k
            work .= view(x,:,ix) .- view(means,:,ik)
            main_term[ik] = -0.5*sum(abs2, mul!(workB, view(Qs,:,:,ik), work))
    	end
    	slse += logsumexp(alphas .+ sum_qs .+ main_term)
  	end
  	slse
end

function pack(alphas,means,icf)
  	[alphas[:];means[:];icf[:]]
end

function unpack(d,k,packed)
  	alphas = reshape(packed[1:k],1,k)
  	off = k
  	means = reshape(packed[(1:d*k) .+ off],d,k)
  	icf_sz = div(d*(d + 1),2)
  	off += d*k
  	icf = reshape(packed[off+1:end],icf_sz,k)
  	(alphas,means,icf)
end

function jac_forwarddiff(alphas, means, icf, x, wishart)
  	d = size(x, 1)
  	k = size(alphas, 2)
  	function wrapper_gmm_objective(packed)
    	alphas,means,icf = unpack(d,k,packed)
    	gmm_objective(alphas,means,icf,x,wishart)
  	end
  	params = pack(alphas,means,icf)
  	ForwardDiff.gradient(wrapper_gmm_objective, params)
end
