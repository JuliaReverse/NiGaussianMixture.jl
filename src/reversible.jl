NiLang.AD.GVar(wa::Wishart) = Wishart(GVar(wa.gamma), wa.m)

@nograd log_gamma_distrib(out!::GVar, a, p::Int)
@i function log_gamma_distrib(out!, a::T, p::Int) where T
	@routine @invcheckoff begin
		anc ← zero(T)
  		anc += (0.25 * log(π)) * p
	end
	p -= 1
	out! += anc * p
	p += 1
	@invcheckoff for j in 0:p-1
		@routine begin
			s ← 0.0
			j_half ← zero(T)
			j_half += (-0.5) * j
			j_half += a
		end
    	(out!, s) += logabsgamma(j_half)
		s -= sign(j_half)
		~@routine
  	end
	~@routine
end

@i function ltri_unpack(res, D, LT::AbstractArray{T}) where T
  	d ← length(D)
  	@invcheckoff for r = 1:d
		row_start ← (r-1) * (r-2) ÷ 2
		@inbounds for i=1:r-1
    		res[r,i] += LT[row_start + i]
		end
		@inbounds res[r, r] += D[r]
  	end
end

# input should be 1 dimensional
@i function logsumexp(logout!, out!, xs!, inds!, x::AbstractArray{T}) where T
	mx ← zero(T)
  	i_ascending!(xs!, inds!, x)
	mx += xs![end]
	@invcheckoff @inbounds for i=1:length(x)
		x[i] -= mx
		out! += exp(x[i])
		x[i] += mx
	end
  	logout! += log(out!)
	logout! += mx
	mx -= xs![end]
end

function NiLang.chfield(a::AbstractArray{T}, ::typeof(diag), val::AbstractVector{T}) where T
	for i=1:min(length(val), size(a)...)
		@inbounds a[i,i] = val[i]
	end
	return a
end

@i function log_wishart_prior(out!, wishart::Wishart{T}, sum_qs, Qs, icf) where T
	@routine @invcheckoff begin
	  	p ← size(Qs, 1)
	  	k ← size(Qs, 3)
	  	n ← 1
		np ← 0
		halfn ← 0.0
        @zeros T frobenius logg sq_gamma half_sq_gamma C sum_sum_qs

		sq_gamma += wishart.gamma^2
		half_sq_gamma += 0.5 * sq_gamma
		logg += log(wishart.gamma)
		logg -= 0.5 * log(2)
		n += p + wishart.m
		halfn += 0.5 * n
		np += n * p
  		C += np * logg
		log_gamma_distrib(-C, halfn, p)

  		for iq = 1:size(Qs,3)
			for l=1:p
    			@inbounds frobenius += abs2(Qs[l,l,iq])
			end
  		end
  		i_sum(frobenius, (@skip! abs2),icf[p+1:end,:])
		i_sum(sum_sum_qs, sum_qs)
	end

	out! += half_sq_gamma * frobenius
	out! -= wishart.m * sum_sum_qs
	out! -= k * C
	~@routine
end

@i function get_Q(res::AbstractMatrix, icf::AbstractArray{T}) where T
	@routine @invcheckoff begin
		d ← size(res,1)
	  	expd ← zeros(T, d)
	  	expd .+= exp.(icf[1:d])
	end
  	ltri_unpack(res, expd, icf[d+1:end])
	~@routine
end

@i function sum_row(sum_qs,icf,d)
	@invcheckoff @inbounds for ic = 1:size(icf, 2)
		for j=1:d
			sum_qs[ic] += icf[j,ic]
		end
	end
end

@i function gmm_objective(loss::T, alphas, means::AbstractMatrix{T}, icf::AbstractMatrix, xs::AbstractMatrix{XT}, wishart::Wishart) where {T,XT}
	@routine @invcheckoff begin
	  	d ← size(xs,1)
	  	n ← size(xs,2)
	  	m ← size(means,1)
	  	k ← size(means,2)
		sum_qs ← zeros(T,1,size(icf, 2))
		Qs ← zeros(T,d,d,k)
        @zeros T loss1 loss2 salpha out_anc
		xs_anc ← T[]
		inds_anc ← Int[]
	end
	@invcheckoff main_term ← zeros(T,1,k)
	@invcheckoff @routine begin
		sum_row(sum_qs, icf, @keep d)
		for ik = 1:k
	  		get_Q(view(Qs,:,:,ik), view(icf,:,ik))
		end
	end
	loop!(loss1, main_term, Qs, xs, means, alphas, sum_qs, n)
	log_wishart_prior(loss2, wishart, sum_qs, Qs, icf)
	(~log_wishart_prior)(loss2, wishart, sum_qs, Qs, icf)
	logsumexp(salpha, out_anc, xs_anc, inds_anc, alphas)
	~@routine

    loss -= identity(n*d*0.5*log(2 * pi))
	loss += loss1 + loss2
  	loss -= n * salpha
	PUSH!(salpha)
	PUSH!(loss1)
	PUSH!(loss2)
	PUSH!(inds_anc)
	PUSH!(out_anc)
	PUSH!(xs_anc)
	PUSH!(main_term)
	@invcheckoff main_term → zeros(T,0,0)
end

@i function loop!(slse::T, main_term, Qs, x::AbstractArray, means, alphas::AbstractArray, sum_qs, n) where T
  	@invcheckoff begin
	  	k ← size(means, 2)
		d ← size(Qs, 1)
		Outs ← zeros(T, d, k)
		logout ← zero(T)
		local_out ← zeros(T, k)
		out ← zero(T)
		xs ← T[]
		inds ← Int[]
		for ix=1:n
			@routine begin
    			@inbounds for ik=1:k
					for i=1:d
						x[i,ix] -= means[i,ik]
					end

					# gemv
					for j=1:d
						xj ← zero(T)
						xj += x[j,ix]
						for i=1:d
							Outs[i, ik] += Qs[i,j,ik] * xj
						end
						xj -= x[j,ix]
					end
					oi ← zero(T)
					for l=1:d
						 oi += abs2(Outs[l,ik])
					end
					SWAP(local_out[ik], oi)
	      			main_term[ik] -= 0.5 * local_out[ik]
	      			main_term[ik] += sum_qs[ik] + alphas[ik]
					for i=1:size(x, 1)
						x[i,ix] += means[i,ik]
					end
				end
				logsumexp(logout, out, xs, inds, main_term)
    		end
			slse += logout
			~@routine
  		end
	end
end
