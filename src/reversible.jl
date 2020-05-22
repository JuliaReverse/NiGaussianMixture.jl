NiLang.AD.GVar(wa::Wishart) = Wishart(GVar(wa.gamma), wa.m)

@nograd log_gamma_distrib(out!::GVar, a, p::Int)
@i function log_gamma_distrib(out!, a::T, p::Int) where T
	@routine @invcheckoff begin
		anc ← zero(T)
  		anc += (0.25 * log(π)) * p
	end
	p -= identity(1)
	out! += anc * p
	p += identity(1)
	@invcheckoff for j in 0:p-1
		@routine begin
			s ← 0.0
			j_half ← zero(T)
			j_half += (-0.5) * j
			j_half += identity(a)
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
    		res[r,i] += identity(LT[row_start + i])
		end
		@inbounds res[r, r] += identity(D[r])
  	end
end

@i function ascending!(xs!::AbstractVector{T}, inds!, arr::AbstractArray{T}) where T
	@invcheckoff if (length(arr) > 0, ~)
		y ← zero(T)
		y += identity(arr[1])
		ipush!(xs!, y)
		anc ← 1
		ipush!(inds!, anc)
		anc → 0
		@inbounds for i = 2:length(arr)
			if (arr[i] > xs![end], i==inds![end])
				ind ← i
				x ← zero(T)
				x += identity(arr[i])
				ipush!(xs!, x)
				ipush!(inds!, ind)
				ind → 0
			end
		end
	end
end

@i function isum(out!, x::AbstractArray)
	@invcheckoff for i=1:length(x)
		@inbounds out! += identity(x[i])
	end
end

@i function isum(out!, f, x::AbstractArray)
	@invcheckoff for i=1:length(x)
		@inbounds out! += f(x[i])
	end
end

@i function igemm!(out!::AbstractMatrix{T}, x::AbstractMatrix{T}, y::AbstractMatrix{T}) where T
	@safe size(x, 2) == size(y, 1) || throw(DimensionMismatch())
	@invcheckoff @inbounds for k=1:size(y,2)
		for i=1:size(x,1)
			for j=1:size(x,2)
				out![i,k] += x[i,j] * y[j,k]
			end
		end
	end
end

@i function igemv!(out!::AbstractVector{T}, x::AbstractMatrix, y::AbstractVector) where T
	@safe size(x, 2) == size(y, 1) || throw(DimensionMismatch())
	@invcheckoff @inbounds for j=1:size(x,2)
		yj ← zero(T)
		yj += identity(y[j])
		for i=1:size(x,1)
			out![i] += x[i,j] * yj
		end
		yj -= identity(y[j])
	end
end

# input should be 1 dimensional
@i function logsumexp(logout!, out!, xs!, inds!, x::AbstractArray{T}) where T
	mx ← zero(T)
  	ascending!(xs!, inds!, x)
	mx += identity(xs![end])
	@invcheckoff @inbounds for i=1:length(x)
		x[i] -= identity(mx)
		out! += exp(x[i])
		x[i] += identity(mx)
	end
  	logout! += log(out!)
	logout! += identity(mx)
	mx -= identity(xs![end])
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
	  	frobenius ← zero(T)
		logg ← zero(T)
		sq_gamma ← zero(T)
		half_sq_gamma ← zero(T)
		C ← zero(T)
		sum_sum_qs ← zero(T)

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
  		isum(frobenius, (@skip! abs2),icf[p+1:end,:])
		isum(sum_sum_qs, sum_qs)
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
			sum_qs[ic] += identity(icf[j,ic])
		end
	end
end

@i function gmm_objective(loss, alphas, means::AbstractMatrix{T}, icf::AbstractMatrix, xs::AbstractMatrix{XT}, wishart::Wishart) where {T,XT}
	@routine @invcheckoff begin
	  	d ← size(xs,1)
	  	n ← size(xs,2)
	  	m ← size(means,1)
	  	k ← size(means,2)
		sum_qs ← zeros(T,1,size(icf, 2))
		Qs ← zeros(T,d,d,k)
		loss1 ← zero(loss)
		loss2 ← zero(loss)
		salpha ← zero(loss)
		out_anc ← zero(loss)
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
	ipush!(salpha)
	ipush!(loss1)
	ipush!(loss2)
	ipush!(inds_anc)
	ipush!(out_anc)
	ipush!(xs_anc)
	ipush!(main_term)
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
						x[i,ix] -= identity(means[i,ik])
					end
					#igemv!(view(Outs, :, ik), view(Qs,:,:,ik), view(x,:,ix))

					# gemv
					for j=1:d
						xj ← zero(T)
						xj += identity(x[j,ix])
						for i=1:d
							Outs[i, ik] += Qs[i,j,ik] * xj
						end
						xj -= identity(x[j,ix])
					end
					#isum(local_out[ik], (@skip! abs2), view(Outs,:,ik))
					oi ← zero(T)
					for l=1:d
						 oi += abs2(Outs[l,ik])
					end
					SWAP(local_out[ik], oi)
	      			main_term[ik] -= 0.5 * local_out[ik]
	      			main_term[ik] += sum_qs[ik] + alphas[ik]
					for i=1:size(x, 1)
						x[i,ix] += identity(means[i,ik])
					end
				end
				logsumexp(logout, out, xs, inds, main_term)
    		end
			slse += identity(logout)
			~@routine
  		end
	end
end
