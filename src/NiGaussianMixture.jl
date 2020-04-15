module NiGaussianMixture

using SpecialFunctions
using ForwardDiff
using StaticArrays
using LinearAlgebra
using NiLang, NiLang.AD

struct Wishart{GT}
  	gamma::GT
  	m::Int
end

export gmm_objective

include("reversible.jl")
include("forwarddiff.jl")
include("load.jl")

end # module
