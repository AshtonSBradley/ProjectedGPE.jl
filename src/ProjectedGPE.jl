__precompile__()

module ProjectedGPE

#import
using FastGaussQuadrature, DifferentialEquations, ApproxFun, Parameters

@with_kw type CInfo
  basis::String="Hermite"
  Ω::Vector=[1]
  ecut::Number=10
  e0::Number=0.5
  Mmax::Number=10
  M::Vector=[10]
  en::Array=collect(1:10)-0.5
  P::BitArray=[true]
end

include("eigmat.jl")
include("nfieldtrans.jl")
include("anisotrans.jl")
include("growthrate.jl")
include("scatteringkernel.jl")
include("makescatteringtrans.jl")
include("makescatteringnoisetrans.jl")
include("makecinfo.jl")
include("maketransinfo.jl")
include("timeevolution.jl")
#include("evalues.jl")
#putting these here for now (=> package VortexDistributions.jl)
include("unwrap.jl")
include("findvortices.jl")

export eigmat, nfieldtrans, anisotrans, growthrate,
scatteringkernel, makescatteringtrans, makescatteringnoisetrans,
makecinfo, maketransinfo, gausshermite, @pack, @unpack, timeevolution,
findvortices, unwrap

end # module
