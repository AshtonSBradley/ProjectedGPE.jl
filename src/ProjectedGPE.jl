__precompile__()

module ProjectedGPE

#import
using FastGaussQuadrature, DifferentialEquations, Parameters

@with_kw type cfieldinfo
  basis::String="hermite"
  ecut::Number
  e0::Number
  Espec::Vector
  N::Number
  N1::Number
  ω1::Number
  Px::BitArray
end

include("eigmat.jl")
include("nfieldtrans.jl")
include("anisotrans.jl")
include("growthrate.jl")
include("scatteringkernel.jl")
include("makescatteringtrans.jl")
include("makescatteringnoisetrans.jl")
include("makecfieldinfo.jl")
#include("maketransinfo.jl")
#include("evalues.jl")
#putting these here for now (=> package VortexDistributions.jl)
include("unwrap.jl")
include("findvortices.jl")

export eigmat, nfieldtrans, anisotrans, growthrate,
scatteringkernel, makescatteringtrans, makescatteringnoisetrans,
makecfieldinfo, gausshermite, findvortices, unwrap


end # module
