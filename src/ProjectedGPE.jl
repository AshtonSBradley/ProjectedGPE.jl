__precompile__()

module ProjectedGPE

#import
using FastGaussQuadrature, DifferentialEquations

abstract cField

include("eigmat.jl")
include("nfieldtrans.jl")
include("unwrap.jl") #putting this here for now (=> package VortexDetect.jl)
include("scatteringkernel.jl")
include("findvortices.jl")

export eigmat, nfieldtrans, unwrap, scatteringkernel, findvortices


end # module
