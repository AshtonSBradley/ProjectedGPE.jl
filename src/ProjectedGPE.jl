__precompile__()

module ProjectedGPE

#import
using FastGaussQuadrature, DifferentialEquations

abstract cField

include("eigmat.jl")
include("nfieldtrans.jl")
include("unwrap.jl") #putting this here for now (=> package VortexDetect.jl)

export eigmat, nfieldtrans, unwrap


end # module
