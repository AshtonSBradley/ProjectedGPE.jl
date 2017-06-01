__precompile__()

module ProjectedGPE

#import
using FastGaussQuadrature, DifferentialEquations

abstract cField

include("eigmat.jl")
include("nfieldtrans.jl")
include("sgpegamma.jl")
include("scatteringkernel.jl")
include("makeVtrans.jl")

#putting this here for now (=> package VortexDetect.jl)
include("unwrap.jl")
include("findvortices.jl")

export eigmat, nfieldtrans, sgpegamma, scatteringkernel, makeVtrans, findvortices, unwrap


end # module
