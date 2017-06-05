__precompile__()

module ProjectedGPE

#import
using FastGaussQuadrature, DifferentialEquations, Parameters

#abstract CField

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

export eigmat, nfieldtrans, anisotrans, growthrate, scatteringkernel, makescatteringtrans, makescatteringnoisetrans, makecfieldinfo, findvortices, unwrap


end # module
