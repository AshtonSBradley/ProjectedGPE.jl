__precompile__()

module ProjectedGPE

#using Reexport
#import
using DifferentialEquations
using FastGaussQuadrature
using ApproxFun
using Parameters

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

@with_kw type SimParams @deftype Float64
  ωx = 2π
  ωy = 0.
  ωz = 0.
  t0 = 1/ωx
  Γ̄  = 1e-4; @assert Γ̄<=1.0
  M̄  = 0.0;  @assert M̄<=1.0
  g  = 0.1
  μ  = 10.0
  ti = 0.0
  tf = 100t0
  Nt::Int64 = 50
  t::Vector = collect(linspace(ti,tf,Nt))
  dt = 0.1π/μ
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
include("timeevolution2.jl")
#include("evalues.jl")
#putting these here for now (=> package VortexDistributions.jl)
include("unwrap.jl")
include("findvortices.jl")

export eigmat, nfieldtrans, anisotrans, growthrate,
scatteringkernel, makescatteringtrans, makescatteringnoisetrans,
makecinfo, maketransinfo, gausshermite, @pack, @unpack, timeevolution,
timeevolution2, findvortices, unwrap

end # module
