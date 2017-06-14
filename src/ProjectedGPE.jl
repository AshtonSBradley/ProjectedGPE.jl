__precompile__()

module ProjectedGPE

using Reexport
#import
@reexport using DifferentialEquations
@reexport using FastGaussQuadrature
@reexport using ApproxFun
@reexport using Parameters

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

@with_kw type Params @deftype Float64
#== paste this template into timeevolution and modify ==#
#Rb87 mass and scattering length
  ħ = 1.05457e-34
  m = 1.419e-25
  a = 5e-9
#trap frequencies
  ωx = 2π
  ωy = 0.
  ωz = 0.
#choice of time and length units
  t0 = 1.0/ωx
  x0 = sqrt(ħ*t0/m)
  E0 = ħ/t0
#interactions
  g  = (4*pi*ħ^2*a/m)*x0^3/E0 #dimensionless 3D
#damping parameters (dimensionless)
  Γ̄  = 1e-4; @assert Γ̄<=1.0
  M̄  = 0.0;  @assert M̄<=1.0
#chemical potential (dimensionless)
  μ  = 12.0
#time evolution parameters
  ti = 0.0
  tf = 100t0
  Nt::Int64 = 50
  t::Vector = collect(linspace(ti,tf,Nt))
  dt = 0.01π/μ #integrate step size [ - should have dt ≪ 2π/μ]
#== end template ==#
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
timeevolution2, findvortices, unwrap, CInfo, Params

end # module
