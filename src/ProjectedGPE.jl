__precompile__()

module ProjectedGPE

using Reexport, Combinatorics
#import

@reexport using DifferentialEquations
@reexport using FastGaussQuadrature
#@reexport using ApproxFun
@reexport using Parameters

@with_kw type CInfo
  Γ̄::Float64=1e-4
  M̄::Float64=0.0
  g::Float64=0.1
  μ::Float64=12.0
  basis::String="Hermite"
  Ω::Vector=[1]
  ecut::Number=10
  e0::Number=0.5
  Mmax::Number=10
  Mult::Number=10
  M::Vector=[10]
  en::Array=collect(1:10)-0.5
  P::BitArray=[true]
  Tx::Array{Float64,2}=eye(2,2)
  Ty::Array{Float64,2}=eye(2,2)
  Tz::Array{Float64,2}=eye(2,2)
  W::Array=randn(2,2,2)
  ψ::Array=randn(2,2,2)
end

@with_kw type Params @deftype Float64
#== paste this template into timeevolution and modify ==#
#fundamental constants/units
  ħ = 1.0545718e-34
  kB = 1.38064852e-23
  amu = 1.660339040e-27
  Bohr = 5.29e-11
#Rb87 mass and scattering length
  m = 86.909180527*amu
  as = 100*Bohr
#trap frequencies
  ωx = 2π
  ωy = 0.
  ωz = 0.
#choice of time and length units
  t0 = 1.0/ωx
  x0 = sqrt(ħ*t0/m)
  E0 = ħ/t0
#interactions
  g  = (4*pi*ħ^2*as/m)*x0^3/E0 #dimensionless 3D
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
include("makespecops.jl")
include("nenergy.jl")
include("ladderops.jl")
include("nfieldtrans.jl")
include("anisotrans.jl")
include("growthrate.jl")
include("scatteringkernel.jl")
include("makescatteringtrans.jl")
include("makescatteringnoisetrans.jl")
include("makecinfo.jl")
include("makealltrans.jl")
include("makespecops.jl")
include("pgpe.jl")
include("timeevolution.jl")
include("timeevolution2.jl")
#include("evalues.jl")
#putting these here for now (=> package VortexDistributions.jl)
include("unwrap.jl")
include("findvortices.jl")

export eigmat, nfieldtrans, anisotrans!, growthrate,
scatteringkernel, makescatteringtrans, makescatteringnoisetrans,
makecinfo, makealltrans, makespecops, ladderops, nenergy,
gausshermite, timeevolution, timeevolution2, CInfo, Params,
nlin!, Lgp!,findvortices, unwrap

end # module
