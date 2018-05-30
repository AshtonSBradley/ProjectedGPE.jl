using Revise, ProjectedGPE

ħ = 1.0545718e-34
kB = 1.38064852e-23
amu = 1.660339040e-27
Bohr = 5.29e-11
#Rb87 mass and scattering length
m = 86.909180527*amu
as = 100*Bohr
#trap frequencies
ωx = 2π
ωy = 4ωx
ωz = 0.
#choice of time, length, energy units
t0 = 1.0/ωy
x0 = sqrt(ħ*t0/m)
E0 = ħ/t0
#interactions
#g  = (4*pi*ħ^2*as/m)*x0^3/E0 #dimensionless 3D
g = 0.1 #test 2D
#damping parameters (dimensionless)
γ = 0.05
ℳ  = 0.0
#chemical potential (dimensionless)
μ  = 12.0
#time evolution parameters
ti = 0.0
tf = 1.0/γ  #evolve for 2 damping times
Nt = 50
t = collect(linspace(ti,tf,Nt))
dt = 0.01π/μ #integrate step size [ - sh

basis = "Hermite"
ecut = 30*ħ*ωy/E0
Ω = [ωx; ωy]*t0
cinfo = makecinfo(ecut,Ω,basis)
tinfo = maketinfo(cinfo,4)
@unpack en,P,M = cinfo ;Mx,My = M
 c0   = randn(Mx,My)+im*randn(Mx,My); c0=P.*c0
ψ = tinfo.Tx*c0*tinfo.Ty'

function nlin!(dc,c,tinfo)
    #ψ = Tx*c*Ty'
    #dc.= Tx'*(W.*abs2.(ψ).*ψ)*Ty
    c2x!(ψ,c,tinfo)
    x2c!(dc,abs2.(ψ).*ψ,tinfo)
end

function Lgp!(dc,c,p,t)
    nlin!(dc,c,tinfo)
    dc .= P.*(-im*(1-im*γ)*((en - μ).*c .+ g*dc))
end
tinfo

dc=c0
nlin!(dc,c0,tinfo)
c2x!(ψ,c0,tinfo)
x2c!(dc,abs2.(ψ).*ψ,tinfo)
