#test time evolution (DifferentialEquations 4.0)
using Revise, ProjectedGPE

function timeevolution2()
siminfo = Params()

#== start template ==#
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
  tf = 1.0/Γ̄  #evolve for 2 damping times
  Nt = 50
  t = collect(linspace(ti,tf,Nt))
  dt = 0.01π/μ #integrate step size [ - should have dt ≪ 2π/μ]
#== end template ==#

@pack siminfo = ωx,ωy,ωz,γ,ℳ,g,t0,x0,E0,μ,ti,tf,Nt,t,dt

#Initialize CField (dimensionless)
  basis = "Hermite"
  ecut = 30*ħ*ωy/E0
  Ω = [ωx; ωy]*t0
  cinfo = makecinfo(ecut,Ω,basis)
  @unpack en,P,M = cinfo ;Mx,My = M
  x,wx,Tx,y,wy,Ty = makealltrans(M,4,Ω)
  Wxy = wx.*wy'
#test transform
  c0   = randn(Mx,My)+im*randn(Mx,My); c0=P.*c0
  ψ0   = Tx*c0*Ty' #initial condition
  ψ    = Tx*c0*Ty' #a field to write to in place

#PGPE time evolution
#out of place
function nlin(c)
    ψ = Tx*c*Ty'
    Tx'*(Wxy.*abs2.(ψ).*ψ)*Ty
end

#in place
function nlin!(dc,c)
    ψ = Tx*c*Ty'
    dc.= Tx'*(Wxy.*abs2.(ψ).*ψ)*Ty
end

#dPGPE in reservoir "frame"
#out of place
function Lgp(c,p,t)
     P.*(-im*(1-im*γ)*((en - μ).*c .+ g*nlin(c)))
end

#in place
function Lgp!(dc,c,p,t)
    nlin!(dc,c)
    dc .= P.*(-im*(1-im*γ)*((en - μ).*c .+ g*dc))
end

c0    = P.*(randn(Mx,My) + im*randn(Mx,My)) #create random initial state
tspan = (t[1],t[end])
prob = ODEProblem(Lgp!,c0,tspan)
#alg = Tsit5()
alg =    DP5()
#abstol = 1e-3
#reltol = 1e-3
println("Started evolution ...")
#@time sol = solve(prob,alg,dt=dt,saveat=t,alg_hints=[:stiff],save_everystep=false,dense=false);
@time sol = solve(prob,alg,dt=dt,saveat=t);
println("... Finished.")
return siminfo,cinfo,sol
end

siminfo,cinfo,sol = timeevolution2()

#plot solution for 2D
## Transform to cartesian grid
@unpack ħ,m,ωx,ωy,ωz,γ,ℳ,g,x0,t0,E0,μ,ti,tf,Nt,t,dt = siminfo
@unpack M,Ω,ecut,P,en = cinfo; Mx,My = M;

#examine solution
using Interact, PyPlot

#lets be explicit about units:
Rx = sqrt(2μ*E0/m/ωx^2)/x0
Ry = sqrt(2μ*E0/m/ωy^2)/x0
yMax=1.5Ry
xMax=1.5Rx

Nx = 400
Ny = Nx
x = collect(linspace(-xMax,xMax,Nx))
y = collect(linspace(-yMax,yMax,Ny))
Tx = eigmat("Hermite",Mx,x,ωx/ωy)
Ty = eigmat("Hermite",My,y,1); #units of ωy
#θ = unwrap(angle(ψ));

#Plot
f=figure(figsize=(12,3))
@manipulate for i=1:length(t) withfig(f,clear=true) do
    ψ = Tx*sol[i]*Ty';
    pcolormesh(x,y,g*abs2.(ψ'))
    xlabel("x/x0")
    ylabel("y/x0")
    colorbar()
    end
end
