# create ground state
# imprint vortex
# evolve in real time
using Revise, ProjectedGPE, VortexDistributions

using PyPlot

function groundstate(γ,tf,Ncut)
#== start parameters ==#
sinfo = Params()
#fundamental constants/units
  ħ = 1.0545718e-34; kB = 1.38064852e-23; amu = 1.660339040e-27; Bohr = 5.29e-11
#Rb87 mass and scattering length
  m = 86.909180527*amu; as = 100*Bohr
#trap frequencies
  ωx = 2π; ωy = ωx; ωz = 0.
#choice of time, length, energy units
  t0 = 1.0/ωy; x0 = sqrt(ħ*t0/m); E0 = ħ/t0
#interactions and chemical potential, e.g., g  = (4*pi*ħ^2*as/m)*x0^3/E0 (dimensionless 3D)
  g = 0.1; μ  = 10.0
#damping parameters (dimensionless)
  #γ = 0.5
  ℳ  = 0.0
#time evolution parameters
  ti = 0.0; #tf = 10.0/γ  #evolve for 2 damping times
  Nt = 50; t = collect(linspace(ti,tf,Nt)); dt = 0.01π/μ #integrate step size [ - should have dt ≪ 2π/μ]
#
@pack sinfo = ωx,ωy,ωz,γ,ℳ,g,t0,x0,E0,μ,ti,tf,Nt,t,dt
#== end parameters ==#

#= Initialize cfield (dimensionless) =#
  basis = "Hermite"
  ecut = Ncut*ħ*ωy/E0
  Ω = [ωx; ωy]*t0
  cinfo = makecinfo(ecut,Ω,basis)
  tinfo = maketinfo(cinfo,4)
  @unpack en,P,M = cinfo

#= create initial state and extra fields =#
  c0  = randn(M...)+im*randn(M...); c0=P.*c0
  ψ0  = c2x(c0,tinfo) #initial condition
  ψ = similar(ψ0)  #a field to write to in place

#= PGPE (in place)       =#
#add spectral terms into Lgp!, or x-space terms into nlin!
  function Lgp!(dc,c,p,t)
      nlin!(dc,c,tinfo)
      dc .= P.*(-im*(1-im*γ)*((en - μ).*c .+ g*dc))
  end
  function nlin!(dc,c,tinfo)
      c2x!(ψ,c,tinfo)
      x2c!(dc,abs2.(ψ).*ψ,tinfo)
  end
#========================#

tspan = (t[1],t[end])
prob = ODEProblem(Lgp!,c0,tspan)
alg =    DP5() #abstol = 1e-3, reltol = 1e-3

println("Started evolution ...")
#@time sol = solve(prob,alg,dt=dt,saveat=t,alg_hints=[:stiff],save_everystep=false,dense=false);
@time sol = solve(prob,alg,dt=dt,saveat=t);
println("... Finished.")
return sinfo,cinfo,sol
end

#EVOLVE
sinfo,cinfo,sol = groundstate(0.3,10/γ,30.5)

function addvortex()
@unpack ħ,m,ωx,ωy,ωz,γ,ℳ,g,x0,t0,E0,μ,ti,tf,Nt,t,dt = sinfo
#@unpack M,Ω,ecut,P,en = cinfo

#= Plot solution =#
#lets be explicit about units:
Rx = sqrt(2μ*E0/m/ωx^2)/x0
Ry = sqrt(2μ*E0/m/ωy^2)/x0
yMax=1.5Ry
xMax=1.5Rx

Nx = 600
Ny = Nx
x = collect(linspace(-xMax,xMax,Nx))
y = collect(linspace(-yMax,yMax,Ny))
tinfop = maketinfoplot(cinfo,x,y)
ψ = c2x(sol[end],tinfop)

#= imprint vortex =#
xv = 0.2*Rx; yv = 0.0; σv = 1
testvort = [xv yv σv]
dx = x[2]-x[1];dy = y[2]-y[1]

ix = find(abs.(x-xv).<=dx);ix = ix[1]
iy = find(abs.(y-yv).<=dy);iy = iy[1]
#ξ = x0/(sqrt(μ))
ξ = 1/(sqrt(g*abs2.(ψ)[ix,iy]))

makevortex!(ψ,testvort,x,y,ξ)
#project back to cfield
return xgrid2c(ψ,x,y,cinfo)
end

c0 = addvortex()

function realtime(tf,Ncut)

#== start parameters ==#
sinfo = Params()
#fundamental constants/units
  ħ = 1.0545718e-34; kB = 1.38064852e-23; amu = 1.660339040e-27; Bohr = 5.29e-11
#Rb87 mass and scattering length
  m = 86.909180527*amu; as = 100*Bohr
#trap frequencies
  ωx = 2π; ωy = ωx; ωz = 0.
#choice of time, length, energy units
  t0 = 1.0/ωy; x0 = sqrt(ħ*t0/m); E0 = ħ/t0
#interactions and chemical potential, e.g., g  = (4*pi*ħ^2*as/m)*x0^3/E0 (dimensionless 3D)
  g = 0.1; μ  = 10.0
#damping parameters (dimensionless)
  γ = 0.0
  ℳ  = 0.0
#time evolution parameters
  ti = 0.0; #tf = 10.0/γ  #evolve for 2 damping times
  Nt = 50; t = collect(linspace(ti,tf,Nt)); dt = 0.01π/μ #integrate step size [ - should have dt ≪ 2π/μ]
#
@pack sinfo = ωx,ωy,ωz,γ,ℳ,g,t0,x0,E0,μ,ti,tf,Nt,t,dt
#== end parameters ==#

#= Initialize cfield (dimensionless) =#
  basis = "Hermite"
  ecut = Ncut*ħ*ωy/E0
  Ω = [ωx; ωy]*t0
  cinfo = makecinfo(ecut,Ω,basis)
  tinfo = maketinfo(cinfo,4)
  @unpack en,P,M = cinfo

#= create initial state and extra fields =#
  ψ0  = c2x(c0,tinfo) #initial condition
  ψ = similar(ψ0)  #a field to write to in place

#= PGPE (in place)       =#
#add spectral terms into Lgp!, or x-space terms into nlin!
  function Lgp!(dc,c,p,t)
      nlin!(dc,c,tinfo)
      dc .= P.*(-im*(1-im*γ)*((en - μ).*c .+ g*dc))
  end
  function nlin!(dc,c,tinfo)
      c2x!(ψ,c,tinfo)
      x2c!(dc,abs2.(ψ).*ψ,tinfo)
  end
#========================#

tspan = (t[1],t[end])
prob = ODEProblem(Lgp!,c0,tspan)
alg =    DP5() #abstol = 1e-3, reltol = 1e-3

println("Started evolution ...")
#@time sol = solve(prob,alg,dt=dt,saveat=t,alg_hints=[:stiff],save_everystep=false,dense=false);
@time sol = solve(prob,alg,dt=dt,saveat=t);
println("... Finished.")
return sinfo,cinfo,sol
end

#should probably be careful not to write over initial state here
sinfo2,cinfo2,sol2 = realtime(60.,30.5)

function findvortex()
## Transform to cartesian grid
@unpack ħ,m,ωx,ωy,ωz,γ,ℳ,g,x0,t0,E0,μ,ti,tf,Nt,t,dt = sinfo2
#@unpack M,Ω,ecut,P,en = cinfo

#lets be explicit about units:
Rx = sqrt(2μ*E0/m/ωx^2)/x0
Ry = sqrt(2μ*E0/m/ωy^2)/x0
yMax=1.5Ry
xMax=1.5Rx

Nx = 600
Ny = Nx
x = collect(linspace(-xMax,xMax,Nx))
y = collect(linspace(-yMax,yMax,Ny))
tinfop = maketinfoplot(cinfo,x,y)
ψ = c2x(sol2[1],tinfop)
mask = complex(Float64.(sqrt.(x.^2.+y'.^2).<0.7*Rx*ones(real(ψ))))

xt = zeros(t)
yt = zeros(t)
for j in eachindex(t)
    c2x!(ψ,sol2[j],tinfop)
    xt[j], yt[j], _ = findvortices(x,y,ψ.*mask)
end
return xt, yt
end

xt,yt = findvortex()

@unpack ħ,m,ωx,ωy,ωz,γ,ℳ,g,x0,t0,E0,μ,ti,tf,Nt,t,dt = sinfo2
plot(t,xt)
plot(t,yt)
#=
figure(figsize=(15,5))
subplot(1,2,1)
pcolormesh(x,y,g*abs2.(ψ'))
xlabel("x/x0")
ylabel("y/x0")
colorbar()
subplot(1,2,2)
pcolormesh(x,y,angle.(ψ'))
xlabel("x/x0")
ylabel("y/x0")
colorbar()

figure(figsize=(15,5))
subplot(1,2,1)
pcolormesh(x,y,g*abs2.(ψ'))
xlabel("x/x0")
ylabel("y/x0")
colorbar()
subplot(1,2,2)
pcolormesh(x,y,angle.(ψ'))
xlabel("x/x0")
ylabel("y/x0")
colorbar()

figure()
plot(x,g*abs2.(ψ[:,300]))
plot(x,μ*(1-x.^2/Rx^2).*vortexcore(x-xv,ξ).^2)
=#
