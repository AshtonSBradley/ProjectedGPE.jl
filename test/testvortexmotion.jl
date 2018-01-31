# create ground state
# imprint vortex
# evolve in real time
using Revise, ProjectedGPE, VortexDistributions

ENV["MPLBACKEND"]="tkagg"
using PyPlot

function timeevolution2()
sinfo = Params()

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
  ωy = ωx
  ωz = 0.
#choice of time, length, energy units
  t0 = 1.0/ωy
  x0 = sqrt(ħ*t0/m)
  E0 = ħ/t0
#interactions
  #g  = (4*pi*ħ^2*as/m)*x0^3/E0 #dimensionless 3D
  g = 0.1 #test 2D
#damping parameters (dimensionless)
  γ = 0.5
  ℳ  = 0.0
#chemical potential (dimensionless)
  μ  = 10.0
#time evolution parameters
  ti = 0.0
  tf = 10.0/γ  #evolve for 2 damping times
  Nt = 50
  t = collect(linspace(ti,tf,Nt))
  dt = 0.01π/μ #integrate step size [ - should have dt ≪ 2π/μ]
#== end template ==#

@pack sinfo = ωx,ωy,ωz,γ,ℳ,g,t0,x0,E0,μ,ti,tf,Nt,t,dt

#Initialize CField (dimensionless)
  basis = "Hermite"
  ecut = 30.5*ħ*ωy/E0
  Ω = [ωx; ωy]*t0
  cinfo = makecinfo(ecut,Ω,basis)
  tinfo = maketinfo(cinfo,4)
  @unpack en,P,M = cinfo
  #x,wx,Tx,y,wy,Ty = makealltrans(M,4,Ω)
  #W = wx.*wy'
#test transform
  c0  = randn(M...)+im*randn(M...); c0=P.*c0
  ψ0  = c2x(c0,tinfo) #initial condition
  ψ = similar(ψ0)  #a field to write to in place

#equation of motion (in place)
  function nlin!(dc,c,tinfo)
      c2x!(ψ,c,tinfo)
      x2c!(dc,abs2.(ψ).*ψ,tinfo)
  end

  function Lgp!(dc,c,p,t)
      nlin!(dc,c,tinfo)
      dc .= P.*(-im*(1-im*γ)*((en - μ).*c .+ g*dc))
  end

tspan = (t[1],t[end])
prob = ODEProblem(Lgp!,c0,tspan)
alg =    DP5() #abstol = 1e-3, reltol = 1e-3

println("Started evolution ...")
#@time sol = solve(prob,alg,dt=dt,saveat=t,alg_hints=[:stiff],save_everystep=false,dense=false);
@time sol = solve(prob,alg,dt=dt,saveat=t);
println("... Finished.")
return sinfo,cinfo,sol
end

sinfo,cinfo,sol = timeevolution2()

#plot solution for 2D
## Transform to cartesian grid
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
#= imprint vortex =#
xv = 0.2*Rx; yv = 0.0; σv = 1
testvort = [xv yv σv]
Δx = x[2]-x[1];Δy = y[2]-y[1]

ix = find(abs.(x-xv).<=Δx);ix = ix[1]
iy = find(abs.(y-yv).<=Δy);iy = iy[1]
#ξ = x0/(sqrt(μ))
ξ = 1/(sqrt(g*abs2.(ψ)[ix,iy]))
#g*abs2.(ψ)[ix,iy]

makevortex!(ψ,testvort,x,y',ξ)

#project back to cfield
c0 = xgrid2c(ψ,x,y,cinfo)

sum(abs2.(sol[end]))

sum(abs2.(c0))

ψ0 = c2x(c0,tinfop)
figure(figsize=(15,5))
subplot(1,2,1)
pcolormesh(x,y,g*abs2.(ψ0'))
xlabel("x/x0")
ylabel("y/x0")
colorbar()
subplot(1,2,2)
pcolormesh(x,y,angle.(ψ0'))
xlabel("x/x0")
ylabel("y/x0")
colorbar()

figure()
plot(x,g*abs2.(ψ0[:,300]))
#= TIME EVOLVE (real time) =#
plot(x,μ*(1-x.^2/Rx^2).*vortexcore(x-xv,ξ).^2)

function timeevolution2()
sinfo = Params()

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
  ωy = ωx
  ωz = 0.
#choice of time, length, energy units
  t0 = 1.0/ωy
  x0 = sqrt(ħ*t0/m)
  E0 = ħ/t0
#interactions
  #g  = (4*pi*ħ^2*as/m)*x0^3/E0 #dimensionless 3D
  g = 0.1 #test 2D
#damping parameters (dimensionless)
  γ = 0.0
  ℳ  = 0.0
#chemical potential (dimensionless)
  μ  = 10.0
#time evolution parameters
  ti = 0.0
  tf = 10.0  #evolve for 10 oscillator units
  Nt = 50
  t = collect(linspace(ti,tf,Nt))
  dt = 0.01π/μ #integrate step size [ - should have dt ≪ 2π/μ]
#== end template ==#

@pack sinfo = ωx,ωy,ωz,γ,ℳ,g,t0,x0,E0,μ,ti,tf,Nt,t,dt

#Initialize CField (dimensionless)
  basis = "Hermite"
  ecut = 30.5*ħ*ωy/E0
  Ω = [ωx; ωy]*t0
  cinfo = makecinfo(ecut,Ω,basis)
  tinfo = maketinfo(cinfo,4)
  @unpack en,P,M = cinfo
  #x,wx,Tx,y,wy,Ty = makealltrans(M,4,Ω)
  #W = wx.*wy'
#test transform
  #c0 set above for vortex
  ψ0  = c2x(c0,tinfo) #initial condition
  ψ = similar(ψ0)  #a field to write to in place

#equation of motion (in place)
  function nlin!(dc,c,tinfo)
      c2x!(ψ,c,tinfo)
      x2c!(dc,abs2.(ψ).*ψ,tinfo)
  end

  function Lgp!(dc,c,p,t)
      nlin!(dc,c,tinfo)
      dc .= P.*(-im*(1-im*γ)*((en - μ).*c .+ g*dc))
  end

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
sinfo2,cinfo2,sol2 = timeevolution2()

#plot solution for 2D
## Transform to cartesian grid
@unpack ħ,m,ωx,ωy,ωz,γ,ℳ,g,x0,t0,E0,μ,ti,tf,Nt,t,dt = sinfo2
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
tinfop = maketinfoplot(cinfo2,x,y)

ψ = c2x(sol2[10],tinfop)
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
#= TIME EVOLVE (real time) =#
plot(x,μ*(1-x.^2/Rx^2).*vortexcore(x-xv,ξ).^2)
