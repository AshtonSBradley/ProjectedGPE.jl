function evolve2()
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
  Γ̄  = 0.1
  M̄  = 0.0
#chemical potential (dimensionless)
  μ  = 12.0
#time evolution parameters
  ti = 0.0
  tf = 2.0/Γ̄  #evolve for 2 damping times
  Nt = 50
  t = collect(linspace(ti,tf,Nt))
  dt = 0.01π/μ #integrate step size [ - should have dt ≪ 2π/μ]
#== end template ==#

  @pack siminfo = ωx,ωy,ωz,Γ̄,M̄,g,t0,E0,x0,μ,ti,tf,Nt,t,dt

  #Initialize CField
  basis = "Hermite"
  ecut = 30 #units of ħ*ωy
  Ω = [ωx; ωy]*t0
  cinfo = makecinfo(ecut,Ω,basis)
  @unpack en,P,M = cinfo ;Mx=M[1];My=M[2]
  x,wx,Tx,y,wy,Ty = makealltrans(M,Ω,n=4,basis=basis)
  Wxy = wx.*wy'
#test transform
  c0   = randn(Mx,My)+im*randn(Mx,My);
  c0   = P.*c0;   #Project
  ψ0   = Tx*c0*Ty'
  ψ    = Tx*c0*Ty' #a field to write to in place if needed

#PGPE nonlinearity 1D
#out of place
function nlin(c)
    ψ = Tx*c*Ty'
    Tx'*(Wxy.*abs.(ψ).^2.*ψ)*Ty
end

function nlin!(c,dc)
    ψ = Tx*c*Ty'
    dc.= Tx'*(Wxy.*abs.(ψ).^2.*ψ)*Ty
end

#dPGPE in reservoir frame
#out of place
function Lgp(t,c)
 -im*(1-im*Γ̄)*((en - μ).*c .+ g*nlin(c))
end

#in place
function Lgp!(t,c,dc)
    dc = nlin!(c,dc)
    dc.= -im*(1-im*Γ̄)*((en - μ).*c .+ g*dc)
    dc.= P.*dc
end

c0    = randn(Mx,My) + im*randn(Mx,My); #create new random initial state
tspan = (t[1],t[end])
prob = ODEProblem(Lgp!,c0,tspan)
alg = Tsit5()
#abstol = 1e-3!
#reltol = 1e-3
println("Started evolution ...")
@time sol = solve(prob,alg,dt=dt,saveat=t,alg_hints=[:stiff],save_everystep=false,dense=false);
println("... Finished.")
return siminfo,cinfo,sol
end
