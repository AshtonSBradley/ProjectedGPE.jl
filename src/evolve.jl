function evolve()
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
  γ  = 0.1
  ℳ  = 0.0
#chemical potential (dimensionless)
  μ  = 12.0
#time evolution parameters
  ti = 0.0
  tf = 1.0/γ  #evolve for 2 damping times
  Nt = 50
  t = collect(linspace(ti,tf,Nt))
  dt = 0.01π/μ #integrate step size [ - should have dt ≪ 2π/μ]
#== end template ==#

  @pack sinfo = ωx,ωy,ωz,γ,ℳ,g,t0,E0,x0,μ,ti,tf,Nt,t,dt

  #Initialize CField
  basis = "Hermite"
  ecut = 30 #units of ħ*ωy
  Ω = [ωx; ωy]*t0
  cinfo = makecinfo(ecut,Ω,basis)
  tinfo = maketinfo(cinfo,4)
  @unpack en, P, M = cinfo; Mx, My = M;
  #x,wx,Tx,y,wy,Ty = makealltrans(M,Ω,n=4,basis=basis)
  #W = wx.*wy'
#test transform
  c0   = randn(Mx,My)+im*randn(Mx,My);
  c0   = P.*c0;   #Project
  ψ0   = c2x(c0,tinfo)
  ψ    = simlar(ψ0) #a field to write to in place if needed



c0    = randn(Mx,My) + im*randn(Mx,My); #create new random initial state
tspan = (t[1],t[end])
prob = ODEProblem(Lgp!,c0,tspan)
alg = DP5()
#abstol = 1e-3!
#reltol = 1e-3
println("Started evolution ...")
@time sol = solve(prob,alg,dt=dt,saveat=t);
println("... Finished.")
return sinfo,cinfo,sol
end
