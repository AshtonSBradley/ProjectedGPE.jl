function timeevolution2()
  #definition of the C-region
  siminfo = SimParams()
  #system params

  ωx = 2π
  ωy = 4ωx
  ωz = 0.0
  t0 = 1/ωy  #choose time unit as largest ω
  g = 0.1
  Γ̄ = 0.07
  M̄ = 0.0
  μ  = 20.0

  #Time grid
  ti = 0.0              #initial time
  tf = 80*t0          #final time
  Nt = 40             #size of time vector
  t  = collect(linspace(ti,tf,Nt))

  dt = 0.1*t0    #integrate step size [ - should have dt ≪ 2π/μ]

  @pack siminfo = ωx,ωy,ωz,Γ̄,M̄,g,t0,μ,ti,tf,Nt,t,dt

  #Initialize CField
  basis = "Hermite"
  ecut = 30 #units of ħ*ωy
  Ω = [ωx; ωy]*t0
  cinfo = makecinfo(basis,ecut,Ω)
  @unpack en,P,M = cinfo ;Mx=M[1];My=M[2]
  x,wx,Tx,y,wy,Ty = maketransinfo(basis,M,Ω)
  Wxy = wx.*wy'
#test transform
  c0   = randn(Mx,My)+im*randn(Mx,My);
  ψ0   = Tx*c0*Ty'
  ψ    = Tx*c0*Ty' #a field to write to in place if needed

#PGPE nonlinearity 1D
#out of place
function nlin(c)
    ψ = Tx*c*Ty'
    Tx'*(Wxy.*abs(ψ).^2.*ψ)*Ty
end

function nlin!(c,dc)
    ψ = Tx*c*Ty'
    dc[:,:] = Tx'*(Wxy.*abs(ψ).^2.*ψ)*Ty
end

#dPGPE in reservoir frame
#out of place
function Lgp(t,c)
 -im*(1-im*Γ̄)*((en - μ).*c .+ g*nlin(c))
end

#in place
function Lgp!(t,c,dc)
    dc[:] = nlin!(c,dc)
    dc[:] = -im*(1-im*Γ̄)*((en - μ).*c .+ g*dc)
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
