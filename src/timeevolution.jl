function timeevolution()
  #definition of the C-region

  #system params
  g = 0.1
  ωx = 2*pi
  γ = 0.07
  t0 = 1/ωx
  μ = 20

  #Time grid
  ti = 0.              #initial time
  tf = 80*t0          #final time
  Nt = 40             #size of time vector
  t  = collect(linspace(ti,tf,Nt))

  dt = 0.1*t0    #integrate step size [ - should have dt ≪ 2π/μ]

  #Initialize CField
  basis = "Hermite"
  ecut = 30ωx
  Ω = [ωx]
  cinfo = makecinfo(basis,ecut,Ω)
  @unpack en,P,M = cinfo; Mx=M[1]
  x,wx,Tx = maketransinfo(basis,M,Ω)

#test transform
  c0   = randn(Mx)+im*randn(Mx);
  ψ0   = Tx*c0

#PGPE nonlinearity 1D
#out of place
function nlin(c)
    ψ = Tx*c
    sum(Tx'*(wx.*abs(ψ).^2.*ψ),2)
end

#in place (mutating)
function nlin!(c,dc)
    #ψ = anisotrans(c,Tx)
    #ψ[:] = wx.*abs(ψ).^2.*ψ
    #dc[:] = sum(anisotrans(ψ,Tx'),2)
    ψ = Tx*c
    dc[:] = sum(Tx'*(wx.*abs(ψ).^2.*ψ),2)
end

#dPGPE in reservoir frame
#out of place
function Lgp(t,c)
 -im*(1-im*γ)*((en - μ).*c .+ g*nlin(c))
end

#in place
function Lgp!(t,c,dc)
    dc[:] = nlin!(c,dc)
    dc[:] .= -im*(1-im*γ)*((en - μ).*c .+ g*dc)
end

c0    = randn(Mx) + im*randn(Mx); #create new random initial state
tspan = (t[1],t[end])
prob = ODEProblem(Lgp!,c0,tspan)
alg = Tsit5()
#abstol = 1e-3!
#reltol = 1e-3
println("Started evolution ...")
@time sol = solve(prob,alg,dt=dt,saveat=t,alg_hints=[:stiff],save_everystep=false,dense=false);
println("... Finished.")
end
