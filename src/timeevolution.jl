#function timeevolution()
  #definition of the C-region
  #M = 30
  #x,wx,Tx=nfieldtrans("Hermite",M,4)

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
  @unpack Espec,P,N = cinfo
  x,wx,Tx = maketransinfo(basis,N,Ω)
  
#test transform
  c0   = randn(M)+im*randn(M);
  ψ0   = Tx*c0

#PGPE nonlinearity 1D
#out of place
function nlin(c)
    ψ = anisotrans(c,Tx)
    sum(Tx'*(wx.*abs(ψ).^2.*ψ),2)
end

#in place (mutating)
function nlin!(c,dc)
    ψ = anisotrans(c,Tx)
    dc[:] = sum(Tx'*(wx.*abs(ψ).^2.*ψ),2)
end

#dPGPE in reservoir frame
#out of place
function Lgp(t,c)
 -im*(1-im*γ)*((Espec - μ).*c .+ g*nlin(c))
end

#in place
function Lgp!(t,c,dc)
    dc[:] = nlin!(c,dc)
    dc[:] .= -im*(1-im*γ)*((Espec - μ).*c .+ g*dc)
end

c0    = randn(M) + im*randn(M); #create new random initial state
tspan = (t[1],t[end])
prob = ODEProblem(Lgp!,c0,tspan)
alg = Tsit5()
#abstol = 1e-3!
#reltol = 1e-3
println("Started evolution...")
sol = solve(prob,alg,dt=dt,saveat=t,alg_hints=[:stiff],save_everystep=false,dense=false);
println("... Finished.")
#end
