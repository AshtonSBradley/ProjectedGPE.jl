#test new formulation of algorithm, for adaptive sde solving

using ProjectedGPE, DifferentialEquations

using PyPlot, Interact

import ProjectedGPE.Lgp!
#definition of the C-region
M = 30

#system params
g = 0.1
ω = 1.0
γ = 0.1
t0 = 1/ω
μ = 25.0

#Time grid
ti = 0.              #initial time
tf = 80*t0         #final time
Nt = 50             #size of time vector
t  = collect(linspace(ti,tf,Nt))

dt = 0.01*t0    #initial step size [ - should have dt ≪ 2π/μ]

function makexktrans2(basis,M,ω=1.0)
#n-field to x-space
#x,wx,Tx = nfieldtrans(basis,M,nx,ω)
x,wx = gausshermite(2M) #rule order should be set based on terms in eq
wx = wx.*exp.(x.^2)/√ω
Tx = eigmat("Hermite",M,x/√ω,ω) #take grid to spectral
#n-field to k-space
k,wk = gausshermite(2M) #set rule order by terms in eq
wk = wk.*exp.(k.^2)*√ω
Tk = eigmat("Hermite",M,k*√ω,1/ω)
Tk = Tk*diagm((-im).^(0:M-1)) #take grid to spectral
#Transform from x --> k via auxiliary states.
Txk = Tk*Tx'
    return x,wx,Tx,k,wk,Tk,Txk
end

x,wx,Tx,k,wk,Tk,Txk = makexktrans2("Hermite",M)

#M modes on 4-field grid, gives 2M points
x4,wx4,Tx4 = nfieldtrans("Hermite",M,4)

Tx24 = Tx4*(Tx'.*wx')

#PGPE kinetic term
#out of place
function kinetic(ψ)
    return Txk'*(0.5*wk.*k.^2.*(Txk*(wx.*ψ)))
end

#in place (mutating)
function kinetic!(ψ,dψ)
    dψ .= Txk'*(0.5*wk.*k.^2.*(Txk*(wx.*ψ)))
end

#ψ0=randn(size(x))+im*randn(size(x)) #create new random initial state
c0=zeros(M)+im*zeros(M) #zero pad to length 2M
c0[25]=sqrt(100) # initial state
ψ0=Tx*c0
ϕ = Tx24*ψ0

function nlin(ψ)
    #Tx*(P.*Tx'*(wx.*abs.(ψ).^2.*ψ))  #right weight??
    ϕ = Tx24*ψ
    return Tx24'*(wx4.*abs.(ϕ).^2.*ϕ)
end

#Time evolution

#dPGPE in reservoir frame (μ is energy zero)

#out of place
function Lgp(t,ψ)
    #return -im*(1-im*γ)*((g*abs.(ψ).^2 .+ 0.25.*x.^2.*exp(-.5*x.^2) .- μ).*ψ .+ kinetic(ψ))
    ψ = -im*(1-im*γ)*((g*abs.(ψ).^2 .- μ .+ 0.5*(x/sqrt(2)).^2).*ψ .+ kinetic(ψ)) #4 field to 2 field
    #ψ = Tx*(Tx'*(wx.*ψ)) #project
    #ψ = wx.*ψ #project
end

#in place
function Lgp!(t,ψ,dψ)
    kinetic!(ψ,dψ)
    #dψ[:] = -im*(1-im*γ)*((g*abs.(ψ).^2  .+ 0.5.*x.^2 .- μ).*ψ .+ dψ) #2 field to 2 field (non-linear works...?)
    #maybe we just need a function for nonlinear term specific to grids (aswell as potential)
    #result: at least nonlin will be a quadrature
    #get all quaratures right:
    #4 field n for P.E., K.E.
    #same trans, different weight for nlin (?doesn't work if make different):
    dψ .= -im*(1-im*γ)*((0.5.*x.^2 .- μ).*ψ .+ g*nlin(ψ) .+ dψ)
end



tspan = (t[1],t[end])
prob = ODEProblem(Lgp!,ψ0,tspan,callback=CallbackSet(),mass_matrix=I)
alg = Tsit5()
#alg = DP5()
#abstol = 1e-12
#reltol = 1e-12
println("Started ...")
@time sol = solve(prob,alg,dt=dt,saveat=t,alg_hints=[:stiff],save_everystep=false,dense=false);
println("... Finished.")


## Transform to cartesian grid for plotting

R = sqrt(2μ)
xMax = 2*R
Nx = 200
xp = collect(linspace(-xMax,xMax,Nx))
Txp = eigmat("Hermite",M,xp)
Ntot=zeros(length(t))
for i=1:length(t)
    Ntot[i]=sum(abs.(Tx'*(wx.*sol[i])).^2)
end
figure(figsize=(6,2))
plot(t,Ntot./Ntot[1]);#ylim(0,1.2)


#Plot
f=figure(figsize=(9,2))
@manipulate for i=1:length(t) withfig(f,clear=true) do
        ψ = Txp*Tx'*(wx.*sol[i]);
        θ = unwrap(angle(ψ));
        subplot(1,2,1);plot(xp/R,g*abs(ψ).^2,xp/R,ones(xp)*μ,"g:");
        xlim(-xMax/R,xMax/R);#ylim(-1,1.5*μ);
        xlabel(L"x/R");ylabel(L"g|ψ|^2")
        subplot(1,2,2);plot(xp/R,θ);
        xlim(-xMax/R,xMax/R);#ylim(-pi*1.1,pi*1.1);
        xlabel(L"x/R");ylabel(L"angle(ψ)");
    end
end
