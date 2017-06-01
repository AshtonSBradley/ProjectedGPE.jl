"""

`sgpegamma(τ,μ,ω,ϵ)`

Returns the dimensionless damping rate in the simple-growth stochastic projected Gross-Pitaevskii equation.
Arguments are in units of ħω, where ω is reference trap frequency.

## Arguments
* `τ` - Reservoir temperature.
* `μ` - Reservoir chemical potential.
* `ω` - Frequency unit of simulation.
* `ϵ` - Cutoff energy above which modes comprise an incoherent reservoir.
[Simple growth SPGPE, PRA 77, 033616, 2008](https://arxiv.org/abs/1507.02023)

"""

function sgpegamma(τ,μ,ω,ϵ)
## bare rate and full cutoff dependent rates.
# requires lerch.m

#opts = F;
MHz     =1e6;
kHz     =1e3;
cm      =1e-2;
mm      =1e-3;
mum     =1e-6;
nm      =1e-9;
muK     =1e-6;
nK      =1e-9;
kB      =1.38e-23;
ħ       =1.034e-34;

amu     = 1.6605402e-27;
m       = 87*amu; #Rb87
bohR    = 5.29177249e-11;
as      = 100*bohR;
T       = τ*ħ*ω/kB;         #T in Kelvin

γ0      = 4*m/π/ħ^3*(as*kB*T)^2;   #All with dimensions

μ       = μ*ħ*ω;            #make μ, Ecut have dimensions
ϵ       = ϵ*hbar*w0;

μ       = μ/kB/T;          #in units of kB*T (for exp(beta...))
ϵ       = ϵ/kB/T;
G1i     = log(1-exp(μ-ϵ))^2;

Nmax2 = 1000;
NmaxLerch = 1000;
R = 1:Nmax2;
temp =  exp(2*(μ-ϵ))*exp(R*(μ-2*ϵ)).*lerch(exp(μ-ϵ),1,R+1,NmaxLerch).^2;
G2i = sum(temp);

Gi = G1i + G2i;

g = γ0*Gi;

#Rates in computational units:
gbareComp  = γ0*ħ/kB/T;
gComp      = g*ħ/kB/T;

#change sinfo.gamma here:

#sinfo.gamma =
G = gComp;
end
