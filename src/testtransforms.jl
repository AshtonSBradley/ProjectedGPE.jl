using ProjectedGPE, ApproxFun, FastGaussQuadrature, PyPlot, Combinatorics
M=350
tol = 1e-14
#c=Vector{BigFloat}(randn(M))+im*Vector{BigFloat}(randn(M))
#c=randn(M)+im*randn(M)
c=zeros(M)
c[end]=1.
N1=sum(abs(c).^2)

#evaluate modes using BigFloats
#can then use ApproxFun.jl for polynomials, applying correct normalization.
#Combinatorics.jl needed to evaluate factorial of Vector.

n=0:BigFloat(M)-1;
Hnorm = sqrt(BigFloat(π))*BigFloat(2).^n.*factorial.(n)

#Inserting norm.
f=Fun(Hermite(),c./sqrt(Hnorm))
f2=Fun(GaussWeight(Hermite()),c./sqrt(Hnorm))
x=linspace(-1.5sqrt(2M),1.5sqrt(2M),1000)
plot(x,abs(f.(x)).^2.*exp.(-x.^2))

x,w=gausshermite(M)

N2=sum(w.*abs(f.(x)).^2)
N3=sum(w.*abs(f2.(x).*exp(x.^2)).^2)

Float64(N2)

abs((N1-N2)/N1)


3*M*tol

#block to create T matrix
M=10
x,w=gausshermite(M)
T=zeros(x*ones(1,M))
n=0:BigFloat(M)-1;
Hnorm = sqrt(BigFloat(π))*BigFloat(2).^n.*factorial.(n)

for j = 0:M-1
  c=zeros(M)
  c[j+1]=1.
  #Hn=Fun(Hermite(),c./sqrt(Hnorm));
  #T[:,j+1]=Hn.(x).*exp(-x.^2/2)
  Hn=Fun(GaussWeight(Hermite()),c)
  T[:,j+1]=Hn.(x).*exp(x.^2/2)
end

x1,y1,T1=nfieldtrans("Hermite",M,2,2)

mind = M
maximum(abs(T1[:,mind]-T[:,mind]))

T2=eigmat("Hermite",M,x,2)
maximum(abs(T1[:,mind]-T2[:,mind]))

#-------------------------------
#Make and test Laguerre transform and quadrature
M=27
α = 0.
ω = 1.
x,w=gausslaguerre(M,α)
T=zeros(x*ones(1,M))
n=0:M-1;
Lnorm = exp.(lgamma.(n+α+1)-lgamma.(n+1))/ω

for j = 0:M-1
  c=zeros(M)
  c[j+1]=1.
  Lna=Fun(Laguerre(α),c./sqrt(Lnorm));
  T[:,j+1]=Lna.(sqrt(ω)*x).*exp(-sqrt(ω)*x/2).*(sqrt(ω)*x).^(α/2)
end

c=zeros(M);
c[end]=BigFloat(1.)

#Gauss Laguerre quadrature can absorb the factor x^α/2 using associated Laguerre
#c=randn(M)+im*rand(M)
N1 = sum(abs(c).^2)
L=Fun(Laguerre(α),c./sqrt(Lnorm))




L2=Fun(WeightedLaguerre(α),c./sqrt(Lnorm))
N2=sum(w.*abs(L.(sqrt(ω)*x)).^2)



N3=sum(w.*abs(L2.(sqrt(ω)*x).*exp(sqrt(ω)*x)./(sqrt(ω)x).^(α/2)).^2)
x=linspace(0,20sqrt(M),1000)
plot(x,abs(L.(x)).^2.*exp.(-x).*x.^α)

Lnorm

#------------------------------------
#test Laguerre
M=27
α = 0.
ω = 1.
#x,w=gausslaguerre(M,α)

#T = eigmat("Laguerre",M,x/ω,ω,α)

x,w,T=nfieldtrans("Laguerre",M,2,ω,α)
#norm test

c=zeros(M);
c[end]=1.
N1 = sum(abs(c).^2)

Ln = T*c


N2=sum(w.*abs(Ln).^2)
f1 = Fun(GaussWeight(Hermite()),[1])

#------------------------------------
#test Hermite
M=21

ω = 10
#x,w=gausslaguerre(M,α)

#T = eigmat("Laguerre",M,x/ω,ω,α)

x,w,T=nfieldtrans("Hermite",M,2,ω)
#norm test

c=zeros(M);
c[end]=1.
N1 = sum(abs(c).^2)

Hn = T*c


N2=sum(w.*abs(Hn).^2)
