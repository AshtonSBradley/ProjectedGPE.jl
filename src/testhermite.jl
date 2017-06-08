using ProjectedGPE, ApproxFun, PyPlot, Combinatorics
M=100
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
f=Fun(Hermite(),c./sqrt(Hnorm));

x=linspace(-1.5sqrt(2M),1.5sqrt(2M),1000)
plot(x,abs(f.(x)).^2.*exp.(-x.^2))

x,w=gausshermite(M)

N2=sum(w.*abs(f.(x)).^2)

Float64(N2)

abs((N1-N2)/N1)


3*M*tol

#block to create T matrix
M=100
x,w=gausshermite(M)
T=zeros(x*ones(1,M))
n=0:BigFloat(M)-1;
Hnorm = sqrt(BigFloat(π))*BigFloat(2).^n.*factorial.(n)

for j = 0:M-1
  c=zeros(M)
  c[j+1]=1.
  Hn=Fun(Hermite(),c./sqrt(Hnorm));
  T[:,j+1]=Hn.(x).*exp(-x.^2/2)
end

x1,y1,T1=nfieldtrans("hermite",M,2)

mind = M
maximum(abs(T1[:,mind]-T[:,mind]))
