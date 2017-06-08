using ProjectedGPE, ApproxFun, PyPlot, Combinatorics
M=50
tol = 1e-14
#c=Vector{BigFloat}(randn(M))+im*Vector{BigFloat}(randn(M))
#c=randn(M)+im*randn(M)
c=zeros(M)
c[end]=1.
N1=sum(abs(c).^2)
n=0:BigFloat(M)-1;
Hnorm = sqrt(BigFloat(π))*BigFloat(2).^n.*factorial.(n)

#use ApproxFun to evaluate hermite polynomials, inserting norm.
f=Fun(Hermite(),c./sqrt(Hnorm));

#put onto specific x grid
#S=Hermite(-10..10);
#x=points(S,100);

x=linspace(-1.5sqrt(2M),1.5sqrt(2M),500)
plot(x,abs(f.(x)).^2.*exp.(-x.^2))

x,w=gausshermite(M)

N2=sum(w.*abs(f.(x)).^2)

abs((N1-N2)/N1)

3*M*tol

typeof(f)

g=Fun(Fourier())
