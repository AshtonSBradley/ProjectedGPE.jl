M=30
c=randn(M)+im*randn(M)
N=sum(abs(c).^2)

using ApproxFun, PyPlot
f=Fun(Hermite(),c);

#put onto specific x grid 
S=Hermite(-10..10);
x=points(S,100);
plot(x,abs(f.(x)).^2.*exp.(-x.^2))
