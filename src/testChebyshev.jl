#test fast transform via mapped Chebyschev
using ApproxFun, FastTransforms, BenchmarkTools

import FastTransforms: Butterfly

M=128
N=Int(M/2)
# Get the actual matrix from the plan in ApproxFun.plan_transform(::Hermite(),rand(256)), say.
#PlanHerm = ApproxFun.plan_transform(Hermite(),rand(N))
#Px,Py=PlanHerm.plan
#P=Px*Py'

PlanHerm = ApproxFun.plan_transform(Hermite(),rand(M));
x,w = PlanHerm.plan;

V=ApproxFun.hermitep(0:M-1,x)';
nrm=(V.^2)*w;
Q = Diagonal(inv.(nrm))*V*Diagonal(w)

P=Q[1:N,:]
BF = Butterfly(P, floor(Int, log2(size(P, 1))-6); sketch = :none, atol = 1e-14);

r = FastTransforms.allranks(BF) # just diagnostic info.

mean(r), std(r)

x0 = randn(M)
x = randn(M)+im*randn(M);
y = zeros(x[1:N]);
zr = zeros(N);
zi = zeros(N)
z=zeros(x)
@time A_mul_B!(y, P, x)

@time A_mul_B!(zr, BF, real(x)); A_mul_B!(zi, BF, imag(x));#z=zr+im*zi

norm(y-z) # Did it work?


#test for transform to quadrature grid
N=128
x,wx,Tx = maketransinfo("Hermite",128,[1])

BFT = Butterfly(Tx, floor(Int, log2(size(Tx, 1))-6); sketch = :none, atol = 1e-24);

r = FastTransforms.allranks(BFT) # just diagnostic info.

mean(r), std(r)

x0 = randn(N)
x = randn(N)+im*randn(N);
y = zero(x);
zr = zero(x0);
zi = zero(x0)
z=zero(x)
@time A_mul_B!(y, Tx, x)

@time A_mul_B!(zr, BFT, real(x)); A_mul_B!(zi, BFT, imag(x));#z=zr+im*zi

norm(y-z) # Did it work?
