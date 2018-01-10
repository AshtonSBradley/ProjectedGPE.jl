#test fast transform via mapped Chebyschev
using ApproxFun, FastTransforms, BenchmarkTools

import FastTransforms: Butterfly

M=128

# Get the actual matrix from the plan in ApproxFun.plan_transform(::Hermite(),rand(256)), say.

PlanHerm = ApproxFun.plan_transform(Hermite(),rand(M));
x,w = PlanHerm.plan;

V=ApproxFun.hermitep(0:M-1,x)';
nrm=(V.^2)*w;
P= Diagonal(inv.(nrm))*V*Diagonal(w)

#P=Q[1:N,:]
BF = Butterfly(P, floor(Int, log2(size(P, 1))-6); sketch = :none, atol = 1e-14);

r = FastTransforms.allranks(BF) # just diagnostic info.

mean(r), std(r)
y = zeros(M)+im*zeros(M)
x = randn(M)+im*randn(M);
zr = zeros(M);
zi = zeros(M)
z = zeros(x)
@time A_mul_B!(y, P, x)

@time A_mul_B!(zr, BF, real(x)); A_mul_B!(zi, BF, imag(x));z=zr+im*zi

norm(y-z) # Did it work?


#test for transform to quadrature grid
N=256
x,wx,Tx = nfieldtrans("Hermite",N,2)

BFT = Butterfly(Tx, floor(Int, log2(size(Tx, 1))-6); sketch = :none, atol = 1e-32)

r = FastTransforms.allranks(BFT) # just diagnostic info.

mean(r), std(r)

M=N
x0 = randn(N)
x = randn(N)+im*randn(N)
x1 = randn(M)+im*randn(M)
y = randn(M)
zr = randn(N)
zi = randn(N)
z = randn(N)+im*randn(N)

@time A_mul_B!(x1, Tx, x)

@time A_mul_B!(zr, BFT, real(x)); A_mul_B!(zi, BFT, imag(x));z=zr+im*zi

norm(x-z) # Did it work?
