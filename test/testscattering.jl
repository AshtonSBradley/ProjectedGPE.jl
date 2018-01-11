using Revise, ProjectedGPE, PyPlot
ecut = 20
basis = "Hermite"

ωx = 1
ωy = π
ωz = 5

#test 3D
Ω=[ωx, ωy, ωz]
M=[25,17,33]
x,wx,Tx,y,wy,Ty,z,wz,Tz=makealltrans(M,Ω,basis)
x,wx,Tx,k,wk,Txk = makescatteringtrans(ecut,4,1.0,basis)

k = collect(linspace(0,10,200))
plot(k,scatteringkernel(k))

σ = 1.0

f1(k) = @. erfcx(abs(k)*σ/sqrt(2))/sqrt(8*π*σ^2)
f2(k) = @. exp(abs(k*σ/2)^2)*besselk(0,abs(k*σ/2)^2)/(2*π)
f3(k) = @. 1/abs(k)

semilogy(k,f1(k))
semilogy(k,f2(k))
semilogy(k,f3(k))
xlabel("k")
legend(["1D","2D","3D"])

M=20
n=4
ω=1.
basis="Hermite"
Na=43

x,wx,Tx,k,wk,Txk=makescatteringtrans(20,n=2)

c=randn(M)+im*randn(M)
ψ=Tx*c
ϕ=Txk*ψ
ψinv=Txk'*ϕ
