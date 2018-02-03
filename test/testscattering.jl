using Revise, ProjectedGPE, PyPlot
ecut = 20
basis = "Hermite"

ωx = 1
ωy = π
ωz = 5

#test 3D
Ω=[ωx, ωy, ωz]
M=[25,17,33]
x,wx,Tx,y,wy,Ty,z,wz,Tz = makealltrans(M,4,Ω,basis)
x,wx,Tx,k,wk,Txk = makescatteringtrans(ecut,n=4)

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

#test derivatives and energies
M=20
n=4
ω=2.             #test freq. dependence
basis="Hermite"
Na=43

x,wx,Tx,k,wk,Txk = makescatteringtrans(20,n=2)

c=randn(20)+im*randn(20)
ψ=Tx*c
ϕ=Txk*(wx.*ψ)
ψinv=Txk'*(wk.*ϕ)
#uniterity test passes.

#test derivatives for 2-field product
c=complex(zeros(20))
c[1] = 1.0
ψ=Tx*c
ϕ=Txk*(wx.*ψ)
ψx = Txk'*(wk.*(im*k.*ϕ))
ψxx= Txk'*(wk.*(-k.^2).*ϕ)
V = 0.5*sum(Tx'*(wx.*x.^2.*ψ).*conj(c))
K = -0.5*sum(Tx'*(wx.*ψxx).*conj(c))

E = V+K

figure()
plot(x,abs2.(ψ))
plot(x,real.(ψx))
plot(x,imag.(ψx))

#test 4-field term ∇⋅j
#3-field transform to k space
x,wx,Tx,k,wk,Txk = makescatteringtrans(20,n=3)

c=complex(zeros(20))
c[1] = 1.0
ψ = Tx*c
ϕ = Txk*(wx.*ψ)
ψx = Txk'*(wk.*(im*k.*ϕ))
ψxx= Txk'*(wk.*(-k.^2).*ϕ)

divj = im*0.5(ψ.*conj(ψxx)-conj(ψ).*ψxx)
divjk = Txk*(wx.*divj)
# then multiply by scattering scatteringkernel
divjkS = f1(k).*divjk
# then invert
divjkSx = Txk'*(wk.*divjkS)
# then × by ψ
Vψ = divjkSx.*ψ
# then project back to c-space to get term in EOM
dψᵥ = Tx'*(wx.*Vψ)


figure()
plot(x,abs2.(ψ))
plot(x,real.(ψx))
plot(x,real.(divj))
