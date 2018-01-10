"""
    k, wk, Tkx = makescatteringtrans(M,basis,Nk,Na,ω)

Construct transforms and associated arrays for precise numerical quadrature
evaluation of the scattering potential ``V_\\epsilon``

 `M`: number of modes in the c-field.

 `ω`: trap frequency relative to the chosen reference frequency.

 `basis`: name of orthonormal eigenbasis representing c-field state. Currently `"Hermite"` is implemented.

 `n`: field product that will be integrated to machine precision.

 `Nk`: number of `k` points.

 `Na`: number of auxiliary oscillator states.

 External links

 [SPGPE](https://arxiv.org/abs/1507.02023)

"""

function makescatteringtrans(M,ω=1,basis="Hermite",n=2,Nk=M,Na=M)
#4-field to x-space
x,wx,Tx = nfieldtrans(M,n,ω,basis)
#4-field to k-space
k,wk,Tk = nfieldtrans(M,n,1/ω,basis)
#4-field aux transforms to x
xa,wxa,Txa = nfieldtrans(Na,n,ω,basis)
#4-field aux transform to k
ka,wka,Tka = nfieldtrans(Na,n,1/ω,basis)
Tka = Tka*diagm((-im).^(0:Na-1))
#Transform from x --> k via auxiliary states.
Txk = conj(Tka'*Txa)

return x,wx,Tx,k,wk,Txk
end
