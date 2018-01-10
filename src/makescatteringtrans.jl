"""

k,wk,Tkx = makescatteringtrans(basis,M,Nk,Na,ω)

Constructs transforms and associated arrays for precise numerical quadrature
evaluation of the scattering potential Vϵ

[SPGPE](https://arxiv.org/abs/1507.02023)

 `basis` the basis of eigenstates representing c-field state. Default is "Hermite".

 `M` number of modes in the c-field.

 `ω` is the trap frequency relative to the chosen reference frequency.

 `n` is the field product that will be integrated exactly.

 `Nk` number of `k` points.

 `Na` number of auxiliary oscillator states.

"""

function makescatteringtrans(basis,M,ω=1,n=2,Nk=M,Na=M)
#4-field to x-space
x,wx,Tx = nfieldtrans(basis,M,n,ω)
#4-field to k-space
k,wk,Tk = nfieldtrans(basis,M,n,1/ω)
#4-field aux transforms to x
xa,wxa,Txa = nfieldtrans(basis,Na,n,ω)
#4-field aux transform to k
ka,wka,Tka = nfieldtrans(basis,Na,n,1/ω)
Tka = Tka*diagm((-im).^(0:Na-1))
#Transform from x --> k via auxiliary states.
Txk = conj(Tka'*Txa)

return x,wx,Tx,k,wk,Txk
end
