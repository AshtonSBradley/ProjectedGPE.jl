"""

k,wk,Tkx = makescatteringtrans(basis,M,Nk,Na,ω)

Constructs transforms and associated arrays for precise numerical quadrature
evaluation of the scattering potential Vϵ

[SPGPE](https://arxiv.org/abs/1507.02023)

# Arguments
 - `basis` the basis of eigenstates representing c-field state. Default is "hermite"
 - `M` number of modes in the c-field.
 - `Nk` number of `k` points.
 - `Na` number of auxiliary oscillator states.
 - `ω` is the trap frequency relative to the chosen reference frequency.
"""

function makescatteringtrans(basis,M,Nk=M,Na=M,ω=1)
#4-field to x-space
x,wx,Tx = nfieldtrans(basis,M,4,ω)
#4-field to k-space
k,wk,Tk = nfieldtrans(basis,Nk,4,1/ω)
#4-field aux transforms to x
xa,wxa,Txa = nfieldtrans(basis,Na,4,ω)
#4-field aux transform to k
ka,wka,Tka = nfieldtrans(basis,Na,4,1/ω)
Tka = Tka*diagm((-im).^(0:2Na-1))
#Transform from x --> k via auxiliary states.
Txk = conj(Tka'*Txa)

return x,wx,Tx,k,wk,Txk
end
