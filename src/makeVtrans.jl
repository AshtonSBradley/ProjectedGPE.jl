"""

k,wk,Tkx = makeVtrans(basis,M,f)

Constructs transforms and associated arrays for precise numerical quadrature
evaluation of the scattering potential Vϵ

[SPGPE](https://arxiv.org/abs/1507.02023)

# Arguments
 - `M` number of modes in the c-field.
 - `basis` is an optional string argument denoting the basis of eigenstates representing c-field state. Default is "hermite"
 - `f` is the trap frequency relative to the chosen reference frequency
"""

function makeVtrans(basis,M,Nk=M,Na=M,f=1)
#4-field to x-space
x,wx,Tx = nfieldtrans(basis,M,4,f)

#4-field to k-space
k,wk,Tk = nfieldtrans(basis,Nk,4,1/f)

#4-field aux transforms to x
xa,wxa,Txa = nfieldtrans(basis,Na,4,f)

#4-field aux transform to k
ka, wka, Tka = nfieldtrans(basis,Na,4,1/f)
Tka = diagm((-im).^(0:2Na-1))*Tka

#Transform from x --> k via auxiliary states.
Txk = conj(Tka'*Txa)

return x,w,T,
end
