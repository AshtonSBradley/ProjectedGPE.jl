"""

x,w,T = makeVtrans(M,basis)

Constructs transforms and associated arrays for precise numerical quadrature evalution
of the scattering potential Vϵ

[SPGPE](https://arxiv.org/abs/1507.02023)

# Arguments
 - `M` number of modes in the c-field.
 - `basis` is an optional string argument denoting the basis of eigenstates representing c-field state. Default is "hermite"

"""

function makeVtrans(M,ω)
x,w,T = nfieldtrans(4,M,basis)
