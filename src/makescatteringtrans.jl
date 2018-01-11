"""
    k,wk,Tkx = makescatteringtrans(M;n=4,ω=1.0,basis="Hermite",Nk=M,Na=M)

Construct transforms and associated arrays for precise numerical quadrature
evaluation of the scattering potential ``V_\\epsilon(x,t)`` of the energy-damped SPGPE.

 `M`: number of modes in the c-field.

 `n`: order of field product that will be integrated to machine precision.

 `ω`: trap frequency relative to the chosen reference frequency.

 `basis`: name of orthonormal eigenbasis representing c-field state. Currently `"Hermite"` is implemented.

 `Nk`: number of `k`-space points.

 `Na`: number of auxiliary oscillator states.

 External links

 [Numerical method for the stochastic projected Gross-Pitaevskii equation, Physical Review E 89, 013302 (2014)](https://arxiv.org/abs/1310.0161)

"""

function makescatteringtrans(M;n=4,ω=1.0,basis="Hermite",Nk=M,Na=M)

    #4-field to x-space
    x,wx,Tx = nfieldtrans(M,n,ω=ω,basis=basis)

    #4-field to k-space. To include Nk need low leavel call to eigmat
    k,wk,Tk = nfieldtrans(M,n,ω=1/ω,basis=basis)

    #4-field aux transforms to x
    xa,wax,Tax = nfieldtrans(Na,n,ω=ω,basis=basis)

    #4-field aux transform to k
    ka,wak,Tak = nfieldtrans(Na,n,ω=1/ω,basis=basis)
    Tak = Tak*diagm((-im).^(0:Na-1))

    #Transform from x --> k via auxiliary states.
    Txk = conj(Tak*(Tax'*diagm(wax)))

return x,wx,Tx,k,wk,Txk
end
