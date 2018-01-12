function makescatteringnoisetrans(M,ω=1,basis="Hermite",Nk=M)
#Noise transform: 2-field quadrature to go from basis to k-space
kn,wkn,Tkn = nfieldtrans(Nk,3,ω=1/ω,basis=basis)
Tkn = eigmat(M,kn,ω=1/ω,basis=basis)*diagm((-im).^(0:M-1))

#Noise: position space roots and weights for 3-field quadrature
xn,wxn,Txn = nfieldtrans(M,3,ω=ω,basis=basis)
Txkn = conj(Tkn'*Txn)
return xn,wxn,Txn,kn,wkn,Tkn
end
