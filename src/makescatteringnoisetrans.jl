function makescatteringnoisetrans(basis,M,Nk=M,ω=1)
#Noise transform: 2-field quadrature to go from basis to k-space
kn,wkn,Tkn = nfieldtrans(basis,Nk,3,1/ω)
Tkn = eigmat(basis,M,kn,1/ω)*diagm((-im).^(0:M-1))

#Noise: position space roots and weights for 3-field quadrature
xn,wxn,Txn = nfieldtrans(basis,M,3,ω)
Txkn = conj(Tkn'*Txn)
return xn,wxn,Txn,kn,wkn,Tkn
end
