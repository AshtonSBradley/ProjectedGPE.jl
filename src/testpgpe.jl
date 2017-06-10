ecut = 20

ωx = 1
ωy = π
ωz = 5
Ω = [ωx,ωy,ωz]

#check 1D
cinfo, Espec, P = makecinfo("Hermite",ecut,[ωx])

#pack a 3d vector of frequencies?
@pack cinfo = Ω
