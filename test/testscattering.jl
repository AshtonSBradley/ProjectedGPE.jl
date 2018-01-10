using ProjectedGPE
ecut = 50
basis = "Hermite"

ωx = 1
ωy = π
ωz = 5
Ω = [ωx;ωy;ωz]

cinfo = makecinfo(ecut,Ω,basis)

#test 3D
Ω=[ωx, ωy, ωz]
M=[25,17,33]
x,wx,Tx,y,wy,Ty,z,wz,Tz=makealltrans(M,Ω,basis)
