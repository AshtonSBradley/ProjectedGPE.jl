using ProjectedGPE
ecut = 40
basis = "Hermite"

ωx = 1
ωy = π
ωz = 5
Ω = [ωx;ωy;ωz]

#check 1D
cinfo = makecinfo(ecut,Ω,basis)
#@step makecinfo("Hermite",ecut,[ωx,ωz])

#test 1D
Ω=[ωx]
M=[5]
x,wx,Tx=makealltrans(M,4,Ω,basis)

#test 2D
Ω=[ωx, ωy]
M=[5,10]
x,wx,Tx,y,wy,Ty=makealltrans(M,4,Ω,basis)

#test 3D
Ω=[ωx, ωy, ωz]
M=[25,17,33]
x,wx,Tx,y,wy,Ty,z,wz,Tz=makealltrans(M,4,Ω,basis)
