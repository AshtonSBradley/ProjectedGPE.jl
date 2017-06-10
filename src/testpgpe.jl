using ProjectedGPE
ecut = 50

ωx = 1
ωy = π
ωz = 5
Ω = [ωx,ωy,ωz]

#check 1D
cinfo, Espec, P = makecinfo("Hermite",ecut,[ωx,ωy])
#@step makecinfo("Hermite",ecut,[ωx,ωz])

#test 1D
Ω=[ωx]
M=[5]
x,wx,Tx=maketransinfo("Hermite",M,Ω)

#test 2D
Ω=[ωx, ωy]
M=[5,10]
x,wx,Tx,y,wy,Ty=maketransinfo("Hermite",M,Ω)

#test 2D
Ω=[ωx, ωy, ωz]
M=[25,17,33]
x,wx,Tx,y,wy,Ty,z,wz,Tz=maketransinfo("Hermite",M,Ω)
