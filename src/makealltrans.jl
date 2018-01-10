function makealltrans(N,Ω,basis="Hermite",n=4)
dim=length(N)
if dim==1
  Nx = N[1]
  ωx = Ω[1]
  #n-field transforms for PGPE
  x,wx,Tx = nfieldtrans(Nx,n,ωx,basis)
  return x,wx,Tx
elseif dim==2
  Nx,Ny = N
  ωx,ωy = Ω
  #n-field transforms for PGPE
  x,wx,Tx = nfieldtrans(Nx,n,ωx,basis)
  y,wy,Ty = nfieldtrans(Ny,n,ωy,basis)
  return x,wx,Tx,y,wy,Ty
elseif dim==3
  Nx,Ny,Nz = N
  ωx,ωy,ωz = Ω
  #n-field transforms for PGPE
  x,wx,Tx = nfieldtrans(Nx,n,ωx,basis)
  y,wy,Ty = nfieldtrans(Ny,n,ωy,basis)
  z,wz,Tz = nfieldtrans(Nz,n,ωz,basis)
  return x,wx,Tx,y,wy,Ty,z,wz,Tz
end
end
