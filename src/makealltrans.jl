function makealltrans(N,Ω;n=4,basis="Hermite")
dim=length(N)
if dim==1
  Nx = N[1]
  ωx = Ω[1]
  #n-field transforms for PGPE
  x,wx,Tx = nfieldtrans(Nx,ωx,n=n,basis=basis)
  return x,wx,Tx
elseif dim==2
  Nx,Ny = N
  ωx,ωy = Ω
  #n-field transforms for PGPE
  x,wx,Tx = nfieldtrans(Nx,ωx,n=n,basis=basis)
  y,wy,Ty = nfieldtrans(Ny,ωy,n=n,basis=basis)
  return x,wx,Tx,y,wy,Ty
elseif dim==3
  Nx,Ny,Nz = N
  ωx,ωy,ωz = Ω
  #n-field transforms for PGPE
  x,wx,Tx = nfieldtrans(Nx,ωx,n=n,basis=basis)
  y,wy,Ty = nfieldtrans(Ny,ωy,n=n,basis=basis)
  z,wz,Tz = nfieldtrans(Nz,ωz,n=n,basis=basis)
  return x,wx,Tx,y,wy,Ty,z,wz,Tz
end
end
