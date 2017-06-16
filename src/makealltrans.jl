function makealltrans(basis,N,Ω,n=4)
dim=length(N)
if dim==1
  Nx = N[1]
  ωx = Ω[1]
  #n-field transforms for PGPE
  x,wx,Tx = nfieldtrans(basis,Nx,n,ωx)
  return x,wx,Tx
elseif dim==2
  Nx,Ny = N
  ωx,ωy = Ω
  #n-field transforms for PGPE
  x,wx,Tx = nfieldtrans(basis,Nx,n,ωx)
  y,wy,Ty = nfieldtrans(basis,Ny,n,ωy)
  return x,wx,Tx,y,wy,Ty
elseif dim==3
  Nx,Ny,Nz = N
  ωx,ωy,ωz = Ω
  #n-field transforms for PGPE
  x,wx,Tx = nfieldtrans(basis,Nx,n,ωx)
  y,wy,Ty = nfieldtrans(basis,Ny,n,ωy)
  z,wz,Tz = nfieldtrans(basis,Nz,n,ωz)
  return x,wx,Tx,y,wy,Ty,z,wz,Tz
end
end
