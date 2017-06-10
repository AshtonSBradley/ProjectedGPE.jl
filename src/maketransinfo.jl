function maketransinfo(basis,M,Ω)
dim=length(Ω)
  if dim==3
    ωx,ωy,ωz = Ω
    Mx,My,Mz = M
    #4-field transforms for PGPE
    x,wx,Tx = nfieldtrans(basis,Mx,4,ωx)
    y,wy,Ty = nfieldtrans(basis,My,4,ωy)
    z,wz,Tz = nfieldtrans(basis,Mz,4,ωz)
  end

return x,wx,Tx,y,wy,Ty,z,wz,Tz
end
