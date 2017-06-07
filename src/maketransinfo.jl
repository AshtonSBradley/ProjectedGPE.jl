function maketransinfo(basis,M,ecut,Ω)
dim=length(Ω)
  if dim==3
    ωx,ωy,ωz = Ω
    Mx,My,Mz = M
    #4-field transforms
    x,wx,Tx = nfieldtrans(basis,Mx,4,ωx)
    y,wy,Ty = nfieldtrans(basis,My,4,ωy)
    z,wz,Tz = nfieldtrans(basis,Mz,4,ωz)
  end

return TransInfo
end
