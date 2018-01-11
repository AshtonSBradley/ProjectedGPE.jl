function makespecops(ecut,Ω,basis)
if basis=="Hermite"
  dim=length(Ω)
  if dim==1
    ωx = Ω[1]
    e0 = 0.5*ωx
    ecut < e0 && error("ecut must exceed the zero point energy.")
    Mx,nx,en = nenergy(ecut,e0,ωx,basis)
    P  = en .< ecut
    en = P.*en
    ax = sqrt(1/ωx) #in dimensionless units
    X,Px = ladderops(Mx,ax)
    return P,en,X,Px
  elseif dim==2
    ωx,ωy = Ω
    e0 = 0.5(ωx+ωy)
    ecut < e0 && error("ecut must exceed the zero point energy.")
    Mx,nx,ex = nenergy(ecut,e0,ωx,basis)
    My,ny,ey = nenergy(ecut,e0,ωy,basis)
    en = [ex[i+1]+ey[j+1] for i in nx, j in ny]
    P  = en .< ecut
    en = P.*en
    ax = sqrt(1/ωx); ay = sqrt(1/ωy)
    X,Px = ladderops(Mx,ax)
    Y,Py = ladderops(My,ay)
    return P,en,X,Px,Y,Py
  elseif dim==3
    ωx,ωy,ωz = Ω
    e0 = 0.5(ωx+ωy+ωz)
    ecut < e0 && error("ecut must exceed the zero point energy.")
    Mx,nx,ex = nenergy(ecut,e0,ωx,basis)
    My,ny,ey = nenergy(ecut,e0,ωy,basis)
    Mz,nz,ez = nenergy(ecut,e0,ωz,basis)
    en = [ex[i+1]+ey[j+1]+ez[k+1] for i in nx, j in ny, k in nz]
    P  = en .< ecut
    en = P.*en
    ax = sqrt(1/ωx); ay = sqrt(1/ωy); az = sqrt(1/ωz)
    X,Px = ladderops(Mx,ax)
    Y,Py = ladderops(My,ay)
    Z,Pz = ladderops(Mz,az)
    return P,en,X,Px,Y,Py,Z,Pz
  end

end

end
