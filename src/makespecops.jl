function makespecops(basis,ecut,Ω)
if basis=="Hermite"
  dim=length(Ω)
  if dim==1
    ωx = Ω[1]
    e0 = 0.5*ωx
    ecut < e0 && error("ecut is smaller than the zero point energy!")
    Mx,nx,en = nenergy("Hermite",ecut,e0,ωx)
    P  = en .< ecut
    en = P.*en
    ax = sqrt(1/ωx) #in dimensionless units
    X,Px = makeladderops(Mx,ax)
    return P,en,X,Px
  elseif dim==2
    ωx,ωy = Ω
    e0 = 0.5(ωx+ωy)
    ecut < e0 && error("ecut is smaller than the zero point energy!")
    Mx,nx,ex = nenergy("Hermite",ecut,e0,ωx)
    My,ny,ey = nenergy("Hermite",ecut,e0,ωy)
    en = [ex[i+1]+ey[j+1] for i=nx, j=ny]
    P  = en .< ecut
    en = P.*en
    ax = sqrt(1/ωx);ay = sqrt(1/ωy)
    X,Px = makeladderops(Mx,ax)
    Y,Py = makeladderops(My,ay)
    return P,en,X,Px,Y,Py
  elseif dim==3
    ωx,ωy,ωz = Ω
    e0 = 0.5(ωx+ωy+ωz)
    ecut < e0 && error("ecut is smaller than the zero point energy!")
    Mx,nx,ex = nenergy("Hermite",ecut,e0,ωx)
    My,ny,ey = nenergy("Hermite",ecut,e0,ωy)
    Mz,nz,ez = nenergy("Hermite",ecut,e0,ωz)
    en = [ex[i+1]+ey[j+1]+ez[k+1] for i=nx, j=ny, k=nz]
    P  = en .< ecut
    en = P.*en
    ax = sqrt(1/ωx);ay = sqrt(1/ωy);az = sqrt(1/ωz)
    X,Px = makeladderops(Mx,ax)
    Y,Py = makeladderops(My,ay)
    Z,Pz = makeladderops(Mz,az)
    return P,en,X,Px,Y,Py,Z,Pz
  end

end

end
