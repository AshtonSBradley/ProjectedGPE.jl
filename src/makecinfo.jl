function makecinfo(ecut,Ω,basis="Hermite")
cinfo = Cinfo()
if basis=="Hermite"
  dim=length(Ω)
  if dim==1
    ωx = Ω[1]
    e0 = 0.5*ωx
    ecut < e0 && error("ecut must exceed the zero point energy.")
    Mx,nx,en = nenergy(ecut,e0,ωx,basis)
    Mmax = Mx
    P  = en .< ecut; Mult = sum(P)
    en = P.*en
    M  = [Mx]
  elseif dim==2
    ωx,ωy = Ω
    e0 = 0.5(ωx+ωy)
    ecut < e0 && error("ecut must exceed the zero point energy.")
    Mx,nx,ex = nenergy(ecut,e0,ωx,basis)
    My,ny,ey = nenergy(ecut,e0,ωy,basis)
    M  = [Mx,My]; Mmax = maximum(M)
    en = [ex[i+1]+ey[j+1] for i=nx, j=ny]
    P  = en .< ecut; Mult = sum(P)
    en = P.*en
    M  = [Mx,My]
  elseif dim==3
    ωx,ωy,ωz = Ω
    e0 = 0.5(ωx+ωy+ωz)
    ecut < e0 && error("ecut must exceed the zero point energy.")
    Mx,nx,ex = nenergy(ecut,e0,ωx,basis)
    My,ny,ey = nenergy(ecut,e0,ωy,basis)
    Mz,nz,ez = nenergy(ecut,e0,ωz,basis)
    M  = [Mx,My,Mz]; Mmax = maximum(M)
    en = [ex[i+1]+ey[j+1]+ez[k+1] for i=nx, j=ny, k=nz]
    P  = en .< ecut; Mult = sum(P)
    en = P.*en
else error("Spatial dimension must be 1, 2, or 3.")
  end
end
  @pack cinfo = basis, Ω, ecut, e0, Mmax, M, en, P, Mult
  return cinfo
end
