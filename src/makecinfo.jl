function makecinfo(basis,ecut,Ω)
cinfo = CInfo()
if basis=="Hermite"
  dim=length(Ω)
  if dim==1
    ωx = Ω[1]
    e0 = 0.5*ωx
    ecut < e0 && error("ecut is smaller than the zero point energy!")
    Mx = floor((ecut - e0)/ωx)+1;Mx=Int(Mx)
    M=[Mx]
    nx = collect(0:(Mx-1))
    ex = ωx*(nx+0.5);en=ex
    P = en .< ecut
    en = P.*en
    Mmax = Mx

  elseif dim==2
    ωx,ωy = Ω
    e0 = 0.5(ωx+ωy)
    ecut < e0 && error("ecut is smaller than the zero point energy!")
    Mx = floor((ecut - e0)/ωx)+1;Mx=Int(Mx)
    My = floor((ecut - e0)/ωy)+1;My=Int(My)
    M  = [Mx,My]
    nx = collect(0:(Mx-1));ny = collect(0:(My-1))
    ex = (nx+0.5)ωx;ey = (ny+0.5)ωy
    en = [ex[i+1]+ey[j+1] for i=nx, j=ny]
    P = en .< ecut
    en = P.*en
    Mmax = maximum(M)
  elseif dim==3
    ωx,ωy,ωz = Ω
    e0 = 0.5(ωx+ωy+ωz)
    ecut < e0 && error("ecut is smaller than the zero point energy!")
    Mx = floor((ecut - e0)/ωx)+1;Mx=Int(Mx)
    My = floor((ecut - e0)/ωy)+1;My=Int(My)
    Mz = floor((ecut - e0)/ωz)+1;Mz=Int(Mz)
    M  = [Mx,My,Mz]
    nx = collect(0:(Mx-1));ny = collect(0:(My-1));nz = collect(0:(Mz-1))
    ex = (nx+0.5)ωx;ey = (ny+0.5)ωy;ez = (nz+0.5)ωz
    en = [ex[i+1]+ey[j+1]+ez[k+1] for i=nx, j=ny, k=nz]
    P = en .< ecut
    en = P.*en
    Mmax = maximum(M)
else error("Spatial dimension must be 1, 2, or 3.")
  end
end
  @pack cinfo = basis, Ω, ecut, e0, Mmax, M, en, P
  return cinfo
end
