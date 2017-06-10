function makecinfo(basis,ecut,Ω)
cinfo = CInfo()
if basis=="Hermite"
  dim=length(Ω)
  if dim==1
    ωx = Ω[1]
    e0 = 0.5*ωx
    ecut < e0 && error("ecut is smaller than the zero point energy!")
    Nx = floor((ecut - e0)/ωx)+1;N=[Nx]
    nx = collect(0:(Nx-1))
    ex = ωx*(nx+0.5)
    P = ex .< ecut
    Espec = P.*ex
    Nmax = length(ex)

  elseif dim==2
    ωx,ωy = Ω
    e0 = 0.5(ωx+ωy)
    ecut < e0 && error("ecut is smaller than the zero point energy!")
    Nx = floor((ecut - e0)/ωx)+1;Nx=Int(Nx)
    Ny = floor((ecut - e0)/ωy)+1;Ny=Int(Ny)
    N  = [Nx,Ny]
    nx = collect(0:(Nx-1));ny = collect(0:(Ny-1))
    ex = (nx+0.5)ωx;ey = (ny+0.5)ωy
    Espec = [ex[i+1]+ey[j+1] for i=nx, j=ny]
    P = Espec .< ecut
    Espec = P.*Espec
    Nmax = maximum(N)
  elseif dim==3
    ωx,ωy,ωz = Ω
    e0 = 0.5(ωx+ωy+ωz)
    ecut < e0 && error("ecut is smaller than the zero point energy!")
    Nx = floor((ecut - e0)/ωx)+1;Nx=Int(Nx)
    Ny = floor((ecut - e0)/ωy)+1;Ny=Int(Ny)
    Nz = floor((ecut - e0)/ωz)+1;Nz=Int(Nz)
    N  = [Nx,Ny,Nz]
    nx = collect(0:(Nx-1));ny = collect(0:(Ny-1));nz = collect(0:(Nz-1))
    ex = (nx+0.5)ωx;ey = (ny+0.5)ωy;ez = (nz+0.5)ωz
    Espec = [ex[i+1]+ey[j+1]+ez[k+1] for i=nx, j=ny, k=nz]
    P = Espec .< ecut
    Espec = P.*Espec
    Nmax = maximum(N)
else error("Spatial dimension must be 1, 2, or 3.")
  end
end
  @pack cinfo = basis, Ω, ecut, e0, Nmax, N, Espec, P
  return cinfo
end
