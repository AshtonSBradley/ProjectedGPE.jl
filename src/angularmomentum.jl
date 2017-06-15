function angularmomentum(basis,cfield,Ω)
if basis=="Hermite"
  dim=length(Ω)
  if dim==1
    error("L not defined in 1D")
  elseif dim==2
    Lz = -im*(X.*Py - Y.*Px)
  elseif dim=3
    Lx = -im*(Ŷ.*P̂z - Ẑ.*P̂y)
    Ly = -im*(Ẑ.*P̂x - X̂.*P̂z)
    Lz = -im*(X̂.*P̂y - Ŷ.*P̂x)
  end
end
