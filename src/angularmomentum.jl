"""
    Lx,Ly,Lz = angularmomentum(c,M,a,basis="Hermite")

Construct a compact representation of angular momentum operators in the given basis using step operators.

`basis`: name of the orthonomal mode space. Currently only `"Hermite"` is implemented.

`c`: cfield array on which `L_j` are to act.

`M`: vector of cutoff quantum numbers, i.e. `Mx,My = M` in 2D.

`a`: vector of oscillator lengths, i.e. `ax,ay = a` in 2D.

"""
function angularmomentum(c,M,a,basis="Hermite")
if basis=="Hermite"
  dim=length(size(c))
  if dim==1
    error("L not defined in 1D")
  elseif dim==2
    ax,ay = a
    Mx,My = M
    X,Px = ladderops(Mx,ax,basis)
    Y,Py = ladderops(My,ay,basis)
    Lz = -im*(X*c*Py - Px*c*Y)
    return Lz
  elseif dim==3
    ax,ay,az = a
    Mx,My,Mz = M
    X,Px = ladderops(Mx,ax,basis)
    Y,Py = ladderops(My,ay,basis)
    Z,Pz = ladderops(My,ay,basis)
    Lx = -im*(Y.*Pz - Z.*P̂y)
    Ly = -im*(Ẑ.*P̂x - X̂.*P̂z)
    Lz = -im*(X̂.*P̂y - Ŷ.*P̂x)
    return Lx,Ly,Lz
  end
else error("basis not implemented")
end
end
