function anisotrans(A,Tx)
return Tx*A
end

function anisotrans(A,Tx,Ty)
return Tx*A*Ty.'
end

function anisotrans(A,Tx,Ty,Tz)
sa = size(A);sx = size(Tx);sy = size(Ty);sz = size(Tz)
#!isdefined(Ax) ? Ax = complex(zeros(sx[1],sa[2],sa[3])) : nothing
  Ax = reshape(Tx*A[:,:],(sx[1],sa[2],sa[3]))
  Ax = permutedims(Ax,(2,3,1))
  sa = size(Ax)
  #!isdefined(Ay) ? Ay = complex(zeros(sy[1],sa[2],sa[3])) : nothing
  Ay = reshape(Ty*Ax[:,:],(sy[1],sa[2],sa[3]))
  Ay = permutedims(Ay,(2,3,1))
  sa = size(Ay)
  #!isdefined(Az) ? Az = complex(zeros(sz[1],sa[2],sa[3])) : nothing
  Az = reshape(Tz*Ay[:,:],(sz[1],sa[2],sa[3]))
  Az = permutedims(Az,(2,3,1))

#todo: preallocate sized arrays to avoid any memory allocations
return Az
end
