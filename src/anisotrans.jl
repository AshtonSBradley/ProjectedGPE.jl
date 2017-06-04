function anisotrans(A,Tx)
return Tx*A
end

function anisotrans(A,Tx,Ty)
return Tx*A*Ty.'
end

function anisotrans(A,Tx,Ty,Tz)
sa = size(A);sx = size(Tx);sy = size(Ty);sz = size(Tz)

  #B = complex(zeros(sx[1],sy[1],sz[1]))
  A = reshape(Tx*A[:,:],(sx[1],sa[2],sa[3]))
  A = permutedims(A,(2,3,1))
  sa = size(A);
  A = reshape(Ty*A[:,:],(sy[1],sa[2],sa[3]))
  A = permutedims(A,(2,3,1))
  sa = size(A);
  A = reshape(Tz*A[:,:],(sz[1],sa[2],sa[3]))
  A = permutedims(A,(2,3,1))
  
#=
  @inbounds for j=1:sz[1]
    @inbounds for i = 1:sz[2]
      B[:,:,j] = Tz[j,i]*anisotrans(A[:,:,i],Tx,Ty)
    end
  end
=#

return A
end
