function anisotrans(A,Tx)
return Tx*A
end

function anisotrans(A,Tx,Ty)
return Tx*A*Ty.'
end

function anisotrans(A,Tx,Ty,Tz)
  sa = size(A);  sx = size(Tx);  sy = size(Ty);  sz = size(Tz)

  B = zeros(sx[1],sy[1],sz[1])

  for j=1:sz[1]
    for i = 1:sz[2]
      B[:,:,j] = Tz[j,i]*anisotrans(A[:,:,i],Tx,Ty)
    end
  end

return B
end
