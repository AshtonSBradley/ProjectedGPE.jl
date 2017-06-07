function makecfieldinfo(basis,ecut,Ω)

if basis=="hermite"
  dim=length(Ω)
  if dim==1
  ω1 = Ω
  e0 = 0.5*ω1
  ecut < e0 && error("ecut is smaller than the zero point energy!")
  N1 = floor((ecut - e0)/ω1)+1
  ex = ω1*(collect(0:(N1-1))+0.5);ex=ex[:]
  Px = ex .< ecut
  Espec= ex.*Px;
  Espec = Espec'
  Px = Px.'
  N = length(ex)
  cinfo = cfieldinfo(basis=basis,ecut=ecut,e0=e0,Espec=Espec[:],N=N,N1=N1,ω1=ω1,Px=Px[:])
  return cinfo
  #return basis, ecut, e0, Espec, N, N1, ω1, Px
  elseif dim==2
  throw(DomainError())
  elseif dim==3
  throw(DomainError())
  else error("dim must be <=3.")
  end
end
#out = struct('basis',basis,'NN',NN,'Ecut',Ecut,'E0',E0,...
#    'f1',f1,'f2',f2,'f3',f3,'N1',N1,'N2',N2,'N3',N3,'dim',dim,'eta',eta,'Stype',Stype);
end
#=
function evalues(ecut,ω1,ω2,basis="hermite")

end

function evalues(ecut,ω1,ω2,ω3,basis="hermite")

end
=#
