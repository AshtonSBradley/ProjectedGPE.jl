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
  Espec = ex.*P
  Nmax = length(ex)
  @pack cinfo = basis, Ω, ecut, e0, Nmax, N, Espec, P

  return cinfo,Espec,P
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
