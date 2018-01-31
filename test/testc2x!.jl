#Initialize CField (dimensionless)
  basis = "Hermite"
  ecut = 30*ħ*ωy/E0
  Ω = [ωx; ωy]*t0
  cinfo = makecinfo(ecut,Ω,basis)
  tinfo = maketinfo(cinfo,4)
  @unpack en,P,M = cinfo ;Mx,My = M
  #x,wx,Tx,y,wy,Ty = makealltrans(M,4,Ω)
  #W = wx.*wy'
#test transform
  c0   = randn(Mx,My)+im*randn(Mx,My); c0=P.*c0
  #ψ0   = Tx*c0*Ty' #initial condition
  ψ0  = tinfo.Tx*c0*tinfo.Ty'
  #ψ    = Tx*c0*Ty' #a field to write to in place
  ψ = tinfo.Tx*c0*tinfo.Ty'
