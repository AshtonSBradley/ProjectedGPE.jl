function nenergy(ecut,e0,ω=1.0,basis="Hermite")
  if basis=="Hermite"
    M = floor((ecut - e0)/ω)+1; M=Int(M)
    n = collect(0:(M-1))
    en = ω*(n+0.5)
  return M,n,en
  else
    error("basis not implemented")
  end
end
