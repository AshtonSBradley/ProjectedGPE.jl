function maketinfoplot(cinfo,x,y)
tinfo = Tinfo()

@unpack en, P, M, Ω, basis = cinfo
ωx,ωy=Ω; Mx,My=M;
Tx = eigmat(Mx,x,ω=ωx/ωy)
Ty = eigmat(My,y,ω=ωy/ωy)

@pack tinfo = Tx, Ty

return tinfo
end
