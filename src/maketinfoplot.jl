function maketinfoplot(cinfo,x,y)
tinfo = Tinfo()

@unpack en,P,M,Ω,basis = cinfo
ωx,ωy=Ω
Tx = eigmat(Mx,x/x0,ω=ωx/ωy)
Ty = eigmat(My,y/x0,ω=ωy/ωy)

@pack tinfo = Tx,Ty

return tinfo
end
