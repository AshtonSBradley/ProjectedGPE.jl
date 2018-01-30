function maketinfo(cinfo,n)
tinfo = Tinfo()

@unpack en,P,M,Ω,basis = cinfo
x,wx,Tx,y,wy,Ty = makealltrans(M,n,Ω,basis)
W = wx.*wy'
@pack tinfo = Tx,Ty,W

return tinfo
end
