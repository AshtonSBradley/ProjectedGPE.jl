function maketinfo(cinfo,n)
tinfo = Tinfo()

@unpack en,P,M = cinfo 
x,wx,Tx,y,wy,Ty = makealltrans(M,Ω,n=n,basis=basis)
W = wx.*wy'
@pack tinfo = Tx,Ty,W

return tinfo
end
