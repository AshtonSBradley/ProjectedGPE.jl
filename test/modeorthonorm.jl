function modeorthonorm(basis,N,err = 1e-9,ω=1.)

x,w,T = nfieldtrans(N,2,ω=ω)
Normtest = T'*(T.*(w.*ones(N)))

#check orthonormality is close to identity matrix for the basis
maxerr = maximum(abs.(Normtest-eye(Normtest)))
maxerr < err
end
