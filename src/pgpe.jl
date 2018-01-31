#PGPE nonlinearity
#out of place
function nlin(c)
    ψ = c2x(c,tinfo)
    return x2c(abs2.(ψ).*ψ,tinfo)
end

#in place
function nlin!(dc,c,tinfo)
    c2x!(ψ,c,tinfo)
    x2c!(dc,abs2.(ψ).*ψ,tinfo)
end

#dPGPE in reservoir frame
#out of place
function Lgp(c,p,t)
     return P.*(-im*(1-im*γ)*((en - μ).*c .+ g*nlin(c)))
end

#in place
function Lgp!(dc,c,p,t)
    nlin!(dc,c,tinfo)
    dc .= P.*(-im*(1-im*γ)*((en - μ).*c .+ g*dc))
end
