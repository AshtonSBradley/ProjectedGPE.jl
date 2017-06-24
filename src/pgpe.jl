#PGPE nonlinearity
#1D
function nlin!(c::Array{Complex{Float64},1},dc::Array{Complex{Float64},1})
    #ψ .= Tx*c
    #dc.= Tx'*(wx.*abs.(ψ).^2.*ψ)
    anisotrans!(c,Tx,ψ)
    anisotrans!(W.*abs.(ψ).^2.*ψ,Tx',dc)
end

#2D
function nlin!(c::Array{Complex{Float64},2},dc::Array{Complex{Float64},2})
    anisotrans!(c,Tx,Ty,ψ)
    anisotrans!(W.*abs.(ψ).^2.*ψ,Tx',Ty',dc)
end

#3D
function nlin!(c::Array{Complex{Float64},3},dc::Array{Complex{Float64},3})
    anisotrans!(c,Tx,Ty,Tz,ψ)
    anisotrans!(W.*abs.(ψ).^2.*ψ,Tx',Ty',Tz',dc)
end

#dPGPE in reservoir "frame"
function Lgp!(t,c,dc)
    nlin!(c,dc)
    dc .= P.*(-im*(1-im*Γ̄)*((en - μ).*c .+ g*dc))
end
