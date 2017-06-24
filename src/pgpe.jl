#PGPE nonlinearity
#1D
function nlin!(c::Array{Complex{Float64},1},cinfo,dc::Array{Complex{Float64},1})
    #ψ .= Tx*c
    #dc.= Tx'*(wx.*abs.(ψ).^2.*ψ)
    anisotrans!(c,cinfo.Tx,ψ)
    anisotrans!(W.*abs.(ψ).^2.*ψ,cinfo.Tx',dc)
end

#2D
function nlin!(c::Array{Complex{Float64},2},cinfo,dc::Array{Complex{Float64},2})
    anisotrans!(c,cinfo.Tx,cinfo.Ty,ψ)
    anisotrans!(W.*abs.(ψ).^2.*ψ,cinfo.Tx',cinfo.Ty',dc)
end

#3D
function nlin!(c::Array{Complex{Float64},3},cinfo,dc::Array{Complex{Float64},3})
    anisotrans!(c,cinfo.Tx,cinfo.Ty,cinfo.Tz,ψ)
    anisotrans!(W.*abs.(ψ).^2.*ψ,cinfo.Tx',cinfo.Ty',cinfo.Tz',dc)
end

#dPGPE in reservoir "frame"
function Lgp!(t,c,cinfo,dc)
    nlin!(c,cinfo,dc)
    dc .= P.*(-im*(1-im*Γ̄)*((en - μ).*c .+ g*dc))
end
