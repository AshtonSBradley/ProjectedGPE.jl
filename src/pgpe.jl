#PGPE nonlinearity
#1D
function nlin!(c::Array{Complex{Float64},1},cinfo::ProjectedGPE.CInfo,dc::Array{Complex{Float64},1})
    #ψ .= Tx*c
    #dc.= Tx'*(wx.*abs.(ψ).^2.*ψ)
    anisotrans!(c,cinfo.Tx,cinfo.ψ)
    anisotrans!(cinfo.W.*abs.(cinfo.ψ).^2.*cinfo.ψ,cinfo.Tx',dc)
end

#2D
function nlin!(c::Array{Complex{Float64},2},cinfo::ProjectedGPE.CInfo,dc::Array{Complex{Float64},2})
    anisotrans!(c,cinfo.Tx,cinfo.Ty,cinfo.ψ)
    anisotrans!(cinfo.W.*abs.(cinfo.ψ).^2.*cinfo.ψ,cinfo.Tx',cinfo.Ty',dc)
end

#3D
function nlin!(c::Array{Complex{Float64},3},cinfo::ProjectedGPE.CInfo,dc::Array{Complex{Float64},3})
    anisotrans!(c,cinfo.Tx,cinfo.Ty,cinfo.Tz,cinfo.ψ)
    anisotrans!(cinfo.W.*abs.(cinfo.ψ).^2.*cinfo.ψ,cinfo.Tx',cinfo.Ty',cinfo.Tz',dc)
end

#dPGPE in "reservoir frame"
function Lgp!(t,c,cinfo,dc)
    nlin!(c,cinfo,dc)
    dc .= cinfo.P.*(-im*(1-im*cinfo.Γ̄)*((cinfo.en - cinfo.μ).*c .+ cinfo.g*dc))
end
