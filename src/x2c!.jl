function x2c!(c::Array{Complex{Float64},1},ψ::Array{Complex{Float64},1},tinfo::ProjectedGPE.Tinfo)
    c2x!(c,ψ.*tinfo.W,tinfo.Tx')
end

function x2c!(c::Array{Complex{Float64},2},ψ::Array{Complex{Float64},2},tinfo::ProjectedGPE.Tinfo)
    c2x!(c,ψ.*tinfo.W,tinfo.Tx',tinfo.Ty')
end

function x2c!(c::Array{Complex{Float64},3},ψ::Array{Complex{Float64},3},tinfo::ProjectedGPE.Tinfo)
    c2x!(c,ψ.*tinfo.W,tinfo.Tx',tinfo.Ty',tinfo.Tz')
end
