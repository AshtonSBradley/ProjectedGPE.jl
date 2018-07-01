function x2c!(c::Array{Complex{Float64},1},ψ::Array{Complex{Float64},1},tinfo::ProjectedGPE.Tinfo)
    c = tinfo.Tx'*(ψ.*tinfo.W)
end

function x2c!(c::Array{Complex{Float64},2},ψ::Array{Complex{Float64},2},tinfo::ProjectedGPE.Tinfo)
    c = tinfo.Tx'*(ψ.*tinfo.W)*tinfo.Ty
end

function x2c!(c::Array{Complex{Float64},3},ψ::Array{Complex{Float64},3},tinfo::ProjectedGPE.Tinfo)
    sa = size(ψ);sx = size(tinfo.Tx');sy = size(tinfo.Ty');sz = size(tinfo.Tz')
    Ax = reshape(ψ.*tinfo.W,(sa[1],sa[2]*sa[3]))
    Ax = reshape(tinfo.Tx'*Ax,(sx[1],sa[2],sa[3]))
    Ax = permutedims(Ax,(2,3,1));sa = size(Ax)

    Ay = reshape(Ax,(sa[1],sa[2]*sa[3]))
    Ay = reshape(tinfo.Ty'*Ay,(sy[1],sa[2],sa[3]))
    Ay = permutedims(Ay,(2,3,1));sa = size(Ay)

    Az = reshape(Ay,(sa[1],sa[2]*sa[3]))
    Az = reshape(tinfo.Tz'*Az,(sz[1],sa[2],sa[3]))
    ψ = permutedims(Az,(2,3,1))

end
