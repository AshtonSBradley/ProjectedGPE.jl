function x2c(ψ::Array{Complex{Float64},1},tinfo::ProjectedGPE.Tinfo)
    return tinfo.Tx'*(tinfo.W.*ψ)
end

function x2c(ψ::Array{Complex{Float64},2},tinfo::ProjectedGPE.Tinfo)
    return tinfo.Tx'*(tinfo.W.*ψ)*tinfo.Ty
end

function x2c(ψ::Array{Complex{Float64},3},tinfo::ProjectedGPE.Tinfo)
    sa = size(ψ);sx = size(tinfo.Tx');sy = size(tinfo.Ty');sz = size(tinfo.Tz')
    Ax = reshape(tinfo.W.*ψ,(sa[1],sa[2]*sa[3]))
    Ax = reshape(tinfo.Tx'*Ax,(sx[1],sa[2],sa[3]))
    Ax = permutedims(Ax,(2,3,1));sa = size(Ax)

    Ay = reshape(Ax,(sa[1],sa[2]*sa[3]))
    Ay = reshape(tinfo.Ty'*Ay,(sy[1],sa[2],sa[3]))
    Ay = permutedims(Ay,(2,3,1));sa = size(Ay)

    Az = reshape(Ay,(sa[1],sa[2]*sa[3]))
    Az = reshape(tinfo.Tz'*Az,(sz[1],sa[2],sa[3]))
    return permutedims(Az,(2,3,1))

end
