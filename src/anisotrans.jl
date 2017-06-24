function anisotrans!(c::Array{Complex{Float64},1},Tx,ψ::Array{Complex{Float64},1})
    ψ .= Tx*c
end

function anisotrans!(c::Array{Complex{Float64},2},Tx,Ty,ψ::Array{Complex{Float64},2})
    ψ .= Tx*c*Ty'
end

function anisotrans!(c::Array{Complex{Float64},3},Tx,Ty,Tz,ψ::Array{Complex{Float64},3})
    sa = size(c);sx = size(Tx);sy = size(Ty);sz = size(Tz)
    Ax = reshape(c,(sa[1],sa[2]*sa[3]))
    Ax = reshape(Tx*Ax,(sx[1],sa[2],sa[3]))
    Ax = permutedims(Ax,(2,3,1));sa = size(Ax)

    Ay = reshape(Ax,(sa[1],sa[2]*sa[3]))
    Ay = reshape(Ty*Ay,(sy[1],sa[2],sa[3]))
    Ay = permutedims(Ay,(2,3,1));sa = size(Ay)

    Az = reshape(Ay,(sa[1],sa[2]*sa[3]))
    Az = reshape(Tz*Az,(sz[1],sa[2],sa[3]))
    ψ .= permutedims(Az,(2,3,1))
end
