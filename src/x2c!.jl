function x2c!(c::Array{Complex{Float64},1},ψ::Array{Complex{Float64},1},Tx)
    c2x!(c,ψ.*W,Tx')
end

function x2c!(c::Array{Complex{Float64},2},ψ::Array{Complex{Float64},2},Tx,Ty)
    c2x!(c,ψ.*W,Tx',Ty')
end

function x2c!(c::Array{Complex{Float64},3},ψ::Array{Complex{Float64},3},Tx,Ty,Tz)
    c2x!(c,ψ.*W,Tx',Ty',Tz')
end
