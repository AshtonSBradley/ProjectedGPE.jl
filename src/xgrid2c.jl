    function xgrid2c(ψ::Array{Complex{Float64},1},cinfo::ProjectedGPE.Cinfo)
        error("implement")
    end

    function xgrid2c(ψ::Array{Complex{Float64},2},cinfo::ProjectedGPE.Cinfo)
        Mx,My = cinfo.M
        Tx = eigmat(Mx,x/x0)
        Ty = eigmat(My,y/x0)
        return Tx'*(ψ*Ty)*Δx*Δy/x0^2
    end

    function x2c(ψ::Array{Complex{Float64},3},cinfo::ProjectedGPE.Cinfo)
        error("impliment")
    end
