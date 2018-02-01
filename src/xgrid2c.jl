    function xgrid2c(ψ::Array{Complex{Float64},1},cinfo::ProjectedGPE.Cinfo)
        error("implement")
    end

    function xgrid2c(ψ::Array{Complex{Float64},2},x,y,cinfo::ProjectedGPE.Cinfo)
        Mx,My = cinfo.M
        dx = x[2]-x[1];dy = y[2]-y[1]
        Tx = eigmat(Mx,x)
        Ty = eigmat(My,y)
        return Tx'*(ψ*Ty)*dx*dy
    end

    function x2c(ψ::Array{Complex{Float64},3},cinfo::ProjectedGPE.Cinfo)
        error("implement")
    end
