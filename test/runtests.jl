using ProjectedGPE
using Base.Test

# Run tests
tic()
include("modeorthonorm.jl")
@testset "PGPE tests" begin include("pgpetests.jl") end
toc()
