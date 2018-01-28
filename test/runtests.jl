using ProjectedGPE
using Base.Test

# Run tests
tic()
include("modeorthonorm.jl")
@testset "Hermite basis tests" begin include("testHermites.jl") end
toc()
