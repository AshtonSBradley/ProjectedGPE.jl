using ProjectedGPE
using Base.Test


# Run tests
tic()
include("makepsi.jl")
include("checkpositions.jl")
include("testallpositions.jl")
include("testallcharges.jl")
include("modeorthonorm.jl")
@testset "Vortex charge and position tests" begin include("vortextests.jl") end
@testset "PGPE tests" begin include("pgpetests.jl") end
toc()
