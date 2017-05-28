using ProjectedGPE
using Base.Test

# write your own tests here
include("testvortexposition.jl")

@time @test testvortexposition(30)
