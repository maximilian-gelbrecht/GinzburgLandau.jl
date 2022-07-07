using GinzburgLandau, OrdinaryDiffEq, StatsBase
using Test

@testset "GinzburgLandau.jl" begin
    include("finite_differences.jl")
    #include("pseudo_spectral.jl")
end
