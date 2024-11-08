
using DynamicSampling

using HypothesisTests
using Random
using StableRNGs
using StatsBase
using Test

@testset "DynamicSampling.jl Tests" begin
    include("package_sanity_tests.jl")
    include("weighted_sampler_tests.jl")
    include("benchmark_tests.jl")
end
