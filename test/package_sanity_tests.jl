
using Aqua

@testset "Code quality" begin
    Aqua.test_all(DynamicSampling, ambiguities = false)
end
