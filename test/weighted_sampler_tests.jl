
@testset "DynamicSampling.jl Tests" begin
    b = 100
    range = 1:b
    weights = [Float64(i) for i in range]
    weights2 = (Float64(i) for i in range)

    s1 = DynamicSampler()
    for i in range
    	push!(s1, (i, weights[i]))
    end

    @test sort!(allvalues(s1)) == collect(range)
    @test all(x -> 1 <= x <= b, [rand(s1) for _ in 1:10^3])

    deleteat!(s1, 1)
    deleteat!(s1, 2)

    @test all(x -> 3 <= x <= b, [rand(s1) for _ in 1:10^3])

    s2 = DynamicSampler()
    append!(s2, (range, weights2))

    @test sort!(allvalues(s2)) == collect(range)
    @test all(x -> 1 <= x <= b, [rand(s2) for _ in 1:10^3])

    e1 = rand(s2; info=true)
    e2 = rand(s2; info=true)

    deleteat!(s2, e1)
    deleteat!(s2, e2)
    @test all(x -> x != e1.idx && x != e2.idx && 1 <= x <= b, [rand(s2) for _ in 1:10^3])

    rng = StableRNG(41)

    s3 = DynamicSampler(rng)
    append!(s3, (range, weights2))

    samples_counts = countmap([rand(s3) for _ in 1:10^5])
    counts_est = [samples_counts[i] for i in 1:b]
    ps_exact = [i/((b ÷ 2)*(b+1)) for i in 1:b]

    chisq_test = ChisqTest(counts_est, ps_exact)
    @test pvalue(chisq_test) > 0.05

    for i in 1:(b ÷ 2)
    	deleteat!(s3, i)
    end

    samples_counts = countmap([rand(s3) for _ in 1:10^5])
    counts_est = [samples_counts[i] for i in (b ÷ 2 + 1):b]
    ps_exact = [i/((b ÷ 2)*(b+1) - (b ÷ 4)*(b ÷ 2 + 1)) for i in 51:b]

    chisq_test = ChisqTest(counts_est, ps_exact)
    @test pvalue(chisq_test) > 0.05

    s4 = DynamicSampler(rng)
    
    for (i, w) in zip(range, weights2)
        push!(s4, (i, w))
    end

    deleteat!(s4, 1)
    deleteat!(s4, 2)

    push!(s4, (2, 200.0))
    push!(s4, (1000, 1000.0))

    samples_counts = countmap([rand(s4) for _ in 1:10^5])
    counts_est = [samples_counts[i] for i in [2:b..., 1000]]
    wsum = (b ÷ 2)*(b+1) - 1 + 1000
    ps_exact = [i == 2 ? 200/wsum : i/wsum for i in [2:b..., 1000]]

    chisq_test = ChisqTest(counts_est, ps_exact)
    @test pvalue(chisq_test) > 0.05
end
