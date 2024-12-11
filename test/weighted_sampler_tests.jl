
@testset "DynamicSampling.jl Tests" begin
    b = 100
    range = 1:b
    weights = [Float64(i) for i in range]
    weights2 = (Float64(i) for i in range)

    s1 = DynamicSampler()
    for i in range
    	push!(s1, i, weights[i])
    end

    @test sort!(allinds(s1)) == collect(range)
    @test all(x -> 1 <= x <= b, [rand(s1) for _ in 1:10^3])

    delete!(s1, 1)
    delete!(s1, 2)

    @test all(x -> 3 <= x <= b, [rand(s1) for _ in 1:10^3])

    s2 = DynamicSampler()
    append!(s2, range, weights2)

    @test sort!(allinds(s2)) == collect(range)
    @test all(x -> 1 <= x <= b, [rand(s2) for _ in 1:10^3])

    e1 = rand(s2)
    e2 = rand(s2)
    while e2 == e1
        e2 = rand(s2)
    end

    delete!(s2, e1)
    delete!(s2, e2)

    @test all(x -> x != e1 && x != e2 && 1 <= x <= b, [rand(s2) for _ in 1:10^3])

    rng = StableRNG(42)

    s3 = DynamicSampler(rng)
    append!(s3, range, weights2)

    samples_counts = countmap([rand(s3) for _ in 1:10^5])
    counts_est = [samples_counts[i] for i in 1:b]
    ps_exact = [i/((b ÷ 2)*(b+1)) for i in 1:b]

    chisq_test = ChisqTest(counts_est, ps_exact)
    @test pvalue(chisq_test) > 0.05

    samples_counts = countmap(rand(s3, 10^5))
    counts_est = [samples_counts[i] for i in 1:b]
    ps_exact = [i/((b ÷ 2)*(b+1)) for i in 1:b]

    chisq_test = ChisqTest(counts_est, ps_exact)
    @test pvalue(chisq_test) > 0.05

    for i in 1:(b ÷ 2)
    	delete!(s3, i)
    end

    samples_counts = countmap([rand(s3) for _ in 1:10^5])
    counts_est = [samples_counts[i] for i in (b ÷ 2 + 1):b]
    ps_exact = [i/((b ÷ 2)*(b+1) - (b ÷ 4)*(b ÷ 2 + 1)) for i in 51:b]

    chisq_test = ChisqTest(counts_est, ps_exact)
    @test pvalue(chisq_test) > 0.05

    s4 = DynamicSampler(rng)
    
    for (i, w) in zip(range, weights2)
        push!(s4, i, w)
    end

    delete!(s4, 1)
    delete!(s4, 2)

    push!(s4, 2, 200.0)
    push!(s4, 1000, 1000.0)

    samples_counts = countmap([rand(s4) for _ in 1:10^5])
    counts_est = [samples_counts[i] for i in [2:b..., 1000]]
    wsum = (b ÷ 2)*(b+1) - 3 + 200 + 1000
    ps_exact = [i == 2 ? 200/wsum : i/wsum for i in [2:b..., 1000]]

    chisq_test = ChisqTest(counts_est, ps_exact)
    @test pvalue(chisq_test) > 0.05

    @test isempty(s4) == false
    empty!(s4)
    @test isempty(s4) == true

    push!(s4, 1, 1.0)
    push!(s4, 2, 2.0)
    samples_counts = countmap([rand(s4) for _ in 1:10^4])
    counts_est = [samples_counts[1], samples_counts[2]]
    ps_exact = [1/3, 2/3]
    @test pvalue(chisq_test) > 0.05

    delete!(s4, 2)
    @test unique([rand(s4) for _ in 1:10^3]) == [1]
end
