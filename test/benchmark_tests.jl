
using BenchmarkTools

function b1(rng, n)
    s = DynamicSampler(rng)
    sizehint!(s, n)
    append!(s, 1:n, (rand(rng) for _ in 1:n))
    t = 0
    for i in 1:n-1
        e = rand(s; info=true)
        t += e.idx
        delete!(s, e)
    end
    return t
end

function b2(rng, n)
    s = DynamicSampler(rng)
    sizehint!(s, n)
    for i in 1:n
        push!(s, i, rand(rng))
    end
    t = 0
    for i in 1:n-1
        idx = rand(s)
        t += idx
        delete!(s, idx)
    end
    return t
end

@testset "Benchmark Tests" begin
    rng = Xoshiro(42)
    b_small_fast = @benchmark b1($rng, 100)
    b_large_fast = @benchmark b1($rng, 10^6)
    b_small_slow = @benchmark b2($rng, 100)
    b_large_slow = @benchmark b2($rng, 10^6)
    
    println("Benchmark results\n")
    println("Benchmark Fast 10^2 items: $(mean(b_small_fast.times)/10^6) ms")
    println("Benchmark Fast 10^6 items: $(mean(b_large_fast.times)/10^6) ms")
    println("Benchmark Slow 10^2 items: $(mean(b_small_slow.times)/10^6) ms")
    println("Benchmark Slow 10^6 items: $(mean(b_large_slow.times)/10^6) ms")
end
