
using BenchmarkTools

function b(rng, n)
    s = DynamicSampler(rng)
    sizehint!(s, n)
    append!(s, (1:n, (rand(rng) for _ in 1:n)))
    t = 0
    for i in 1:n-1
        e = rand(s; info=true)
        t += e.idx
        deleteat!(s, e)
    end
    return t
end

@testset "Benchmark Tests" begin
	rng = Xoshiro(42)
	b_small = @benchmark b($rng, 100)
	b_large = @benchmark b($rng, 10^6)

	println("Benchmark results\n")
	println("Benchmark 10^2 items: $(mean(b_small.times)/10^6) ms")
	println("Benchmark 10^6 items: $(mean(b_large.times)/10^6) ms")
end
