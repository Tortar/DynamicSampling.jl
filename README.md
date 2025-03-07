# DynamicSampling.jl

[![CI](https://github.com/Tortar/DynamicSampling.jl/workflows/CI/badge.svg)](https://github.com/Tortar/DynamicSampling.jl/actions?query=workflow%3ACI)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://tortar.github.io/DynamicSampling.jl/stable/)
[![codecov](https://codecov.io/gh/Tortar/DynamicSampling.jl/graph/badge.svg?token=F8W0MC53Z0)](https://codecov.io/gh/Tortar/DynamicSampling.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This package implements efficient samplers which can be used to sample from
a dynamic discrete distribution, represented by a set of pairs of indices and
weights, supporting removal, addition and sampling of elements in constant time.

# Example

```julia
julia> using DynamicSampling, Random

julia> rng = Xoshiro(42);

julia> sampler = DynamicSampler(rng);

julia> # the sampler contains indices
       for i in 1:10
           push!(sampler, i, Float64(i))
       end

julia> rand(sampler)
7

julia> delete!(sampler, 8)
DynamicSampler(indices = [1, 2, 3, 4, 5, 6, 7, 9, 10], weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 10.0])
```

# Benchmark

We can try to compare this method with the equivalent static methods from StatsBase.jl
to understand how much we pay for dynamicity.

Let's first reimplement the weighted methods of StatsBase.jl using DynamicSampling.jl

```julia
julia> using DynamicSampling

julia> function static_sample_with_replacement(rng, inds, weights, n)
           sp = DynamicSampler(rng)
           append!(sp, inds, weights)
           return rand(sp, n)
       end

julia> function static_sample_without_replacement(rng, inds, weights, n)
           sp = DynamicSampler(rng)
           append!(sp, inds, weights)
           s = Vector{Int}(undef, n)
           for i in eachindex(s)
               idx = rand(sp)
               delete!(sp, idx)
               s[i] = idx
           end
           return s
       end
```

let's look at some benchmarks in respect to `StatsBase.jl` which depict the
worst case scenario where the dynamic methods are used statically. We
will use a small and a big `n` in respect to the number of indices

```julia
julia> using StatsBase, Random, BenchmarkTools

julia> rng = Xoshiro(42);

julia> inds = 1:10^6

julia> weights = Float64.(inds)

julia> n_small = 10^2

julia> n_big = 5*10^5

julia> t1_d = @benchmark static_sample_with_replacement($rng, $inds, $weights, $n_small);

julia> t1_s = @benchmark sample($rng, $inds, $(Weights(weights)), $n_small; replace=true);

julia> t2_d = @benchmark static_sample_with_replacement($rng, $inds, $weights, $n_big);

julia> t2_s = @benchmark sample($rng, $inds, $(Weights(weights)), $n_big; replace=true);

julia> t3_d = @benchmark static_sample_without_replacement($rng, $inds, $weights, $n_small);

julia> t3_s = @benchmark sample($rng, $inds, $(Weights(weights)), $n_small; replace=false);

julia> t4_d = @benchmark static_sample_without_replacement($rng, $inds, $weights, $n_big);

julia> t4_s = @benchmark sample($rng, $inds, $(Weights(weights)), $n_big; replace=false);

julia> using StatsPlots

julia> times_static = mean.([t1_s.times, t2_s.times, t3_s.times, t4_s.times]) ./ 10^6

julia> times_dynamic = mean.([t1_d.times, t2_d.times, t3_d.times, t4_d.times]) ./ 10^6

julia> groupedbar(
           ["small wr", "big wr", "small wor", "big wor"], [times_static times_dynamic], 
           ylabel="time (ms)", labels=["static" "dynamic"], dpi=1200
       )
```

<img src="https://github.com/user-attachments/assets/c654359a-3729-4855-9a94-085ff2d88521" width="500" />

From the figure, we can conclude that the dynamic versions are quite competitive even
in this worst case analysis.

Another insightful benchmark is the time for a single random draw in respect to the number of indices in the sampler:

```julia
using Random, BenchmarkTools, Plots

using DynamicSampling

q = 1
rng = Xoshiro(42)
ts1 = []
ts2 = []
for n in [10^i for i in 1:8]
    sp = DynamicSampler(rng)
    append!(sp, 1:n, 1:1/q:10)
    b1 = @benchmark rand($sp)
    b2 = @benchmark rand($sp, 10^4)
    push!(ts1, mean(b1.times))
    push!(ts2, mean(b2.times)./10^4)
    q += n
end

plot!(1:8, ts1, markershape=:circle, xlabel="number of indices in the sampler", 
     ylabel="time for one draw (ns)", xticks = (1:8, ["\$10^$(i)\$" for i in 1:8]),
     guidefontsize=8, dpi = 1200, label="single", ytick=[10*x for x in 1:10]
)
plot!(1:8, ts2, markershape=:circle, xlabel="number of indices in the sampler", 
     ylabel="time for one draw (ns)", xticks = (1:8, ["\$10^$(i)\$" for i in 1:8]),
     guidefontsize=8, dpi = 1200, label="bulk", ytick=[10*x for x in 1:10]
)
```

<img src="https://github.com/user-attachments/assets/0261b65b-cb58-49c9-bc0a-30f094ae90bb" width="500" />

This hints on the fact that the operation becomes essentially memory bound when the number of indices surpass roughly 1 million elements,
in particular in the case of single-instance random draws.

# References

References for the techniques implemented in this library can be found in
[A constant-time kinetic Monte Carlo algorithm for simulation of large biochemical reaction networks](https://www.researchgate.net/publication/5338017_A_constant-time_kinetic_Monte_Carlo_algorithm_for_simulation_of_large_biochemical_reaction_networks) and [Weighted random sampling with replacement with dynamic weights](https://www.aarondefazio.com/tangentially/?p=58)

