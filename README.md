# DynamicSampling.jl

[![CI](https://github.com/Tortar/DynamicSampling.jl/workflows/CI/badge.svg)](https://github.com/Tortar/DynamicSampling.jl/actions?query=workflow%3ACI)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://tortar.github.io/DynamicSampling.jl/stable/)
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

If you happen to require to remove already sampled elements, 
you can use instead

```julia
julia> i = rand(sampler; info=true)
IndexInfo(idx = 9, weight = 9.0)

julia> delete!(sampler, i)
DynamicSampler(indices = [1, 2, 3, 4, 5, 6, 7, 10], weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0])
```

which is a bit more efficient than removing not sampled indices.

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
               e = rand(sp; info=true)
               s[i] = e.idx
               delete!(sp, e)
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
           ["small sample wr", "big sample wr", "small sample wor", "big sample wor"],
           [times_static times_dynamic], ylabel="time (ms)", labels=["static" "dynamic"]
       )
```

<img src="https://github.com/user-attachments/assets/77c23c51-9452-4e31-ad1c-33bb2aa9c622" width="500" />

From the figure, we can conclude that the dynamic versions are quite competitive even
in this worst case analysis.




