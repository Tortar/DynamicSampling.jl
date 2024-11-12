# DynamicSampling.jl

[![CI](https://github.com/Tortar/DynamicSampling.jl/workflows/CI/badge.svg)](https://github.com/Tortar/DynamicSampling.jl/actions?query=workflow%3ACI)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This package implements efficient samplers which can be used to sample from
a set of indices while being able to remove, add and sample elements from
the sampler in constant time.

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

julia> deleteat!(sampler, 8)
DynamicSampler(indices = [1, 2, 3, 4, 5, 6, 7, 9, 10], weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 10.0])
```

If you happen to require to remove already sampled elements, 
you can use instead

```julia
julia> i = rand(sampler; info=true)
IndexInfo(idx = 9, weight = 9.0)

julia> deleteat!(sampler, i)
DynamicSampler(indices = [1, 2, 3, 4, 5, 6, 7, 10], weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0])
```

which is a bit more efficient than removing not sampled indices.
