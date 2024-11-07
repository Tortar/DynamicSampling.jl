
# DynamicSampling.jl

This package implements efficient samplers which can be used to sample from
(weighted) indices while being able to remove and add elements from the sampler in 
constant time.

# Example

```julia
julia> using DynamicSampling

julia> sampler = DynamicSampler();

julia> for i in 1:10
           push!(sampler, (i, Float64(i)))
       end

julia> rand(sampler)
7

julia> deleteat!(sampler, 7)
DynamicSampler(indices = [1, 2, 3, 4, 5, 6, 8, 9, 10], weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 9.0, 10.0])
```

Importantly using `deleteat!` as above will incur a non-constant overhead. However,
if you happen to require to remove already sampled elements, you can use instead

```julia
julia> i = rand(sampler; info=true)
IndexInfo(idx = 9, weight = 9.0)

julia> deleteat!(sampler, i)
DynamicSampler(indices = [1, 2, 3, 4, 5, 6, 8, 10], weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0])
```

which will be a `O(1)` operation.

