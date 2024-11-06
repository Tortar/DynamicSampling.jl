
# DynamicSampling.jl

This package implements samplers which can be used to sample with replacement from collections
of items while being able to remove and add elements from the sampler in constant time.


# Example

```julia
julia> using DynamicSampling

julia> sampler = DynamicSampler(100);

julia> for i in 1:100
           push!(sampler, (i, Float64(i)))
       end

julia> rand(sampler)
46

julia> deleteat!(sampler, 10);
```

