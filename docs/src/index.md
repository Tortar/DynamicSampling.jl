# API

`DynamicSampler` supports these basic operations:

- `Base.rand(sp::DynamicSampler)`: samples an index according to weights.
- `Base.push!(sp::DynamicSampler, idx, weight)`: adds a new index to the sampler.
- `Base.append!(sp::DynamicSampler, indices, weights)`: adds a collection of indices to the sampler.
- `Base.delete!(sp::DynamicSampler, idx)`: removes an index from the sampler.
- `Base.setindex!(sp::DynamicSampler, idx, weight)`: updates the weight of an index in the sampler.
- `Base.empty!(sp::DynamicSampler)`: removes all indices from the sampler.
- `Base.isempty(sp::DynamicSampler)`: check if the sampler contains no indices.
- `DynamicSampling.allinds(sp::DynamicSampler)`: returns all indices in the sampler.
