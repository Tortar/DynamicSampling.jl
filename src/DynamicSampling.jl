module DynamicSampling

using Random

# Strongly based on https://www.aarondefazio.com/tangentially/?p=58
# and https://github.com/adefazio/sampler
struct DynamicSampler
    N::Int
    toplevel::Int
    totweight::Base.RefValue{Float64}
    weights::Vector{Float64}
    level_weights::Vector{Float64}
    level_buckets::Vector{Vector{Int}}
    level_max::Vector{Float64}
end

function DynamicSampler(N::Int, minweight=1.0, maxweight=100.0)
	return DynamicSampler(Random.default_rng(), N, minweight, maxweight)
end
function DynamicSampler(rng, N::Int, minweight=1.0, maxweight=100.0)
    toplevel = ceil(Int, log2(maxweight))
    bottomlevel = ceil(Int, log2(minweight))
    nlevels = toplevel - bottomlevel + 1
    totweight = 0.0
    weights = zeros(Float64, N)
    level_weights = zeros(Float64, nlevels)
    level_buckets = [Int[] for _ in 1:nlevels]
    level_max = [2.0^(toplevel-i+1) for i in 1:nlevels]
    return DynamicSampler(N, toplevel, Ref(totweight), weights,
        level_weights, level_buckets, level_max)
end

struct DynamicIndex
    idx::Int
    weight::Float64
    level::Int
    idx_in_level::Int
end
DynamicIndex(idx, weight) = DynamicIndex(idx, weight, 0, 0)

function Base.push!(S::DynamicSampler, e::Tuple)
    idx, weight = e
    S.totweight[] += weight
    level = getlevel(S.toplevel, weight)
    S.level_weights[level] += weight
    push!(S.level_buckets[level], idx)
	S.weights[idx] = weight
    return S
end

function Base.append!(S::DynamicSampler, e::Tuple)
    inds, weights = e
    nlevels = zeros(Int, length(S.level_buckets))
    sumweights = 0.0
    levs = zeros(Int8, length(inds))
    for (i, w) in enumerate(weights)
        level = getlevel(S.toplevel, w)
        nlevels[level] += 1
        S.level_weights[level] += w
        S.weights[i] = w
        sumweights += w
        levs[i] = level
    end
    S.totweight[] += sumweights
    for (i, bucket) in enumerate(S.level_buckets)
        resize!(bucket, length(bucket) + nlevels[i])
        nlevels[i] = length(bucket)
    end
    for (i, id) in enumerate(inds)
        level = levs[i]
        bucket = S.level_buckets[level]
        bucket[nlevels[level]] = i
        nlevels[level] -= 1
    end
    return S
end

sample(S::DynamicSampler) = _sample(S).idx
function _sample(S::DynamicSampler)
    local level, idx, idx_in_level, weight  
    # Sample a level using the CDF method
    u = rand() * S.totweight[]
    cumulative_weight = 0.0
    for i in eachindex(S.level_weights)
        cumulative_weight += S.level_weights[i]
        level = i
        u < cumulative_weight && break
    end         
    # Now sample within the level using rejection sampling
    bucket = S.level_buckets[level]
    level_size = length(bucket)
    level_max = S.level_max[level]
    while true
        idx_in_level = rand(1:level_size)
        idx = bucket[idx_in_level]
        weight = S.weights[idx]
        rand() * level_max <= weight && break
    end     
    return DynamicIndex(idx, weight, level, idx_in_level)
end
    
function Base.deleteat!(S::DynamicSampler, idx)
    weight = S.weights[idx]
    level = getlevel(S.toplevel, weight)
    idx_in_level = findfirst(x -> x == idx, S.level_buckets[level])
    isnothing(idx_in_level) && return error()
    _deleteat!(S, idx, weight, level, idx_in_level)
    return S
end
function Base.deleteat!(S::DynamicSampler, e::DynamicIndex)
    idx, weight, level, idx_in_level = e.idx, e.weight, e.level, e.idx_in_level
    _deleteat!(S, idx, weight, level, idx_in_level)
    return S
end
function _deleteat!(S, idx, weight, level, idx_in_level)
    S.weights[idx] = 0.0
    S.totweight[] -= weight
    S.level_weights[level] -= weight
    # Swap with last element for efficent delete
    bucket = S.level_buckets[level]
    bucket[idx_in_level], bucket[end] = bucket[end], bucket[idx_in_level]
    pop!(bucket)
end

getlevel(toplevel, weight) = toplevel - ceil(Int, log2(weight)) + 1

end
