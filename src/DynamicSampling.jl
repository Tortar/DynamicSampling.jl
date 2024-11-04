module DynamicSampling

# Strongly based on https://www.aarondefazio.com/tangentially/?p=58
# and https://github.com/adefazio/sampler
struct DynamicSampler
    N::Int
    max_value::Float64
    min_value::Float64
    top_level::Int
    bottom_level::Int
    nlevels::Int
    total_weight::Base.RefValue{Float64}
    weights::Vector{Float64}
    level_weights::Vector{Float64}
    level_buckets::Vector{Vector{Int}}
    level_max::Vector{Float64}
end
function DynamicSampler(N, max_value=100.0, min_value=1.0)
    top_level = ceil(Int, log2(max_value))
    bottom_level = ceil(Int, log2(min_value))
    nlevels = top_level - bottom_level + 1
    total_weight = 0.0
    weights = zeros(Float64, N)
    level_weights = zeros(Float64, nlevels)
    level_buckets = [Int[] for _ in 1:nlevels]
    level_max = [2^Float64(top_level-i+1) for i in 1:nlevels]
    return DynamicSampler(N, max_value, min_value, top_level,
        bottom_level, nlevels, Ref(total_weight), weights, level_weights,
        level_buckets, level_max)
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
    return _push!(S, idx, weight)
end
function Base.push!(S::DynamicSampler, e::DynamicIndex)
    idx, weight = e.idx, e.weight
    return _push!(S, idx, weight)
end

function _push!(S::DynamicSampler, idx, weight)
    S.total_weight[] += weight
        
    S.weights[idx] = weight
        
    level = getlevel(S, weight)
        
    S.level_weights[level] += weight
    push!(S.level_buckets[level], idx)
    return S
end

function _sample(S::DynamicSampler)
    local level, idx, idx_in_level, idx_weight

    u = rand() * S.total_weight[]
        
    # Sample a level using the CDF method
    cumulative_weight = 0.0
    for i in eachindex(S.level_weights)
        cumulative_weight += S.level_weights[i]
        level = i
        u < cumulative_weight && break
    end
            
    # Now sample within the level using rejection sampling
    level_size = length(S.level_buckets[level])
    level_max = S.level_max[level]
    while true
        idx_in_level = rand(1:level_size)
        idx = S.level_buckets[level][idx_in_level]
        idx_weight = S.weights[idx]
        u_lvl = rand() * level_max
        u_lvl <= idx_weight && break
    end     
    return DynamicIndex(idx, idx_weight, level, idx_in_level)
end

sample(S::DynamicSampler) = _sample(S).idx
    
function Base.pop!(S::DynamicSampler, idx)
    weight = S.weights[idx]
    level = getlevel(S, weight)
    idx_in_level = findfirst(x -> x == idx, S.level_buckets[level])
    isnothing(idx_in_level) && return error()
    _pop!(S, idx, weight, level, idx_in_level)
    return idx
end
function Base.pop!(S::DynamicSampler, e::DynamicIdx)
    idx, weight, level, idx_in_level = e.idx, e.weight, e.level, e.idx_in_level
    _pop!(S, idx, weight, level, idx_in_level)
    return e
end

function _pop!(S, idx, weight, level, idx_in_level)
    S.weights[idx] = 0.0
    S.total_weight[] -= weight
    S.level_weights[level] -= weight
    # Swap with last element for efficent delete
    bucket = S.level_buckets[level]
    bucket[idx_in_level], bucket[end] = bucket[end], bucket[idx_in_level]
    pop!(bucket)
end

function getlevel(S::DynamicSampler, weight)
    raw_level = ceil(Int, log2(weight))
    return S.top_level - raw_level + 1
end

end
