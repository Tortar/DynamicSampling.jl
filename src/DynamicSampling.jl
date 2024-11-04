module DynamicSampling

struct DynamicRangeSampler
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

struct WeightedIdx
    idx::Int
    weight::Float64
    level::Int
    idx_in_level::Int
end

function WeightedIdx(idx, weight)
    return WeightedIdx(idx, weight, 0, 0)
end

function DynamicRangeSampler(N, max_value=100.0, min_value=1.0)
    top_level = ceil(Int, log2(max_value))
    bottom_level = ceil(Int, log2(min_value))
    nlevels = top_level - bottom_level + 1
    total_weight = 0.0
    weights = zeros(Float64, N)
    level_weights = zeros(Float64, nlevels)
    level_buckets = [Int[] for _ in 1:nlevels]
    level_max = [2^Float64(top_level-i+1) for i in 1:nlevels]
    return DynamicRangeSampler(N, max_value, min_value, top_level,
        bottom_level, nlevels, Ref(total_weight), weights, level_weights,
    	level_buckets, level_max)
end

function Base.push!(S::DynamicRangeSampler, e::WeightedIdx)
    S.total_weight[] += e.weight
        
    S.weights[e.idx] = e.weight
        
    raw_level = ceil(Int, log2(e.weight))
    level = S.top_level - raw_level + 1
        
    S.level_weights[level] += e.weight
    push!(S.level_buckets[level], e.idx)
    return S
end

function sample(S::DynamicRangeSampler)
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
    return WeightedIdx(idx, idx_weight, level, idx_in_level)
end
        
function Base.pop!(S::DynamicRangeSampler, e::WeightedIdx)
    # Remove it
    S.weights[e.idx] = 0.0
    S.total_weight[] -= weight
    S.level_weights[e.level] -= weight
    # Swap with last element for efficent delete
    bucket = S.level_buckets[e.level]
    bucket[e.idx_in_level], bucket[end] = bucket[end], bucket[e.idx_in_level]
    pop!(bucket)
    return (e.idx, e.weight)
end

end
