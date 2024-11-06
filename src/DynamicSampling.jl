# Strongly based on https://www.aarondefazio.com/tangentially/?p=58
# and https://github.com/adefazio/sampler

module DynamicSampling

export DynamicSampler

using Random

struct DynamicSampler{R}
    rng::R
    N::Int
    totvalues::Base.RefValue{Int}
    totweight::Base.RefValue{Float64}
    weights::Vector{Float64}
    level_weights::Vector{Float64}
    level_buckets::Vector{Vector{Int}}
    level_max::Vector{Float64}
    level_inds::Vector{Int}
end

DynamicSampler(N::Int) = DynamicSampler(Random.default_rng(), N)
function DynamicSampler(rng, N::Int)
    totweight = 0.0
    weights = zeros(Float64, N)
    level_weights = Float64[0.0]
    level_buckets = [Int[],]
    level_max = Float64[]
    level_inds = Int[]
    return DynamicSampler(rng, N, Ref(0), Ref(totweight), weights, level_weights, 
    	level_buckets, level_max, level_inds)
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
    S.weights[idx] != 0.0 && error()
    S.totweight[] += weight
    S.totvalues[] += 1
    level_raw = ceil(Int, log2(weight))
    createlevel!(S, level_raw)
    level = getlevel(first(S.level_inds), weight)
    S.level_weights[level] += weight
    push!(S.level_buckets[level], idx)
	S.weights[idx] = weight
    return S
end

function Base.append!(S::DynamicSampler, e::Tuple)
    inds, weights = e
    nlevels = zeros(Int, length(S.level_buckets))
    sumweights = 0.0
    levs = zeros(Int16, length(inds))
    for (i, w) in enumerate(weights)
        S.weights[i] != 0.0 && error()
        level_raw = ceil(Int, log2(w))
        createlevel!(S, level_raw, nlevels)
        level = level_raw - first(S.level_inds) + 1
        nlevels[level] += 1
        S.level_weights[level] += w
        S.weights[i] = w
        sumweights += w
        S.totvalues[] += 1
        levs[i] = level_raw
    end
    S.totweight[] += sumweights
    for (i, bucket) in enumerate(S.level_buckets)
        resize!(bucket, length(bucket) + nlevels[i])
        nlevels[i] = length(bucket)
    end
    for (i, id) in enumerate(inds)
        level_raw = levs[i]
        level = level_raw - first(S.level_inds) + 1
        bucket = S.level_buckets[level]
        bucket[nlevels[level]] = i
        nlevels[level] -= 1
    end
    return S
end

function Base.rand(S::DynamicSampler; info = false)
    local level, idx, idx_in_level, weight  
    # Sample a level using the CDF method
    u = rand(S.rng) * S.totweight[]
    cumulative_weight = 0.0
    for i in eachindex(S.level_weights)
        cumulative_weight += S.level_weights[i]
        level = i
        u < cumulative_weight && break
    end         
    # Now sample within the level using rejection sampling
    bucket = S.level_buckets[level]
    level_size = length(bucket)
    if level_size == 0
        notempty = findall(b -> length(b) > 0, S.level_buckets)
        level = rand(S.rng, notempty)
        bucket = S.level_buckets[level]
        level_size = length(bucket)
    end
    level_max = S.level_max[level]
    while true
        idx_in_level = rand(S.rng, 1:level_size)
        idx = bucket[idx_in_level]
        weight = S.weights[idx]
        rand(S.rng) * level_max <= weight && break
    end     
    return info == false ? idx : DynamicIndex(idx, weight, level, idx_in_level)
end
    
function Base.deleteat!(S::DynamicSampler, idx)
    weight = S.weights[idx]
    level = getlevel(first(S.level_inds), weight)
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
    S.totvalues[] -= 1
    S.totweight[] -= weight
    S.level_weights[level] -= weight
    # Swap with last element for efficent delete
    bucket = S.level_buckets[level]
    bucket[idx_in_level], bucket[end] = bucket[end], bucket[idx_in_level]
    pop!(bucket)
end

Base.isempty(S::DynamicSampler) = S.totvalues[] == 0

function createlevel!(S, level_w, nlevels=nothing)
    if isempty(S.level_max)
        push!(S.level_max, 2.0^level_w)
        push!(S.level_inds, level_w)
    elseif level_w > S.level_inds[end]
        for i in S.level_inds[end]+1:level_w
            push!(S.level_max, 2.0^i)
            push!(S.level_inds, i)
            push!(S.level_buckets, Int[])
            push!(S.level_weights, 0.0)
            !isnothing(nlevels) && push!(nlevels, 0)
        end
    elseif level_w < S.level_inds[begin]
        for i in S.level_inds[begin]-1:-1:level_w-1
            pushfirst!(S.level_max, 2.0^i)
            pushfirst!(S.level_inds, i)
            pushfirst!(S.level_buckets, Int[])
            pushfirst!(S.level_weights, 0.0)
            !isnothing(nlevels) && pushfirst!(nlevels, 0)
        end
    end
    return S
end

getlevel(minlevel, weight) = ceil(Int, log2(weight)) - minlevel + 1

end
