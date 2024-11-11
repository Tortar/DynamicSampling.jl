# Based on https://www.aarondefazio.com/tangentially/?p=58
# and https://github.com/adefazio/sampler

module DynamicSampling

export DynamicSampler
export allvalues

using Random

struct DynamicSampler{R}
    rng::R
    totvalues::Base.RefValue{Int}
    totweight::Base.RefValue{Float64}
    weights::Vector{Float64}
    level_weights::Vector{Float64}
    level_buckets::Vector{Vector{Int}}
    level_max::Vector{Float64}
    level_inds::Vector{Int16}
    inds_to_level::Vector{Int}
end

DynamicSampler() = DynamicSampler(Random.default_rng())
function DynamicSampler(rng)
    totweight = 0.0
    weights = Float64[]
    level_weights = Float64[0.0]
    level_buckets = [Int[],]
    level_max = Float64[0.0]
    level_inds = Int16[]
    inds_to_level = Int[]
    return DynamicSampler(rng, Ref(0), Ref(totweight), weights, level_weights, 
        level_buckets, level_max, level_inds, inds_to_level)
end

struct IndexInfo
    idx::Int
    weight::Float64
    level::Int
    idx_in_level::Int
end
IndexInfo(idx, weight) = IndexInfo(idx, weight, 0, 0)

Base.sizehint!(sp::DynamicSampler, N) = resize_w!(sp, N)

@inline function Base.push!(sp::DynamicSampler, idx, weight)
    resize_w!(sp, idx)
    sp.weights[idx] != 0.0 && error()
    sp.totweight[] += weight
    sp.totvalues[] += 1
    level_raw = ceil(Int, log2(weight))
    createlevel!(sp, level_raw)
    level = level_raw - Int(first(sp.level_inds)) + 1
    sp.level_max[level] = max(sp.level_max[level], weight)
    sp.level_weights[level] += weight
    bucket = sp.level_buckets[level]
    push!(bucket, idx)
    !isempty(sp.inds_to_level) && (sp.inds_to_level[idx] = length(bucket))
    sp.weights[idx] = weight
    return sp
end

function Base.append!(sp::DynamicSampler, inds, weights)
    nlevels = zeros(Int, length(sp.level_buckets))
    sumweights = 0.0
    levs = zeros(Int16, length(inds))
    resize_w!(sp, maximum(inds))
    @inbounds for (i, w) in enumerate(weights)
        sp.weights[i] != 0.0 && error()
        level_raw = ceil(Int, log2(w))
        createlevel!(sp, level_raw, nlevels)
        level = level_raw - Int(first(sp.level_inds)) + 1
        sp.level_max[level] = max(sp.level_max[level], w)
        nlevels[level] += 1
        sp.level_weights[level] += w
        sp.weights[i] = w
        sumweights += w
        sp.totvalues[] += 1
        levs[i] = ceil(Int, level_raw)
    end
    sp.totweight[] += sumweights
    @inbounds for (i, bucket) in enumerate(sp.level_buckets)
        resize!(bucket, length(bucket) + nlevels[i])
        nlevels[i] = length(bucket)
    end
    @inbounds for (i, id) in enumerate(inds)
        level_raw = levs[i]
        level = level_raw - Int(first(sp.level_inds)) + 1
        bucket = sp.level_buckets[level]
        bucket[nlevels[level]] = id
        !isempty(sp.inds_to_level) && (sp.inds_to_level[idx] = nlevels[level])
        nlevels[level] -= 1
    end
    return sp
end

@inline Base.@constprop :aggressive function Base.rand(sp::DynamicSampler; info = false)
    local level, idx, idx_in_level, weight
    # Sample a level using the CDF method
    u = rand(sp.rng) * sp.totweight[]
    cumulative_weight = 0.0
    level = 1
    for i in eachindex(sp.level_weights)
        cumulative_weight += sp.level_weights[i]
        if u < cumulative_weight
            level = i
            break
        end
    end   
    bucket = sp.level_buckets[level]
    level_size = length(bucket)
    if level_size == 0
        notempty = findall(b -> length(b) > 0, sp.level_buckets)
        level = rand(sp.rng, notempty)
        bucket = sp.level_buckets[level]
        level_size = length(bucket)
    end
    level_max = sp.level_max[level]      
    # Now sample within the level using rejection sampling
    @inbounds while true
        idx_in_level = rand(sp.rng, 1:level_size)
        idx = bucket[idx_in_level]
        weight = sp.weights[idx]
        rand(sp.rng) * level_max <= weight && break
    end
    return info == false ? idx : IndexInfo(idx, weight, level, idx_in_level)
end

@inline function Base.deleteat!(sp::DynamicSampler, idx)
    weight = sp.weights[idx]
    level = getlevel(Int(first(sp.level_inds)), weight)
    if isempty(sp.inds_to_level)
        resize!(sp.inds_to_level, length(sp.weights))
        for bucket in sp.level_buckets
            for (i, idx) in enumerate(bucket)
                sp.inds_to_level[idx] = i
            end
        end
    end
    idx_in_level = sp.inds_to_level[idx]
    _deleteat!(sp, idx, weight, level, idx_in_level)
    return sp
end
@inline function Base.deleteat!(sp::DynamicSampler, e::IndexInfo)
    idx, weight, level, idx_in_level = e.idx, e.weight, e.level, e.idx_in_level
    _deleteat!(sp, idx, weight, level, idx_in_level)
    return sp
end
@inline function _deleteat!(sp, idx, weight, level, idx_in_level)
    sp.weights[idx] = 0.0
    sp.totvalues[] -= 1
    sp.totweight[] -= weight
    sp.level_weights[level] -= weight
    bucket = sp.level_buckets[level]
    bucket[idx_in_level], bucket[end] = bucket[end], bucket[idx_in_level]
    if !isempty(sp.inds_to_level)
        idx_other = bucket[idx_in_level]
        sp.inds_to_level[idx_other] = idx_in_level
    end
    pop!(bucket)
end

Base.isempty(sp::DynamicSampler) = sp.totvalues[] == 0

allvalues(sp::DynamicSampler) = reduce(vcat, sp.level_buckets)

function Base.show(io::IO, mime::MIME"text/plain", sp::DynamicSampler)
    inds = allvalues(sp)
    print("DynamicSampler(indices = $(inds), weights = $(sp.weights[inds]))")
end

function Base.show(io::IO, mime::MIME"text/plain", di::IndexInfo)
    print("IndexInfo(idx = $(di.idx), weight = $(di.weight))")
end

@inline function resize_w!(sp, N)
    N_curr = length(sp.weights)
    if N > N_curr
        resize!(sp.weights, N)
        !isempty(sp.inds_to_level) && resize!(sp.inds_to_level, N)
        @inbounds @simd for i in N_curr+1:N
            sp.weights[i] = 0.0
        end
    end
    return sp
end

@inline function createlevel!(sp, level_w, nlevels=nothing)
    if isempty(sp.level_inds)
        push!(sp.level_inds, Int16(level_w))
    elseif level_w > sp.level_inds[end]
        for i in sp.level_inds[end]+1:level_w
            push!(sp.level_max, 0.0)
            push!(sp.level_inds, Int16(i))
            push!(sp.level_buckets, Int[])
            push!(sp.level_weights, 0.0)
            !isnothing(nlevels) && push!(nlevels, 0)
        end
    elseif level_w < sp.level_inds[begin]
        for i in sp.level_inds[begin]-1:-1:level_w
            pushfirst!(sp.level_max, 0.0)
            pushfirst!(sp.level_inds, Int16(i))
            pushfirst!(sp.level_buckets, Int[])
            pushfirst!(sp.level_weights, 0.0)
            !isnothing(nlevels) && pushfirst!(nlevels, 0)
        end
    end
    return sp
end

getlevel(minlevel, weight) = ceil(Int, log2(weight)) - minlevel + 1

end
