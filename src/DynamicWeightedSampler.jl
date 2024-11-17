
# Inspired by https://www.aarondefazio.com/tangentially/?p=58
# and https://github.com/adefazio/sampler

mutable struct DynamicInfo
    totvalues::Int
    totweight::Float64
    idx::Int
    weight::Float64
    level::Int16
    idx_in_level::Int
end

struct DynamicSampler{R}
    rng::R
    info::DynamicInfo
    weights::Vector{Float64}
    level_weights::Vector{Float64}
    level_buckets::Vector{Vector{Int}}
    level_max::Vector{Float64}
    level_inds::Vector{Int16}
    inds_to_level::Vector{Int}
end

DynamicSampler() = DynamicSampler(Random.default_rng())
function DynamicSampler(rng)
    weights = Float64[]
    level_weights = Float64[0.0]
    level_buckets = [Int[],]
    level_max = Float64[0.0]
    level_inds = Int16[]
    inds_to_level = Int[]
    return DynamicSampler(rng, DynamicInfo(0, 0.0, 0, 0.0, Int16(0), 0),
        weights, level_weights, level_buckets, level_max, level_inds, 
        inds_to_level)
end

Base.sizehint!(sp::DynamicSampler, N) = resize_w!(sp, N)

@inline function Base.push!(sp::DynamicSampler, idx, w)
    sp.info.idx = 0
    resize_w!(sp, idx)
    sp.weights[idx] != 0.0 && error()
    sp.info.totweight += w
    sp.info.totvalues += 1
    level_raw = fast_ceil_log2(w)
    createlevel!(sp, level_raw)
    level = level_raw - Int(sp.level_inds[1]) + 1
    prev_w_max = sp.level_max[level]
    sp.level_max[level] = w > prev_w_max ? w : prev_w_max
    sp.level_weights[level] += w
    bucket = sp.level_buckets[level]
    push!(bucket, idx)
    !isempty(sp.inds_to_level) && (sp.inds_to_level[idx] = length(bucket))
    sp.weights[idx] = w
    return sp
end

function Base.append!(sp::DynamicSampler, inds, weights)
    sp.info.idx = 0
    nlevels = zeros(Int, length(sp.level_buckets))
    sumweights = 0.0
    sumvalues = 0
    levs = Vector{Int16}(undef, length(inds))
    resize_w!(sp, maximum(inds))
    @inbounds for (i, w) in enumerate(weights)
        sp.weights[i] != 0.0 && error()
        level_raw = fast_ceil_log2(w)
        createlevel!(sp, level_raw, nlevels)
        level = level_raw - Int(sp.level_inds[1]) + 1
        prev_w_max = sp.level_max[level]
        sp.level_max[level] = w > prev_w_max ? w : prev_w_max
        nlevels[level] += 1
        sp.level_weights[level] += w
        sp.weights[i] = w
        sumweights += w
        sumvalues += 1
        levs[i] = Int16(level_raw)
    end
    sp.info.totweight += sumweights
    sp.info.totvalues += sumvalues
    @inbounds @simd for i in 1:length(sp.level_buckets)
        bucket = sp.level_buckets[i]
        bucket_length = length(bucket)
        resize!(bucket, bucket_length + nlevels[i])
        nlevels[i] = bucket_length
    end
    offset_level = Int(sp.level_inds[1]) - 1
    @inbounds for (i, id) in enumerate(inds)
        level = Int(levs[i]) - offset_level
        bucket = sp.level_buckets[level]
        nlevels[level] += 1
        bucket[nlevels[level]] = id
        !isempty(sp.inds_to_level) && (sp.inds_to_level[idx] = nlevels[level])
    end
    return sp
end

Base.rand(sp::DynamicSampler, n::Integer) = [rand(sp) for _ in 1:n]
@inline function Base.rand(sp::DynamicSampler)
    local level, idx, idx_in_level, weight
    # Sample a level using the CDF method
    u = rand(sp.rng) * sp.info.totweight
    cumulative_weight = 0.0
    level = length(sp.level_weights)
    for i in reverse(eachindex(sp.level_weights))
        cumulative_weight += sp.level_weights[i]
        if u < cumulative_weight
            level = i
            break
        end
    end
    bucket = sp.level_buckets[level]
    level_size = length(bucket)
    if level_size == 0
        n_notempty = sum(length(b) > 0 for b in sp.level_buckets)
        rand_notempty = rand(sp.rng, 1:n_notempty)
        notempty = ((i,b) for (i,b) in enumerate(sp.level_buckets) if length(b) > 0)
        level, bucket = first(Iterators.drop(notempty, rand_notempty-1))
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
    sp.info.idx = idx
    sp.info.weight = weight
    sp.info.level = Int16(level)
    sp.info.idx_in_level = idx_in_level
    return idx
end

function Base.delete!(sp::DynamicSampler, indices::Union{UnitRange, Vector{<:Integer}})
    for i in indices
        delete!(sp, i)
    end
    return sp
end
@inline function Base.delete!(sp::DynamicSampler, idx::Integer)
    if sp.info.idx == idx
        _delete!(sp, idx, sp.info.weight, Int(sp.info.level), sp.info.idx_in_level)
    else
        weight = sp.weights[idx]
        level = fast_ceil_log2(weight) - Int(sp.level_inds[1]) + 1
        if isempty(sp.inds_to_level)
            resize!(sp.inds_to_level, length(sp.weights))
            for bucket in sp.level_buckets
                for (i, idx) in enumerate(bucket)
                    sp.inds_to_level[idx] = i
                end
            end
        end
        idx_in_level = sp.inds_to_level[idx]
        _delete!(sp, idx, weight, level, idx_in_level)
    end
    return sp
end
@inline function _delete!(sp, idx, weight, level, idx_in_level)
    sp.info.idx = 0
    sp.weights[idx] = 0.0
    sp.info.totvalues -= 1
    sp.info.totweight -= weight
    sp.level_weights[level] -= weight
    bucket = sp.level_buckets[level]
    bucket[idx_in_level], bucket[end] = bucket[end], bucket[idx_in_level]
    if !isempty(sp.inds_to_level)
        idx_other = bucket[idx_in_level]
        sp.inds_to_level[idx_other] = idx_in_level
    end
    pop!(bucket)
end

function Base.empty!(sp::DynamicSampler)
    @inbounds for (i, bucket) in enumerate(sp.level_buckets)
        sp.level_max[i] = 0.0
        sp.level_weights[i] = 0.0
        @simd for j in bucket
            sp.weights[j] = 0.0
            sp.inds_to_level[j] = 0
        end
        empty!(bucket)
    end
    sp.info.totvalues = 0
    sp.info.totweight = 0.0
    sp.info.idx = 0
end

Base.in(idx::Integer, sp::DynamicSampler) = sp.weights[idx] != 0.0
Base.isempty(sp::DynamicSampler) = sp.info.totvalues == 0

allinds(sp::DynamicSampler) = reduce(vcat, sp.level_buckets)

@inline function Base.setindex!(sp::DynamicSampler, idx, new_weight)
    delete!(sp, idx)
    push!(sp, idx, new_weight)
    return sp
end

function Base.show(io::IO, mime::MIME"text/plain", sp::DynamicSampler)
    inds = allinds(sp)
    print("DynamicSampler(indices = $(inds), weights = $(sp.weights[inds]))")
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

fast_ceil_log2(x) = ceil(Int, log2(x))
function fast_ceil_log2(x::Float64)
    x <= 0.0 && error()
    x_bits = reinterpret(UInt64, x)
    # Extract the exponent part (bits 52-62 in IEEE 754 double-precision)
    exponent = Int((x_bits >> 52) & 0x7FF) - 1023
    return exponent + Int((x_bits & 0xFFFFFFFFFFFFF) != 0)
end
function fast_ceil_log2(x::Float32)
    x <= 0.0f0 && error()
    x_bits = reinterpret(UInt32, x)
    # Extract the exponent part (bits 23-30 in IEEE 754 single-precision)
    exponent = Int((x_bits >> 23) & 0xFF) - 127
    return exponent + Int((x_bits & 0x7FFFFF) != 0)
end
function fast_ceil_log2(x::Float16)
    x <= Float16(0.0) && error()
    x_bits = reinterpret(UInt16, x)
    # Extract the exponent part (bits 10-14 in IEEE 754 half-precision)
    exponent = Int((x_bits >> 10) & 0x1F) - 15
    return exponent + Int((x_bits & 0x3FF) != 0)
end
