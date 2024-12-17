
using Distributions

mutable struct DynamicInfo
    totvalues::Int
    totweight::Float64
    toterror::Float64
    idx::Int
    weight::Float64
    level::Int
    idx_in_level::Int
    level_min::Int
    level_max::Int
    reorder::Int
end

struct DynamicSampler{R}
    rng::R
    info::DynamicInfo
    weights_assigned::BitVector
    raw_levels::Vector{Int16}
    level_weights::Vector{Float64}
    level_werrors::Vector{Float64}
    level_buckets::Vector{Vector{Tuple{Int, Float64}}}
    level_max::Vector{Float64}
    order_level::Vector{Int}
    inds_to_level::Vector{Int}
end

DynamicSampler() = DynamicSampler(Random.default_rng())
function DynamicSampler(rng)
    weights_assigned = BitVector()
    raw_levels = Int16[]
    level_weights = Float64[0.0]
    level_werrors = Float64[0.0]
    order_level = Int[1]
    level_buckets = [Tuple{Int, Float64}[],]
    level_max = Float64[0.0]
    inds_to_level = Int[]
    return DynamicSampler(rng, DynamicInfo(0, 0.0, 0.0, 0, 0.0, 0, 0, typemax(Int), typemin(Int), 0),
        weights_assigned, raw_levels, level_weights, level_werrors, level_buckets, level_max, 
        order_level, inds_to_level)
end

Base.sizehint!(sp::DynamicSampler, N) = resize_w!(sp, N)

@inline function Base.push!(sp::DynamicSampler, idx, w)
    sp.info.idx = 0
    resize_w!(sp, idx)
    sp.weights_assigned[idx] != false && error(lazy"index $(idx) already in the sampler")
    sp.weights_assigned[idx] = true
    sp.info.totvalues += 1
    level_raw = exponent(w) + 1
    createlevel!(sp, level_raw)
    level = level_raw - sp.info.level_min + 1
    prev_w_max = sp.level_max[level]
    sp.level_max[level] = w > prev_w_max ? w : prev_w_max
    sp.info.totweight, sp.info.toterror = two_sum(sp.info.totweight, w + sp.info.toterror)
    sp.level_weights[level], sp.level_werrors[level] = two_sum(sp.level_weights[level], w + sp.level_werrors[level])
    bucket = sp.level_buckets[level]
    push!(bucket, (idx, w))
    if !isempty(sp.inds_to_level)
        sp.raw_levels[idx] = Int16(level_raw)
        sp.inds_to_level[idx] = length(bucket)
    end
    resize_levels!(sp)
    reorder_levels!(sp, 1)
    return sp
end

function Base.append!(sp::DynamicSampler, inds, weights)
    sp.info.idx = 0
    resize_w!(sp, maximum(inds))
    sumweights, sumerrors = sp.info.totweight, sp.info.toterror
    sumvalues = 0
    for (i, w) in zip(inds, weights)
        sp.weights_assigned[i] != false && error(lazy"index $(i) already in the sampler")
        sp.weights_assigned[i] = true
        sumvalues += 1
        level_raw = exponent(w) + 1
        createlevel!(sp, level_raw)
        level = level_raw - sp.info.level_min + 1
        prev_w_max = sp.level_max[level]
        sp.level_max[level] = w > prev_w_max ? w : prev_w_max
        sumweights, sumerrors = two_sum(sumweights, w + sumerrors)
        sp.level_weights[level], sp.level_werrors[level] = two_sum(sp.level_weights[level], w + sp.level_werrors[level])
        bucket = sp.level_buckets[level]
        push!(bucket, (i, w))
        if !isempty(sp.inds_to_level)
            sp.raw_levels[i] = Int16(level_raw)
            sp.inds_to_level[i] = length(bucket)
        end
    end
    sp.info.totweight = sumweights
    sp.info.toterror = sumerrors
    sp.info.totvalues += sumvalues
    resize_levels!(sp)
    reorder_levels!(sp, sumvalues)
    return sp
end
function Base.append!(sp::DynamicSampler, inds::Union{UnitRange, AbstractArray}, 
        weights::Union{UnitRange, AbstractArray})
    @assert length(inds) == length(weights)
    sp.info.idx = 0
    nlevels = zeros(Int, length(sp.level_buckets))
    sumweights, sumerrors = sp.info.totweight, sp.info.toterror
    sumvalues = 0
    levs = Vector{Int16}(undef, length(inds))
    resize_w!(sp, maximum(inds))
    @inbounds for (i, w) in enumerate(weights)
        level_raw = exponent(w) + 1
        createlevel!(sp, level_raw, nlevels)
        level = level_raw - sp.info.level_min + 1
        prev_w_max = sp.level_max[level]
        sp.level_max[level] = w > prev_w_max ? w : prev_w_max
        nlevels[level] += 1
        sumweights, sumerrors = two_sum(sumweights, w + sumerrors)
        sp.level_weights[level], sp.level_werrors[level] = two_sum(sp.level_weights[level], w + sp.level_werrors[level])
        sumvalues += 1
        levs[i] = Int16(level_raw)
        if !isempty(sp.inds_to_level)
            sp.raw_levels[i] = Int16(level_raw)
        end
    end
    sp.info.totweight = sumweights
    sp.info.toterror = sumerrors
    sp.info.totvalues += sumvalues
    @inbounds @simd for i in 1:length(sp.level_buckets)
        bucket = sp.level_buckets[i]
        bucket_length = length(bucket)
        resize!(bucket, bucket_length + nlevels[i])
        nlevels[i] = bucket_length
    end
    offset_level = sp.info.level_min - 1
    @inbounds for (i, id) in enumerate(inds)
        sp.weights_assigned[id] != false && error(lazy"index $(id) already in the sampler")
        sp.weights_assigned[id] = true
        level = Int(levs[i]) - offset_level
        bucket = sp.level_buckets[level]
        nlevels[level] += 1
        bucket[nlevels[level]] = (id, weights[i])
        if !isempty(sp.inds_to_level)
            sp.inds_to_level[i] = nlevels[level]
        end
    end
    resize_levels!(sp)
    reorder_levels!(sp, sumvalues)
    return sp
end

function Base.rand(sp::DynamicSampler, n::Integer)
    sp.info.totweight = sum(sp.level_weights)
    n_each = rand(sp.rng, Multinomial(n, sp.level_weights ./ sp.info.totweight))
    randinds = Vector{Int}(undef, n)
    q = 1
    for (i, k) in enumerate(n_each)
        bucket = sp.level_buckets[i]
        for _ in 1:k
            randinds[q] = extract_rand_idx(sp, i, bucket)[1]
            q += 1
        end
    end
    shuffle!(sp.rng, randinds)
    return randinds
end
@inline function Base.rand(sp::DynamicSampler)
    level, bucket = extract_rand_level(sp)
    idx, weight, level, idx_in_level = extract_rand_idx(sp, level, bucket)
    sp.info.idx = idx
    sp.info.weight = weight
    sp.info.level = level
    sp.info.idx_in_level = idx_in_level
    return idx
end

@inline function extract_rand_level(sp::DynamicSampler)
    u = rand(sp.rng) * sp.info.totweight
    cumulative_weight = 0.0
    level = length(sp.level_weights)
    level_idx = 1
    @inbounds for (i, j) in enumerate(Iterators.reverse(sp.order_level))
        cumulative_weight += sp.level_weights[j]
        if u < cumulative_weight
            level_idx = i
            level = j
            break
        end
    end
    bucket = sp.level_buckets[level]
    if isempty(bucket) || level_idx > 32
        for i in eachindex(sp.level_weights)
            bucket = sp.level_buckets[i]
            sp.level_weights[i] = isempty(bucket) ? 0.0 : sum(x[2] for x in bucket)
            sp.level_werrors[i] = 0.0
        end
        sp.info.totweight = sum(sp.level_weights)
        sp.info.toterror = 0.0
        sortperm!(sp.order_level, sp.level_weights)
        n_notempty = sum(length(b) > 0 for b in sp.level_buckets)
        rand_notempty = rand(sp.rng, 1:n_notempty)
        notempty = Iterators.filter(i -> !isempty(sp.level_buckets[i]), eachindex(sp.level_buckets))
        level = first(Iterators.drop(notempty, rand_notempty-1))
        bucket = sp.level_buckets[level]
    end
    return level, bucket
end

@inline function extract_rand_idx(sp, level, bucket)
    level_max = sp.level_max[level]      
    u = rand(sp.rng) * length(bucket)
    intu = unsafe_trunc(Int, u)
    fracu = u - intu
    idx_in_level = intu + 1
    idx, weight = bucket[idx_in_level]
    @inbounds while fracu * level_max > weight
        u = rand(sp.rng) * length(bucket)
        intu = unsafe_trunc(Int, u)
        fracu = u - intu
        idx_in_level = intu + 1
        idx, weight = bucket[idx_in_level]
    end
    return idx, weight, level, idx_in_level
end

function Base.delete!(sp::DynamicSampler, indices)
    for i in indices
        delete!(sp, i)
    end
    return sp
end
@inline function Base.delete!(sp::DynamicSampler, idx::Integer)
    if sp.info.idx == idx
        _delete!(sp, idx, sp.info.weight, sp.info.level, sp.info.idx_in_level)
    else
        if isempty(sp.inds_to_level)
            resize!(sp.raw_levels, length(sp.weights_assigned))
            resize!(sp.inds_to_level, length(sp.weights_assigned))
            for bucket in sp.level_buckets
                isempty(bucket) && continue
                lraw = exponent(first(bucket)[2]) + 1
                for (i, (idx, w)) in enumerate(bucket)
                    sp.inds_to_level[idx] = i
                    sp.raw_levels[idx] = lraw
                end
            end
        end
        level = sp.raw_levels[idx] - sp.info.level_min + 1
        idx_in_level = sp.inds_to_level[idx]
        weight = sp.level_buckets[level][idx_in_level][2]
        _delete!(sp, idx, weight, level, idx_in_level)
    end
    return sp
end
@inline function _delete!(sp, idx, weight, level, idx_in_level)
    sp.weights_assigned[idx] != true && error(lazy"index $(idx) is not in the sampler")
    sp.info.idx = 0
    sp.weights_assigned[idx] = false
    sp.info.totvalues -= 1
    tw = sp.info.totweight
    sp.info.totweight, sp.info.toterror = two_sum(sp.info.totweight, -weight + sp.info.toterror)
    sp.level_weights[level], sp.level_werrors[level] = two_sum(sp.level_weights[level], -weight + sp.level_werrors[level])
    bucket = sp.level_buckets[level]
    bucket[idx_in_level], bucket[end] = bucket[end], bucket[idx_in_level]
    if !isempty(sp.inds_to_level)
        idx_other, weight_other = bucket[idx_in_level]
        sp.inds_to_level[idx_other] = idx_in_level
    end
    pop!(bucket)
    reorder_levels!(sp, 1)
end

function Base.empty!(sp::DynamicSampler)
    @inbounds for (i, bucket) in enumerate(sp.level_buckets)
        sp.level_max[i] = 0.0
        sp.level_weights[i] = 0.0
        sp.level_werrors[i] = 0.0
        @simd for j in bucket
            sp.weights_assigned[j[1]] = false
            sp.inds_to_level[j[1]] = 0
        end
        empty!(bucket)
    end
    empty!(sp.inds_to_level)
    sp.info.totvalues = 0
    sp.info.totweight, sp.info.toterror = 0.0, 0.0
    sp.info.idx = 0
end

Base.in(idx::Integer, sp::DynamicSampler) = sp.weights_assigned[idx] == true
Base.isempty(sp::DynamicSampler) = sp.info.totvalues == 0

allinds(sp::DynamicSampler) = collect(Iterators.Flatten((Iterators.map(x -> x[1], b) for b in sp.level_buckets)))

@inline function Base.setindex!(sp::DynamicSampler, idx, new_weight)
    delete!(sp, idx)
    push!(sp, idx, new_weight)
    return sp
end

function Base.show(io::IO, mime::MIME"text/plain", sp::DynamicSampler)
    print("DynamicSampler($(reduce(vcat, sp.level_buckets)))")
end

@inline function resize_w!(sp, N)
    N_curr = length(sp.weights_assigned)
    if N > N_curr
        resize!(sp.weights_assigned, N)
        if !isempty(sp.inds_to_level)
            resize!(sp.inds_to_level, N)
            resize!(sp.raw_levels, N)
        end
        @inbounds @simd for i in N_curr+1:N
            sp.weights_assigned[i] = false
        end
    end
    return sp
end

@inline function createlevel!(sp, level_w, nlevels=nothing)
    if sp.info.level_min == typemax(Int)
        sp.info.level_min = level_w
        sp.info.level_max = level_w
    elseif level_w > sp.info.level_max
        for i in sp.info.level_max+1:level_w
            push!(sp.level_max, 0.0)
            push!(sp.level_buckets, Int[])
            push!(sp.level_weights, 0.0)
            push!(sp.level_werrors, 0.0)
            !isnothing(nlevels) && push!(nlevels, 0)
        end
        sp.info.level_max = level_w
    elseif level_w < sp.info.level_min
        for i in sp.info.level_min-1:-1:level_w
            pushfirst!(sp.level_max, 0.0)
            pushfirst!(sp.level_buckets, Int[])
            pushfirst!(sp.level_weights, 0.0)
            pushfirst!(sp.level_werrors, 0.0)
            !isnothing(nlevels) && pushfirst!(nlevels, 0)
        end
        sp.info.level_min = level_w
    end
    return sp
end

@inline function resize_levels!(sp)
    n1, n2 = length(sp.level_weights), length(sp.order_level)
    if n2 < n1
        resize!(sp.order_level, n1)
        @inbounds @simd for i in (n2+1):n1
            sp.order_level[i] = i
        end
    end
    return sp
end

@inline function reorder_levels!(sp, k)
    sp.info.reorder += k
    if sp.info.reorder > length(sp.level_weights)*50
        sp.info.reorder = 0
        sortperm!(sp.order_level, sp.level_weights)
    end
    return sp
end

# From ErrorFreeAritmethic.jl
@inline function two_sum(a::T, b::T) where {T<:Real}
    hi = a + b
    v  = hi - a
    lo = (a - (hi - v)) + (b - v)
    (hi, lo)
end
