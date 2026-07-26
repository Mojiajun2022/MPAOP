# =============================================================================
#  metrics.jl -- quality indicators for multi-objective results
#
#  All functions take the public orientation: one solution per ROW.
#  Minimisation is assumed throughout.
# =============================================================================

"""
    pareto_filter(objs) -> BitVector

`true` for every row of `objs` (`n × nobj`) that is not dominated by another
row.  Handy for post-processing an archive obtained with `mode = :fronts`, or
for merging the results of several independent runs.
"""
function pareto_filter(objs::AbstractMatrix{<:Real})
    n, M = size(objs)
    keep = trues(n)
    @inbounds for i in 1:n
        keep[i] || continue
        for j in 1:n
            i == j && continue
            keep[j] || continue
            # does j dominate i ?
            strict = false
            dom = true
            for k in 1:M
                a = objs[j, k]
                b = objs[i, k]
                if a > b
                    dom = false
                    break
                elseif a < b
                    strict = true
                end
            end
            if dom && strict
                keep[i] = false
                break
            end
        end
    end
    return keep
end

function _hv_recursive(pts::Vector{Vector{Float64}}, ref::Vector{Float64}, M::Int)
    isempty(pts) && return 0.0
    if M == 1
        best = Inf
        for p in pts
            p[1] < best && (best = p[1])
        end
        return max(0.0, ref[1] - best)
    elseif M == 2
        sort!(pts, by=p -> (p[1], p[2]))
        hv = 0.0
        top = ref[2]
        for p in pts
            if p[2] < top
                hv += (ref[1] - p[1]) * (top - p[2])
                top = p[2]
            end
        end
        return hv
    else
        sort!(pts, by=p -> p[M])
        hv = 0.0
        n = length(pts)
        refl = ref[1:M-1]
        for i in 1:n
            upper = i < n ? pts[i+1][M] : ref[M]
            h = upper - pts[i][M]
            h <= 0 && continue
            slice = [p[1:M-1] for p in view(pts, 1:i)]
            hv += h * _hv_recursive(slice, refl, M - 1)
        end
        return hv
    end
end

function _hv_montecarlo(pts::Vector{Vector{Float64}}, ref::Vector{Float64}, M::Int,
    samples::Int, rng::AbstractRNG)
    lo = fill(Inf, M)
    for p in pts, k in 1:M
        p[k] < lo[k] && (lo[k] = p[k])
    end
    box = 1.0
    for k in 1:M
        box *= (ref[k] - lo[k])
    end
    box <= 0 && return 0.0
    x = Vector{Float64}(undef, M)
    hits = 0
    for _ in 1:samples
        for k in 1:M
            x[k] = lo[k] + rand(rng) * (ref[k] - lo[k])
        end
        for p in pts
            covered = true
            for k in 1:M
                if p[k] > x[k]
                    covered = false
                    break
                end
            end
            if covered
                hits += 1
                break
            end
        end
    end
    return box * hits / samples
end

"""
    hypervolume(objs, ref; method = :auto, samples = 200_000, rng = default_rng())

Hypervolume dominated by the rows of `objs` and bounded by the reference point
`ref` (which must be dominated by the solutions).  Larger is better; it is the
only common indicator that is strictly Pareto-compliant, so it is the right
number to report when comparing runs.

`method`
* `:exact` -- dimension-sweep recursion; exact for any `nobj`, cost grows like
  `O(n^{nobj-1})`
* `:mc`    -- Monte-Carlo estimate, cost `O(samples · n · nobj)`
* `:auto`  -- exact for `nobj ≤ 3` (or small sets), Monte-Carlo otherwise
"""
function hypervolume(objs::AbstractMatrix{<:Real}, ref::AbstractVector{<:Real};
    method::Symbol=:auto, samples::Int=200_000,
    rng::AbstractRNG=Random.default_rng())
    n, M = size(objs)
    length(ref) == M || throw(DimensionMismatch("ref has $(length(ref)) entries, objectives have $M"))
    n == 0 && return 0.0

    keep = pareto_filter(objs)
    pts = Vector{Vector{Float64}}()
    @inbounds for i in 1:n
        keep[i] || continue
        ok = true
        for k in 1:M
            if !(objs[i, k] < ref[k])
                ok = false
                break
            end
        end
        ok && push!(pts, Float64[objs[i, k] for k in 1:M])
    end
    isempty(pts) && return 0.0

    use_exact = method === :exact ||
                (method === :auto && (M <= 3 || length(pts) <= 20))
    return use_exact ? _hv_recursive(pts, Float64.(collect(ref)), M) :
           _hv_montecarlo(pts, Float64.(collect(ref)), M, samples, rng)
end

@inline function _min_dist(row::AbstractMatrix{<:Real}, i::Int, other::AbstractMatrix{<:Real}, M::Int)
    best = Inf
    @inbounds for j in axes(other, 1)
        d = 0.0
        for k in 1:M
            t = row[i, k] - other[j, k]
            d += t * t
        end
        d < best && (best = d)
    end
    return sqrt(best)
end

"""
    igd(objs, reference_front)

Inverted Generational Distance: mean distance from every point of the *true*
front to the closest obtained solution.  Smaller is better; measures
convergence **and** coverage.
"""
function igd(objs::AbstractMatrix{<:Real}, reference_front::AbstractMatrix{<:Real})
    M = size(objs, 2)
    nr = size(reference_front, 1)
    (nr == 0 || size(objs, 1) == 0) && return Inf
    s = 0.0
    @inbounds for i in 1:nr
        s += _min_dist(reference_front, i, objs, M)
    end
    return s / nr
end

"""
    gd(objs, reference_front)

Generational Distance: mean distance from every obtained solution to the
closest point of the true front.  Smaller is better; measures convergence only.
"""
function gd(objs::AbstractMatrix{<:Real}, reference_front::AbstractMatrix{<:Real})
    M = size(objs, 2)
    n = size(objs, 1)
    (n == 0 || size(reference_front, 1) == 0) && return Inf
    s = 0.0
    @inbounds for i in 1:n
        s += _min_dist(objs, i, reference_front, M)
    end
    return s / n
end

"""
    spacing_metric(objs)

Schott's spacing: standard deviation of the nearest-neighbour distances inside
the obtained front.  Smaller means a more uniform distribution.
"""
function spacing_metric(objs::AbstractMatrix{<:Real})
    n, M = size(objs)
    n < 2 && return 0.0
    d = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        best = Inf
        for j in 1:n
            i == j && continue
            s = 0.0
            for k in 1:M
                s += abs(objs[i, k] - objs[j, k])     # Manhattan, as in Schott 1995
            end
            s < best && (best = s)
        end
        d[i] = best
    end
    m = sum(d) / n
    return sqrt(sum(x -> (x - m)^2, d) / (n - 1))
end

"""
    max_spread(objs)

Zitzler's maximum spread -- the diagonal of the bounding box of the obtained
front.  Larger means a wider front.
"""
function max_spread(objs::AbstractMatrix{<:Real})
    n, M = size(objs)
    n == 0 && return 0.0
    s = 0.0
    @inbounds for k in 1:M
        lo = Inf
        hi = -Inf
        for i in 1:n
            v = objs[i, k]
            v < lo && (lo = v)
            v > hi && (hi = v)
        end
        s += (hi - lo)^2
    end
    return sqrt(s)
end
