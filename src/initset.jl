# =============================================================================
#  initset.jl -- population initialisation
# =============================================================================

"""
    init_population!(X, lb, ub, p0, rng; method = :uniform)

Fill the internal `dim × nagents` matrix `X` with a starting population.

`method`
* `:uniform` -- independent uniform draws (identical to the original code)
* `:lhs`     -- Latin-hypercube sampling: every dimension is stratified into
                `n` equal bins and each bin is used exactly once.  Gives a much
                more even coverage of the box for the same number of agents and
                usually buys 1-2 iterations of convergence for free.
* `:center`  -- uniform, but agent 1 is placed at the centre of the box

If `p0` has the right length it is copied into agent 1 (as before).
"""
function init_population!(X::Matrix{Float64}, lb::Vector{Float64}, ub::Vector{Float64},
    p0, rng::AbstractRNG; method::Symbol=:uniform)

    dim, n = size(X)

    if method === :lhs
        perm = Vector{Int}(undef, n)
        @inbounds for i in 1:dim
            span = (ub[i] - lb[i]) / n
            @simd for j in 1:n
                perm[j] = j
            end
            shuffle!(rng, perm)
            for j in 1:n
                X[i, j] = lb[i] + span * (perm[j] - 1 + rand(rng))
            end
        end
    else
        @inbounds for j in 1:n, i in 1:dim
            X[i, j] = lb[i] + rand(rng) * (ub[i] - lb[i])
        end
        if method === :center && n >= 1
            @inbounds for i in 1:dim
                X[i, 1] = 0.5 * (lb[i] + ub[i])
            end
        end
    end

    if !isempty(p0)
        if length(p0) == dim
            @inbounds for i in 1:dim
                X[i, 1] = Float64(p0[i])
            end
        else
            @warn "p0_optional length ($(length(p0))) != problem dimension ($dim); ignoring it."
        end
    end
    return X
end

"""
    opposite_population!(Xo, X, lb, ub)

Opposition-Based Learning (Tizhoosh 2005) reflection `x' = lb + ub - x`.
Evaluating the opposite population once at start-up and keeping the better half
is one of the cheapest known accelerators for population metaheuristics: it
costs a single extra generation of function evaluations and typically removes
several percent of the total budget needed to reach a given accuracy.
"""
function opposite_population!(Xo::Matrix{Float64}, X::Matrix{Float64},
    lb::Vector{Float64}, ub::Vector{Float64})
    dim, n = size(X)
    @inbounds for j in 1:n
        @simd for i in 1:dim
            Xo[i, j] = lb[i] + ub[i] - X[i, j]
        end
    end
    return Xo
end

"""
    initialization(SearchAgents_no, dim, ub, lb, p0_optional) -> Matrix

Legacy public helper (kept for backwards compatibility).  Returns the
`SearchAgents_no × dim` matrix expected by user scripts written against
MPAOP ≤ 0.2.
"""
function initialization(SearchAgents_no::Int64, dim::Int64, ub::Vector{Float64},
    lb::Vector{Float64}, p0_optional)
    X = Matrix{Float64}(undef, dim, SearchAgents_no)
    init_population!(X, lb, ub, p0_optional, Random.default_rng())
    P = Matrix{Float64}(undef, SearchAgents_no, dim)
    return transpose_into!(P, X)
end
