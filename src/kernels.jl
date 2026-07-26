# =============================================================================
#  kernels.jl -- allocation-free numerical kernels used by MOMPA
#
#  MEMORY LAYOUT NOTE
#  ------------------
#  Julia stores matrices column-major.  Internally every population is kept
#  **transposed** with respect to the public API:
#
#        internal   X :: Matrix{Float64}   size (dim, nagents)
#        public     P :: Matrix{Float64}   size (nagents, dim)
#
#  so that one agent == one contiguous column.  This makes
#    * copying an agent into the objective-function buffer a plain `memcpy`,
#    * every update loop stride-1,
#    * MPI scatter/gather of agent blocks contiguous (no packing needed).
#  Conversion to the public layout happens only at the boundaries.
# =============================================================================

"""
    sanitize(v) -> Float64

Map `NaN` to `Inf` so that a misbehaving objective function can never poison
`findmin`/dominance comparisons.  Everything else is passed through.
"""
@inline sanitize(v::Real) = (x = Float64(v); isnan(x) ? Inf : x)

# Catches the commonest mistake: a vector-valued objective run with the default
# `num_objectives = 1`.
sanitize(v) = throw(ArgumentError(
    "the objective returned a $(typeof(v)) but `num_objectives = 1` expects a " *
    "single real number. For multi-objective optimisation pass " *
    "`num_objectives = $(applicable(length, v) ? length(v) : "<k>")`."))

# --- box constraint --------------------------------------------------------

"""
    clamp_cols!(X, lb, ub)

In-place box projection of a `dim × nagents` matrix.  Replaces the allocating
`X = clamp.(X, lb', ub')` of the original implementation.
"""
function clamp_cols!(X::AbstractMatrix{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
    dim = size(X, 1)
    @inbounds for j in axes(X, 2)
        @simd for i in 1:dim
            x = X[i, j]
            X[i, j] = ifelse(x < lb[i], lb[i], ifelse(x > ub[i], ub[i], x))
        end
    end
    return X
end

# --- Lévy flight -----------------------------------------------------------

"""
    levy_sigma(beta) -> Float64

Scale parameter of the Mantegna Lévy generator.  Depends only on `beta`, so it
is computed once per run instead of once per iteration.
"""
function levy_sigma(beta::Float64)
    num = gamma(1 + beta) * sin(pi * beta / 2)
    den = gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)
    return (num / den)^(1 / beta)
end

"""
    levy_step(rng, sigma, inv_beta) -> Float64

Single Lévy-distributed increment.  Drawing element-by-element inside the
movement kernel removes the two `n × dim` temporaries the original `levy`
allocated on **every** iteration.
"""
@inline function levy_step(rng::AbstractRNG, sigma::Float64, inv_beta::Float64)
    u = randn(rng) * sigma
    v = randn(rng)
    return u / abs(v)^inv_beta
end

"""
    levy!(Z, rng; beta = 1.5, scale = 1.0) -> Z

Fill `Z` in place with Lévy increments (allocation free).
"""
function levy!(Z::AbstractArray{Float64}, rng::AbstractRNG=Random.default_rng();
    beta::Float64=1.5, scale::Float64=1.0)
    sigma = levy_sigma(beta)
    invb = 1 / beta
    @inbounds for i in eachindex(Z)
        Z[i] = scale * levy_step(rng, sigma, invb)
    end
    return Z
end

"""
    levy(n, m, beta) -> Matrix{Float64}

Backwards-compatible allocating version (kept for user code / legacy scripts).
"""
function levy(n::Integer, m::Integer, beta::Float64)
    Z = Matrix{Float64}(undef, n, m)
    levy!(Z, Random.default_rng(); beta=beta)
    return Z
end

# --- MPA movement ----------------------------------------------------------

# Elite accessor: single global leader (one objective) or per-agent leader (many).
@inline elite_at(E::Vector{Float64}, i::Int, ::Int) = @inbounds E[i]
@inline elite_at(E::Matrix{Float64}, i::Int, j::Int) = @inbounds E[i, j]

"""
    mpa_move!(Xnew, X, E, stage, P, CF, wfac, rng, sigma, inv_beta)

The three-phase Marine-Predators movement operator, fused into one
allocation-free pass over the population.

* `Xnew`, `X` : `dim × n`, may alias-free double buffers
* `E`         : leader(s); `Vector` (one objective) or `dim × n` `Matrix` (many)
* `stage`     : 1 = Brownian / high-velocity, 2 = mixed, 3 = Lévy / low-velocity
* `wfac`      : NMPA inertia weight (`1.0` for `:standard_mpa`)

Random numbers are drawn on the fly; only the draws a given branch actually
consumes are generated (the reference implementation always materialised three
full `n × dim` matrices per half-iteration).
"""
function mpa_move!(Xnew::Matrix{Float64}, X::Matrix{Float64}, E::EA,
    stage::Int, P::Float64, CF::Float64, wfac::Float64,
    rng::AbstractRNG, sigma::Float64, inv_beta::Float64) where {EA}

    dim, n = size(X)
    half = n ÷ 2
    LS = 0.05                     # Lévy scaling used by the reference code

    @inbounds if stage == 1
        for j in 1:n, i in 1:dim
            rb = randn(rng)
            r = rand(rng)
            xp = X[i, j]
            s = rb * (elite_at(E, i, j) - rb * xp)
            Xnew[i, j] = xp + P * r * s
        end
    elseif stage == 2
        for j in 1:n
            if j <= half
                for i in 1:dim
                    rl = LS * levy_step(rng, sigma, inv_beta)
                    r = rand(rng)
                    xp = X[i, j]
                    s = rl * (elite_at(E, i, j) - rl * xp)
                    Xnew[i, j] = wfac * xp + P * r * s
                end
            else
                for i in 1:dim
                    rb = randn(rng)
                    xe = elite_at(E, i, j)
                    s = rb * (rb * xe - X[i, j])
                    Xnew[i, j] = wfac * xe + P * CF * s
                end
            end
        end
    else
        for j in 1:n, i in 1:dim
            rl = LS * levy_step(rng, sigma, inv_beta)
            xe = elite_at(E, i, j)
            s = rl * (rl * xe - X[i, j])
            Xnew[i, j] = xe + P * CF * s
        end
    end
    return Xnew
end

"""
    fads!(Xnew, X, lb, ub, FADs, CF, rng, perm1, perm2)

Fish-Aggregating-Devices / eddy-formation effect.  Both branches write into
`Xnew` (which may be a scratch buffer that is subsequently swapped with `X`),
so no snapshot copy of the population is needed.
"""
function fads!(Xnew::Matrix{Float64}, X::Matrix{Float64},
    lb::Vector{Float64}, ub::Vector{Float64},
    FADs::Float64, CF::Float64, rng::AbstractRNG,
    perm1::Vector{Int}, perm2::Vector{Int})

    dim, n = size(X)
    @inbounds if rand(rng) < FADs
        for j in 1:n, i in 1:dim
            x = X[i, j]
            if rand(rng) < FADs
                x += CF * (lb[i] + rand(rng) * (ub[i] - lb[i]))
            end
            Xnew[i, j] = x
        end
    else
        r = rand(rng)
        s = FADs * (1 - r) + r
        shuffle!(rng, perm1)
        shuffle!(rng, perm2)
        for j in 1:n
            j1 = perm1[j]
            j2 = perm2[j]
            @simd for i in 1:dim
                Xnew[i, j] = X[i, j] + s * (X[i, j1] - X[i, j2])
            end
        end
    end
    return Xnew
end

# --- misc helpers ----------------------------------------------------------

"""
    chunk_ranges(n, k) -> Vector{UnitRange{Int}}

Split `1:n` into at most `k` balanced contiguous ranges (empty ranges dropped).
"""
function chunk_ranges(n::Int, k::Int)
    k = max(1, min(k, n))
    out = Vector{UnitRange{Int}}(undef, k)
    lo = 1
    @inbounds for t in 1:k
        len = cld(n - lo + 1, k - t + 1)
        out[t] = lo:(lo+len-1)
        lo += len
    end
    return out
end

"""
    block_counts(n, nprocs) -> (counts, displs)

Balanced contiguous partition of `n` items over `nprocs` MPI ranks.  Used with
`MPI.Allgatherv!`, which removes the population padding (and therefore the
wasted objective evaluations) of the original MPI path.
"""
function block_counts(n::Int, nprocs::Int)
    counts = Vector{Cint}(undef, nprocs)
    displs = Vector{Cint}(undef, nprocs)
    lo = 0
    @inbounds for r in 1:nprocs
        len = cld(n - lo, nprocs - r + 1)
        counts[r] = Cint(len)
        displs[r] = Cint(lo)
        lo += len
    end
    return counts, displs
end

"""
    transpose_into!(dest, src)

`dest[j, i] = src[i, j]`; converts between the internal (`dim × n`) and public
(`n × dim`) layouts without allocating.
"""
function transpose_into!(dest::AbstractMatrix{Float64}, src::AbstractMatrix{Float64})
    @inbounds for j in axes(src, 2), i in axes(src, 1)
        dest[j, i] = src[i, j]
    end
    return dest
end

"""
    make_rng(rng, seed) -> AbstractRNG

Resolve the user's RNG options.  A `seed` makes an entire run bit-reproducible
(serial and MPI alike, because all stochastic decisions are taken on rank 0).
"""
function make_rng(rng::Union{Nothing,AbstractRNG}, seed::Union{Nothing,Integer})
    rng !== nothing && return rng
    seed !== nothing && return Xoshiro(seed)
    return Random.default_rng()
end
