# =============================================================================
#  MOfun.jl -- multi-objective machinery
#
#  Internally objectives live in a `nobj × nsolutions` matrix (one solution per
#  column) so that a dominance test walks contiguous memory.  All public
#  helpers keep the historical `nsolutions × nobj` orientation.
# =============================================================================

# --- dominance -------------------------------------------------------------

"""
    is_dominated(obj_a, obj_b) -> Bool

`true` when `obj_b` dominates `obj_a` (minimisation).  Semantics identical to
MPAOP ≤ 0.2.
"""
function is_dominated(obj_a::AbstractVector{Float64}, obj_b::AbstractVector{Float64})::Bool
    any_l = false
    @inbounds for i in eachindex(obj_a)
        obj_b[i] > obj_a[i] && return false
        obj_b[i] < obj_a[i] && (any_l = true)
    end
    return any_l
end

"Column `i` of `F` dominates column `j` (both `nobj × n`)."
@inline function _dominates(F::Matrix{Float64}, i::Int, j::Int, M::Int)
    strict = false
    @inbounds for k in 1:M
        a = F[k, i]
        b = F[k, j]
        a > b && return false
        a < b && (strict = true)
    end
    return strict
end

"Lexicographic order on columns -- the pre-sort that makes ENS-SS correct."
@inline function _lexless(F::Matrix{Float64}, i::Int, j::Int, M::Int)
    @inbounds for k in 1:M
        a = F[k, i]
        b = F[k, j]
        a < b && return true
        a > b && return false
    end
    return false
end

@inline function _cols_equal(F::Matrix{Float64}, i::Int, j::Int, M::Int)
    @inbounds for k in 1:M
        F[k, i] == F[k, j] || return false
    end
    return true
end

# --- non-dominated sorting -------------------------------------------------

"""
    nds!(fronts, ranks, F, n, M, order)

Efficient Non-dominated Sort, sequential-search variant (ENS-SS, Zhang et al.
2015).  Solutions are pre-sorted lexicographically, therefore a solution can
only be dominated by one that precedes it and each candidate is compared
against the *last* member of a front first -- the cheapest rejection.

Complexity drops from the textbook `O(M n²)` **with** `O(n²)` bookkeeping
allocations (the previous implementation built an `n`-element `Vector{Vector}`
plus one temporary row per comparison) to `O(M n √n)` in practice with zero
allocation after warm-up.
"""
function nds!(fronts::Vector{Vector{Int}}, ranks::Vector{Int},
    F::Matrix{Float64}, n::Int, M::Int, order::Vector{Int})

    for f in fronts
        empty!(f)
    end
    nf = 0

    resize!(order, n)
    @inbounds for i in 1:n
        order[i] = i
    end
    sort!(order, alg=QuickSort, lt=(a, b) -> _lexless(F, a, b, M))

    @inbounds for s in order
        placed = false
        for k in 1:nf
            f = fronts[k]
            dominated = false
            for t in length(f):-1:1
                if _dominates(F, f[t], s, M)
                    dominated = true
                    break
                end
            end
            if !dominated
                push!(f, s)
                ranks[s] = k
                placed = true
                break
            end
        end
        if !placed
            nf += 1
            if nf > length(fronts)
                push!(fronts, Int[])
            end
            push!(fronts[nf], s)
            ranks[s] = nf
        end
    end
    return nf
end

"""
    non_dominated_sort(objectives) -> (ranks, fronts)

Public wrapper.  `objectives` is `nsolutions × nobj`; `ranks[i]` is the front
index of solution `i` (1 = Pareto front) and `fronts[k]` lists its members.
"""
function non_dominated_sort(objectives::Matrix{Float64})
    n, M = size(objectives)
    n == 0 && return Int[], Vector{Vector{Int}}()
    F = Matrix{Float64}(undef, M, n)
    transpose_into!(F, objectives)          # (n × M) -> (M × n)
    ranks = zeros(Int, n)
    fronts = Vector{Vector{Int}}()
    order = Int[]
    nf = nds!(fronts, ranks, F, n, M, order)
    return ranks, fronts[1:nf]
end

# --- crowding distance -----------------------------------------------------

"""
    crowding!(cd, F, idx, nidx, M, perm)

NSGA-II crowding distance of the solutions listed in `idx[1:nidx]`
(columns of the `nobj × n` matrix `F`), written into `cd` at the *global*
positions.  Uses a preallocated permutation buffer instead of `sortperm`.
"""
function crowding!(cd::Vector{Float64}, F::Matrix{Float64},
    idx::AbstractVector{Int}, nidx::Int, M::Int, perm::Vector{Int})
    nidx == 0 && return cd
    @inbounds for t in 1:nidx
        cd[idx[t]] = 0.0
    end
    if nidx <= 2
        @inbounds for t in 1:nidx
            cd[idx[t]] = Inf
        end
        return cd
    end

    resize!(perm, nidx)
    @inbounds for m in 1:M
        for t in 1:nidx
            perm[t] = idx[t]
        end
        sort!(perm, alg=QuickSort, lt=(a, b) -> F[m, a] < F[m, b])

        lo = F[m, perm[1]]
        hi = F[m, perm[nidx]]
        cd[perm[1]] = Inf
        cd[perm[nidx]] = Inf
        hi == lo && continue
        inv_range = 1.0 / (hi - lo)
        for t in 2:(nidx-1)
            g = perm[t]
            if isfinite(cd[g])
                cd[g] += (F[m, perm[t+1]] - F[m, perm[t-1]]) * inv_range
            end
        end
    end
    return cd
end

"""
    calculate_crowding_distance!(crowding_distances, objectives, front_indices, num_obj)

Public wrapper (legacy signature): `objectives` is `nsolutions × nobj`.
"""
function calculate_crowding_distance!(crowding_distances::Vector{Float64},
    objectives::Matrix{Float64}, front_indices::Vector{Int}, num_obj::Int)
    n = size(objectives, 1)
    F = Matrix{Float64}(undef, num_obj, n)
    transpose_into!(F, objectives)
    crowding!(crowding_distances, F, front_indices, length(front_indices), num_obj, Int[])
    return crowding_distances
end

# =============================================================================
#  Archive
# =============================================================================

"""
    MOArchive

Growable, preallocated external archive.  Columns `1:n` of `prey` / `obj` hold
the current members.

`mode`
* `:pareto` -- the archive is the running **non-dominated set**: a candidate is
  inserted only if no member dominates it, and every member it dominates is
  evicted.  Overflow is resolved by dropping the most crowded solutions.  This
  is the behaviour of the published MO-MPA / MOPSO archives and guarantees the
  returned front contains no dominated point.
* `:fronts` -- legacy MPAOP ≤ 0.2 behaviour: archive ∪ population is sorted into
  fronts and filled front by front, so the archive can (and usually does)
  contain dominated solutions once the first front is smaller than the archive.
"""
mutable struct MOArchive
    prey::Matrix{Float64}     # dim × cap
    obj::Matrix{Float64}      # M   × cap
    cd::Vector{Float64}
    n::Int
    maxsize::Int
    dim::Int
    M::Int
    mode::Symbol
    ranks::Vector{Int}
    fronts::Vector{Vector{Int}}
    order::Vector{Int}
    perm::Vector{Int}
    keep::Vector{Int}
end

function MOArchive(dim::Int, M::Int, maxsize::Int, npop::Int; mode::Symbol=:pareto)
    cap = max(maxsize + npop, 8)
    MOArchive(Matrix{Float64}(undef, dim, cap), Matrix{Float64}(undef, M, cap),
        zeros(Float64, cap), 0, maxsize, dim, M, mode,
        zeros(Int, cap), Vector{Vector{Int}}(), Int[], Int[], Int[])
end

Base.length(A::MOArchive) = A.n
capacity(A::MOArchive) = size(A.prey, 2)

function ensure_capacity!(A::MOArchive, need::Int)
    cap = capacity(A)
    need <= cap && return A
    newcap = max(need, 2cap)
    P = Matrix{Float64}(undef, A.dim, newcap)
    O = Matrix{Float64}(undef, A.M, newcap)
    @inbounds copyto!(view(P, :, 1:A.n), view(A.prey, :, 1:A.n))
    @inbounds copyto!(view(O, :, 1:A.n), view(A.obj, :, 1:A.n))
    A.prey = P
    A.obj = O
    resize!(A.cd, newcap)
    resize!(A.ranks, newcap)
    return A
end

@inline function copy_member!(A::MOArchive, dst::Int, src::Int)
    dst == src && return
    @inbounds for k in 1:A.dim
        A.prey[k, dst] = A.prey[k, src]
    end
    @inbounds for k in 1:A.M
        A.obj[k, dst] = A.obj[k, src]
    end
    return
end

@inline function store_member!(A::MOArchive, dst::Int, X::Matrix{Float64}, jx::Int,
    F::Matrix{Float64}, jf::Int)
    @inbounds for k in 1:A.dim
        A.prey[k, dst] = X[k, jx]
    end
    @inbounds for k in 1:A.M
        A.obj[k, dst] = F[k, jf]
    end
    return
end

"""
    insert_candidate!(A, X, jx, F, jf) -> Bool

Insert column `jx` of `X` (objectives in column `jf` of `F`) into the
non-dominated archive.  `O(n·M)` with early exit and no allocation.
"""
function insert_candidate!(A::MOArchive, X::Matrix{Float64}, jx::Int,
    F::Matrix{Float64}, jf::Int)
    M = A.M
    O = A.obj
    n = A.n

    # 1. rejected if some member is <= candidate in every objective
    #    (that member dominates it, or the two are identical)
    @inbounds for i in 1:n
        worse = false
        for k in 1:M
            if O[k, i] > F[k, jf]
                worse = true
                break
            end
        end
        worse || return false
    end

    # 2. evict every member the candidate dominates, compacting in place
    w = 0
    @inbounds for i in 1:n
        strict = false
        worse = false
        for k in 1:M
            a = F[k, jf]
            b = O[k, i]
            if a > b
                worse = true
                break
            elseif a < b
                strict = true
            end
        end
        if !(strict && !worse)
            w += 1
            copy_member!(A, w, i)
        end
    end

    w += 1
    store_member!(A, w, X, jx, F, jf)
    A.n = w
    return true
end

"Drop the most crowded solutions until `n == maxsize`."
function truncate_archive!(A::MOArchive)
    A.n <= A.maxsize && return A
    n = A.n
    resize!(A.order, n)
    @inbounds for i in 1:n
        A.order[i] = i
    end
    crowding!(A.cd, A.obj, A.order, n, A.M, A.perm)

    resize!(A.keep, n)
    @inbounds for i in 1:n
        A.keep[i] = i
    end
    # largest crowding distance first, then restore ascending order so the
    # in-place compaction below never overwrites a member it still needs
    partialsort!(A.keep, 1:A.maxsize, lt=(a, b) -> A.cd[a] > A.cd[b])
    resize!(A.keep, A.maxsize)
    sort!(A.keep, alg=QuickSort)
    @inbounds for w in 1:A.maxsize
        copy_member!(A, w, A.keep[w])
    end
    A.n = A.maxsize
    return A
end

"Refresh the cached crowding distances (used by `:crowding` elite selection)."
function refresh_crowding!(A::MOArchive)
    n = A.n
    n == 0 && return A
    resize!(A.order, n)
    @inbounds for i in 1:n
        A.order[i] = i
    end
    crowding!(A.cd, A.obj, A.order, n, A.M, A.perm)
    return A
end

"Legacy `:fronts` update: archive ∪ population sorted into fronts, filled in order."
function update_fronts!(A::MOArchive, X::Matrix{Float64}, F::Matrix{Float64})
    npop = size(X, 2)
    ensure_capacity!(A, A.n + npop)
    @inbounds for j in 1:npop
        store_member!(A, A.n + j, X, j, F, j)
    end
    A.n += npop

    # de-duplicate in objective space (lexicographic sweep)
    n = A.n
    resize!(A.order, n)
    @inbounds for i in 1:n
        A.order[i] = i
    end
    sort!(A.order, alg=QuickSort, lt=(a, b) -> _lexless(A.obj, a, b, A.M))
    resize!(A.keep, 0)
    @inbounds for t in 1:n
        i = A.order[t]
        if t == 1 || !_cols_equal(A.obj, i, A.order[t-1], A.M)
            push!(A.keep, i)
        end
    end
    sort!(A.keep, alg=QuickSort)
    @inbounds for (w, i) in enumerate(A.keep)
        copy_member!(A, w, i)
    end
    A.n = length(A.keep)

    n = A.n
    nf = nds!(A.fronts, A.ranks, A.obj, n, A.M, A.order)
    resize!(A.keep, 0)
    @inbounds for k in 1:nf
        f = A.fronts[k]
        if length(A.keep) + length(f) <= A.maxsize
            append!(A.keep, f)
        else
            crowding!(A.cd, A.obj, f, length(f), A.M, A.perm)
            room = A.maxsize - length(A.keep)
            room <= 0 && break
            sel = copy(f)
            partialsort!(sel, 1:room, lt=(a, b) -> A.cd[a] > A.cd[b])
            append!(A.keep, view(sel, 1:room))
            break
        end
    end
    sort!(A.keep, alg=QuickSort)
    @inbounds for (w, i) in enumerate(A.keep)
        copy_member!(A, w, i)
    end
    A.n = length(A.keep)
    return A
end

"""
    update!(A, X, F)

Fold a whole population (`X :: dim × npop`, `F :: nobj × npop`) into the
archive and refresh the cached crowding distances.
"""
function update!(A::MOArchive, X::Matrix{Float64}, F::Matrix{Float64})
    if A.mode === :fronts
        update_fronts!(A, X, F)
    else
        ensure_capacity!(A, A.n + size(X, 2))
        @inbounds for j in axes(X, 2)
            insert_candidate!(A, X, j, F, j)
        end
        truncate_archive!(A)
    end
    refresh_crowding!(A)
    return A
end

"""
    fill_elites!(E, A, rng, strategy)

Pick one leader per agent into the `dim × nagents` matrix `E`.

* `:crowding` -- binary tournament on the crowding distance.  Sparse regions of
  the front win more often, which spreads the population along the front and
  measurably improves the final IGD/spacing.
* `:random`   -- uniform sampling (MPAOP ≤ 0.2 behaviour).
"""
function fill_elites!(E::Matrix{Float64}, A::MOArchive, rng::AbstractRNG, strategy::Symbol)
    dim, nag = size(E)
    if A.n == 0
        fill!(E, 0.0)
        return E
    end
    n = A.n
    P = A.prey
    if strategy === :random || n <= 2
        @inbounds for j in 1:nag
            i = rand(rng, 1:n)
            @simd for k in 1:dim
                E[k, j] = P[k, i]
            end
        end
    else
        cd = A.cd
        @inbounds for j in 1:nag
            a = rand(rng, 1:n)
            b = rand(rng, 1:n)
            i = cd[a] >= cd[b] ? a : b
            @simd for k in 1:dim
                E[k, j] = P[k, i]
            end
        end
    end
    return E
end

# --- public (legacy-signature) wrappers ------------------------------------

"""
    update_archive(archive_prey, archive_objectives, population_prey,
                   population_objectives, max_archive_size, dim, num_obj;
                   mode = :pareto)

Public wrapper keeping the historical `nsolutions × ncols` orientation.
See [`MOArchive`](@ref) for the meaning of `mode`.
"""
function update_archive(archive_prey, archive_objectives,
    population_prey, population_objectives,
    max_archive_size_param, dim, num_obj; mode::Symbol=:pareto)

    npop = size(population_prey, 1)
    A = MOArchive(dim, num_obj, max_archive_size_param, npop; mode=mode)

    na = size(archive_prey, 1)
    if na > 0
        ensure_capacity!(A, na + npop)
        @inbounds for j in 1:na
            for k in 1:dim
                A.prey[k, j] = archive_prey[j, k]
            end
            for k in 1:num_obj
                A.obj[k, j] = archive_objectives[j, k]
            end
        end
        A.n = na
    end

    Xp = Matrix{Float64}(undef, dim, npop)
    Fp = Matrix{Float64}(undef, num_obj, npop)
    transpose_into!(Xp, population_prey)
    transpose_into!(Fp, population_objectives)

    update!(A, Xp, Fp)
    return archive_matrices(A)
end

"Return the archive in the public `n × dim` / `n × nobj` orientation."
function archive_matrices(A::MOArchive)
    P = Matrix{Float64}(undef, A.n, A.dim)
    O = Matrix{Float64}(undef, A.n, A.M)
    transpose_into!(P, view(A.prey, :, 1:A.n))
    transpose_into!(O, view(A.obj, :, 1:A.n))
    return P, O
end

"""
    select_elite_from_archive(archive_prey, archive_objectives, num_elites_needed, dim;
                              strategy = :random, rng = Random.default_rng())

Public wrapper returning a `num_elites_needed × dim` matrix of leaders.
"""
function select_elite_from_archive(archive_prey, archive_objectives, num_elites_needed, dim;
    strategy::Symbol=:random, rng::AbstractRNG=Random.default_rng())
    num_archived = size(archive_prey, 1)
    elites = zeros(Float64, num_elites_needed, dim)
    num_archived == 0 && return elites

    if strategy === :random
        @inbounds for i in 1:num_elites_needed
            idx = rand(rng, 1:num_archived)
            for k in 1:dim
                elites[i, k] = archive_prey[idx, k]
            end
        end
    else
        M = size(archive_objectives, 2)
        F = Matrix{Float64}(undef, M, num_archived)
        transpose_into!(F, archive_objectives)
        cd = zeros(Float64, num_archived)
        crowding!(cd, F, collect(1:num_archived), num_archived, M, Int[])
        @inbounds for i in 1:num_elites_needed
            a = rand(rng, 1:num_archived)
            b = rand(rng, 1:num_archived)
            idx = cd[a] >= cd[b] ? a : b
            for k in 1:dim
                elites[i, k] = archive_prey[idx, k]
            end
        end
    end
    return elites
end
