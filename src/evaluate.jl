# =============================================================================
#  evaluate.jl -- objective-function evaluation back-ends
#
#  Supported execution modes (`parallelism` keyword of MOMPA):
#
#    :serial        one process, one thread
#    :threads       Threads.@spawn over balanced chunks of the population
#    :mpi           one agent block per MPI rank, single Allgatherv per sweep
#    :mpi_threads   hybrid -- MPI across nodes, threads inside each rank
#
#  Every mode reuses preallocated buffers, so a sweep over the population
#  allocates **nothing** beyond whatever the user's `fobj` itself allocates.
# =============================================================================

"""
    EvalCtx

Fully-typed evaluation context.  Being parametric on the objective function and
on the communicator type makes the whole inner loop of MOMPA type
stable (the reference implementation stored `comm = nothing` and reassigned it,
which made every downstream call dynamically dispatched).
"""
struct EvalCtx{F,B,C}
    fobj::F
    fobj_batch::B          # `nothing` or user supplied vectorised objective
    mode::Symbol
    comm::C
    rank::Int
    nprocs::Int
    lo::Int                # first agent owned by this rank
    hi::Int                # last  agent owned by this rank
    counts::Vector{Cint}   # agents per rank
    displs::Vector{Cint}
    chunks::Vector{UnitRange{Int}}   # thread chunks of lo:hi
    bufs::Vector{Vector{Float64}}    # one x-buffer per chunk
    fbuf::Vector{Float64}            # local fitness  block  (MPI)
    obuf::Matrix{Float64}            # local objective block (MPI, nobj × nlocal)
    Xt::Matrix{Float64}              # nlocal × dim scratch for batch mode
    reuse::Bool
    root::Bool
end

is_mpi(ctx::EvalCtx) = ctx.mode === :mpi || ctx.mode === :mpi_threads
uses_threads(ctx::EvalCtx) = ctx.mode === :threads || ctx.mode === :mpi_threads

function build_evalctx(fobj, fobj_batch, mode::Symbol, comm, rank::Int, nprocs::Int,
    n::Int, dim::Int, nobj::Int, nthreads::Int, reuse::Bool,
    chunks_per_thread::Int)

    counts, displs = block_counts(n, nprocs)
    lo = Int(displs[rank+1]) + 1
    hi = lo + Int(counts[rank+1]) - 1
    nlocal = hi - lo + 1

    nchunks = (mode === :threads || mode === :mpi_threads) ?
              max(1, min(nlocal, nthreads * max(1, chunks_per_thread))) : 1
    chunks = nchunks == 1 ? [lo:hi] : [(lo-1) .+ r for r in chunk_ranges(nlocal, nchunks)]

    bufs = [Vector{Float64}(undef, dim) for _ in 1:length(chunks)]
    fbuf = Vector{Float64}(undef, nlocal)
    obuf = Matrix{Float64}(undef, nobj, nlocal)
    Xt = fobj_batch === nothing ? Matrix{Float64}(undef, 0, 0) :
         Matrix{Float64}(undef, nlocal, dim)

    return EvalCtx(fobj, fobj_batch, mode, comm, rank, nprocs, lo, hi,
        counts, displs, chunks, bufs, fbuf, obuf, Xt, reuse, rank == 0)
end

# --- position broadcast ----------------------------------------------------

"""
    sync_positions!(X, ctx)

Make rank 0's population visible to every rank.  Uses the in-place
`MPI.Bcast!` on the raw `Float64` buffer: no serialisation, no allocation --
unlike `MPI.bcast(X, 0, comm)` which round-trips the matrix through
`Serialization` on every call (the original code did this up to six times per
iteration).
"""
@inline function sync_positions!(X::Matrix{Float64}, ctx::EvalCtx)
    if is_mpi(ctx) && ctx.nprocs > 1
        MPI.Bcast!(X, ctx.comm; root=0)
    end
    return X
end

# --- single objective ------------------------------------------------------

@inline function _eval_so_range!(dest::Vector{Float64}, offset::Int, X::Matrix{Float64},
    rows::UnitRange{Int}, buf::Vector{Float64}, fobj::F, reuse::Bool) where {F}
    dim = size(X, 1)
    @inbounds for i in rows
        if reuse
            @simd for k in 1:dim
                buf[k] = X[k, i]
            end
            dest[i-offset] = sanitize(fobj(buf))
        else
            dest[i-offset] = sanitize(fobj(X[:, i]))
        end
    end
    return dest
end

function _eval_so_batch!(dest::Vector{Float64}, offset::Int, X::Matrix{Float64},
    rows::UnitRange{Int}, ctx::EvalCtx)
    dim = size(X, 1)
    Xt = ctx.Xt
    @inbounds for (r, i) in enumerate(rows), k in 1:dim
        Xt[r, k] = X[k, i]
    end
    vals = ctx.fobj_batch(Xt)
    length(vals) == length(rows) ||
        error("fobj_batch returned $(length(vals)) values for $(length(rows)) agents")
    @inbounds for (r, i) in enumerate(rows)
        dest[i-offset] = sanitize(vals[r])
    end
    return dest
end

"""
    evaluate_so!(fitness, X, ctx)

Evaluate every agent (column) of `X` into `fitness`.  On MPI the result is
consistent on all ranks after a single `Allgatherv!`.
"""
function evaluate_so!(fitness::Vector{Float64}, X::Matrix{Float64}, ctx::EvalCtx)
    sync_positions!(X, ctx)

    if is_mpi(ctx) && ctx.nprocs > 1
        rows = ctx.lo:ctx.hi
        if ctx.fobj_batch !== nothing
            _eval_so_batch!(ctx.fbuf, ctx.lo - 1, X, rows, ctx)
        elseif uses_threads(ctx)
            _eval_so_threaded!(ctx.fbuf, ctx.lo - 1, X, ctx)
        else
            _eval_so_range!(ctx.fbuf, ctx.lo - 1, X, rows, ctx.bufs[1], ctx.fobj, ctx.reuse)
        end
        MPI.Allgatherv!(ctx.fbuf, MPI.VBuffer(fitness, ctx.counts, ctx.displs), ctx.comm)
    elseif ctx.fobj_batch !== nothing
        _eval_so_batch!(fitness, 0, X, ctx.lo:ctx.hi, ctx)
    elseif uses_threads(ctx)
        _eval_so_threaded!(fitness, 0, X, ctx)
    else
        _eval_so_range!(fitness, 0, X, ctx.lo:ctx.hi, ctx.bufs[1], ctx.fobj, ctx.reuse)
    end
    return fitness
end

function _eval_so_threaded!(dest::Vector{Float64}, offset::Int, X::Matrix{Float64}, ctx::EvalCtx)
    chunks = ctx.chunks
    bufs = ctx.bufs
    fobj = ctx.fobj
    reuse = ctx.reuse
    @sync for t in eachindex(chunks)
        Threads.@spawn _eval_so_range!(dest, offset, X, chunks[t], bufs[t], fobj, reuse)
    end
    return dest
end

# --- multi objective -------------------------------------------------------

@inline function _eval_mo_range!(dest::Matrix{Float64}, offset::Int, X::Matrix{Float64},
    rows::UnitRange{Int}, buf::Vector{Float64}, fobj::F, reuse::Bool) where {F}
    dim = size(X, 1)
    nobj = size(dest, 1)
    @inbounds for i in rows
        local out
        if reuse
            @simd for k in 1:dim
                buf[k] = X[k, i]
            end
            out = fobj(buf)
        else
            out = fobj(X[:, i])
        end
        length(out) == nobj ||
            error("fobj returned $(length(out)) objectives, expected $nobj " *
                  "(check the `num_objectives` keyword)")
        col = i - offset
        @simd for k in 1:nobj
            dest[k, col] = sanitize(out[k])
        end
    end
    return dest
end

function _eval_mo_batch!(dest::Matrix{Float64}, offset::Int, X::Matrix{Float64},
    rows::UnitRange{Int}, ctx::EvalCtx)
    dim = size(X, 1)
    nobj = size(dest, 1)
    Xt = ctx.Xt
    @inbounds for (r, i) in enumerate(rows), k in 1:dim
        Xt[r, k] = X[k, i]
    end
    vals = ctx.fobj_batch(Xt)
    size(vals, 1) == length(rows) && size(vals, 2) == nobj ||
        error("fobj_batch must return a $(length(rows)) × $nobj matrix, got $(size(vals))")
    @inbounds for (r, i) in enumerate(rows), k in 1:nobj
        dest[k, i-offset] = sanitize(vals[r, k])
    end
    return dest
end

function _eval_mo_threaded!(dest::Matrix{Float64}, offset::Int, X::Matrix{Float64}, ctx::EvalCtx)
    chunks = ctx.chunks
    bufs = ctx.bufs
    fobj = ctx.fobj
    reuse = ctx.reuse
    @sync for t in eachindex(chunks)
        Threads.@spawn _eval_mo_range!(dest, offset, X, chunks[t], bufs[t], fobj, reuse)
    end
    return dest
end

"""
    evaluate_mo!(F, X, ctx)

Evaluate the population into the `nobj × nagents` matrix `F`.  Because agents
are columns, each rank's objective block is contiguous and can be gathered with
one `Allgatherv!` instead of `nagents ÷ nprocs` scatter/gather round trips.
"""
function evaluate_mo!(F::Matrix{Float64}, X::Matrix{Float64}, ctx::EvalCtx)
    sync_positions!(X, ctx)
    nobj = size(F, 1)

    if is_mpi(ctx) && ctx.nprocs > 1
        rows = ctx.lo:ctx.hi
        if ctx.fobj_batch !== nothing
            _eval_mo_batch!(ctx.obuf, ctx.lo - 1, X, rows, ctx)
        elseif uses_threads(ctx)
            _eval_mo_threaded!(ctx.obuf, ctx.lo - 1, X, ctx)
        else
            _eval_mo_range!(ctx.obuf, ctx.lo - 1, X, rows, ctx.bufs[1], ctx.fobj, ctx.reuse)
        end
        counts = ctx.counts .* Cint(nobj)
        displs = ctx.displs .* Cint(nobj)
        MPI.Allgatherv!(ctx.obuf, MPI.VBuffer(F, counts, displs), ctx.comm)
    elseif ctx.fobj_batch !== nothing
        _eval_mo_batch!(F, 0, X, ctx.lo:ctx.hi, ctx)
    elseif uses_threads(ctx)
        _eval_mo_threaded!(F, 0, X, ctx)
    else
        _eval_mo_range!(F, 0, X, ctx.lo:ctx.hi, ctx.bufs[1], ctx.fobj, ctx.reuse)
    end
    return F
end

# --- tiny collective helpers ----------------------------------------------

@inline function bcast_scalar(x::T, ctx::EvalCtx) where {T}
    if is_mpi(ctx) && ctx.nprocs > 1
        r = Ref(x)
        MPI.Bcast!(r, ctx.comm; root=0)
        return r[]
    end
    return x
end

@inline function bcast_vec!(v::Vector{Float64}, ctx::EvalCtx)
    if is_mpi(ctx) && ctx.nprocs > 1
        MPI.Bcast!(v, ctx.comm; root=0)
    end
    return v
end
