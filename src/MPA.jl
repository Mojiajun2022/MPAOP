# =============================================================================
#  MPA.jl -- Marine Predators Algorithm drivers
#
#    SOMPA : single objective
#    MOMPA : multi objective (external archive + Pareto ranking)
#
#  Both keep the keyword interface of MPAOP ≤ 0.2; every new keyword has a
#  default that reproduces the old behaviour unless documented otherwise.
# =============================================================================

"Immutable, fully typed configuration handed to the (type stable) inner loops."
struct MPARun
    dim::Int
    n::Int
    Max_iter::Int
    t1::Int
    t2::Int
    nmpa::Bool
    FADs::Float64
    P::Float64
    Fixbox::Bool
    disp::Bool
    disp_param::Bool
    disp_every::Int
    write_csv_log::Bool
    csv_log_filepath::String
    csv_flush_every::Int
    saveHDF::Bool
    hdf_filepath::String
    history_save_interval::Int
    hdf_compress::Int
    ftol::Float64
    patience::Int
    max_time::Float64
    sigma::Float64
    inv_beta::Float64
end

@inline function _cf(nmpa::Bool, Iter::Int, Max_iter::Int)
    t = Iter / Max_iter
    return nmpa ? abs(2 * (1 - t) - 2) : (1 - t)^(2 * t)
end

@inline function _wfac(nmpa::Bool, Iter::Int, Max_iter::Int)
    return nmpa ? 2 * exp(-(6 * (Iter / Max_iter))^2) : 1.0
end

@inline _stage(Iter::Int, t1::Int, t2::Int) = Iter < t1 ? 1 : (Iter < t2 ? 2 : 3)

function _setup_mpi(parallelism::Symbol)
    if parallelism === :mpi || parallelism === :mpi_threads
        MPI.Initialized() || MPI.Init()
        comm = MPI.COMM_WORLD
        return comm, MPI.Comm_rank(comm), MPI.Comm_size(comm)
    end
    return nothing, 0, 1
end

function _flush_csv!(buf::Vector{NamedTuple{(:text,),Tuple{String}}}, path::String)
    isempty(buf) && return
    try
        CSV.write(path, buf; append=true, header=false)
    catch e
        @warn "Could not append to CSV log $path" exception = e
    end
    empty!(buf)
    return
end

# =============================================================================
#  Single objective
# =============================================================================

"""
    SOMPA(; fobj, lb, ub, SearchAgents_no, Max_iter, kwargs...)
        -> (best_fitness, best_position, convergence_curve)

Single-objective Marine Predators Algorithm.

# Required
* `fobj`             -- `f(x::Vector{Float64}) -> Real`
* `lb`, `ub`         -- box bounds, `Vector{Float64}` of equal length
* `SearchAgents_no`  -- population size
* `Max_iter`         -- iterations (each one costs `2 · SearchAgents_no` evaluations)

# Classic options (unchanged)
`p0_optional`, `variant` (`:standard_mpa` | `:nmpa`), `first_stage_ratio`,
`second_stage_ratio`, `FADs0`, `P0`, `parallelism`, `disp`, `disp_param`,
`Fixbox`, `write_csv_log`, `csv_log_filepath`, `saveHDF`, `hdf_filepath`,
`history_save_interval`.

# Added in 0.3
* `parallelism = :mpi_threads` -- hybrid MPI + threads
* `fobj_batch`   -- `f(X::Matrix) -> Vector`, `X` is `nagents × dim`; lets you
  vectorise / GPU-accelerate the whole population in a single call
* `seed`, `rng`  -- bit-reproducible runs (also under MPI)
* `init`         -- `:uniform` (default), `:lhs`, `:center`
* `opposition`   -- opposition-based learning at start-up (one extra generation)
* `ftol`, `patience`, `max_time` -- early stopping
* `callback(iter, best_fit, best_pos, curve)` -- return `false`/`:stop` to abort
* `disp_every`, `csv_flush_every`, `hdf_compress`
* `nthreads`, `chunks_per_thread`, `reuse_buffer`
* `round_digits` -- `nothing` (default) returns full precision.  MPAOP ≤ 0.2
  silently rounded the returned fitness to 4 digits and the position to 8;
  pass `round_digits = 4` to restore that.

Every rank returns the same triple, so no `if rank == 0` guard is needed.
"""
function SOMPA(;
    fobj::Function,
    lb::Vector{Float64},
    ub::Vector{Float64},
    SearchAgents_no::Int64,
    Max_iter::Int64,
    p0_optional::Vector=[],
    variant::Symbol=:standard_mpa,
    first_stage_ratio::Float64=1 / 3,
    second_stage_ratio::Float64=2 / 3,
    FADs0::Float64=0.2,
    P0::Float64=0.5,
    parallelism::Symbol=:serial,
    disp::Bool=true,
    disp_param::Bool=false,
    Fixbox::Bool=true,
    write_csv_log::Bool=false,
    csv_log_filepath::String="mpa_so_log.csv",
    saveHDF::Bool=false,
    hdf_filepath::String="mpa_so_history.h5",
    history_save_interval::Int=typemax(Int),
    # ---------------- new ----------------
    fobj_batch=nothing,
    seed::Union{Nothing,Integer}=nothing,
    rng::Union{Nothing,AbstractRNG}=nothing,
    init::Symbol=:uniform,
    opposition::Bool=false,
    beta::Float64=1.5,
    ftol::Float64=0.0,
    patience::Int=typemax(Int),
    max_time::Float64=Inf,
    callback=nothing,
    disp_every::Int=1,
    csv_flush_every::Int=1,
    hdf_compress::Int=0,
    nthreads::Int=Threads.nthreads(),
    chunks_per_thread::Int=2,
    reuse_buffer::Bool=true,
    round_digits::Union{Nothing,Int}=nothing,
)
    dim = length(lb)
    length(ub) == dim || throw(DimensionMismatch("lb and ub must have the same length"))
    SearchAgents_no > 0 || throw(ArgumentError("SearchAgents_no must be positive"))
    Max_iter > 0 || throw(ArgumentError("Max_iter must be positive"))
    variant in (:standard_mpa, :nmpa) ||
        throw(ArgumentError("variant must be :standard_mpa or :nmpa"))
    parallelism in (:serial, :threads, :mpi, :mpi_threads) ||
        throw(ArgumentError("parallelism must be :serial, :threads, :mpi or :mpi_threads"))
    @inbounds for i in 1:dim
        lb[i] <= ub[i] || throw(ArgumentError("lb[$i] > ub[$i]"))
    end

    comm, rank, nprocs = _setup_mpi(parallelism)
    n = SearchAgents_no
    ctx = build_evalctx(fobj, fobj_batch, parallelism, comm, rank, nprocs,
        n, dim, 1, nthreads, reuse_buffer, chunks_per_thread)

    if parallelism === :threads || parallelism === :mpi_threads
        rank == 0 && Threads.nthreads() == 1 &&
            @warn "parallelism = :$parallelism but Julia runs with 1 thread; start it with `julia -t auto`."
    end

    cfg = MPARun(dim, n, Max_iter,
        Int(round(Max_iter * first_stage_ratio)),
        Int(round(Max_iter * second_stage_ratio)),
        variant === :nmpa, FADs0, P0, Fixbox,
        disp, disp_param, max(1, disp_every),
        write_csv_log, csv_log_filepath, max(1, csv_flush_every),
        saveHDF, hdf_filepath, history_save_interval, hdf_compress,
        ftol, patience, max_time, levy_sigma(beta), 1 / beta)

    return _sompa_core(ctx, cfg, callback, make_rng(rng, seed), lb, ub,
        p0_optional, init, opposition, variant, parallelism,
        round_digits === nothing ? -1 : round_digits)
end

function _sompa_core(ctx::EvalCtx, cfg::MPARun, cb::CB, rng::AbstractRNG,
    lb::Vector{Float64}, ub::Vector{Float64}, p0,
    init_method::Symbol, opposition::Bool,
    variant::Symbol, parallelism::Symbol, round_digits::Int) where {CB}

    dim, n, Max_iter = cfg.dim, cfg.n, cfg.Max_iter
    root = ctx.root

    X = Matrix{Float64}(undef, dim, n)
    Xn = Matrix{Float64}(undef, dim, n)
    X_old = Matrix{Float64}(undef, dim, n)
    fitness = fill(Inf, n)
    fit_old = fill(Inf, n)
    best_pos = zeros(Float64, dim)
    best_fit = Inf
    curve = zeros(Float64, Max_iter)
    perm1 = collect(1:n)
    perm2 = collect(1:n)

    nsnap = cfg.saveHDF && root ? 2 * Max_iter : 0
    hist_X = Array{Float64,3}(undef, n, dim, nsnap)
    hist_f = Matrix{Float64}(undef, n, nsnap)

    csvbuf = NamedTuple{(:text,),Tuple{String}}[]

    if root
        init_population!(X, lb, ub, p0, rng; method=init_method)
        if cfg.write_csv_log
            push!(csvbuf, (text=string(Dates.now(),
                    " SingleObjectiveMPA ($variant, $parallelism) start running!\n",
                    "______________________________"),))
            _flush_csv!(csvbuf, cfg.csv_log_filepath)
        end
    end
    sync_positions!(X, ctx)

    if opposition
        evaluate_so!(fitness, X, ctx)
        opposite_population!(Xn, X, lb, ub)
        cfg.Fixbox && clamp_cols!(Xn, lb, ub)
        evaluate_so!(fit_old, Xn, ctx)
        @inbounds for j in 1:n
            if fit_old[j] < fitness[j]
                fitness[j] = fit_old[j]
                @simd for i in 1:dim
                    X[i, j] = Xn[i, j]
                end
            end
        end
        fill!(fit_old, Inf)
    end

    Iter = 0
    total_time = 0.0
    stagnation = 0
    prev_best = Inf
    stop = false
    t_run0 = time_ns()

    while Iter < Max_iter && !stop
        iter_t0 = time_ns()
        CF = _cf(cfg.nmpa, Iter, Max_iter)

        # ---------------- sweep A : evaluate current population ------------
        cfg.Fixbox && clamp_cols!(X, lb, ub)
        evaluate_so!(fitness, X, ctx)

        if root && cfg.saveHDF
            transpose_into!(view(hist_X, :, :, 2Iter + 1), X)
            @inbounds copyto!(view(hist_f, :, 2Iter + 1), fitness)
        end

        mn, ind = findmin(fitness)
        if mn < best_fit
            best_fit = mn
            @inbounds for i in 1:dim
                best_pos[i] = X[i, ind]
            end
        end
        if Iter == 0
            copyto!(X_old, X)
            copyto!(fit_old, fitness)
        end
        _memory_step!(X, fitness, X_old, fit_old, dim, n)

        # ---------------- movement (rank 0 only, broadcast on next eval) ---
        if root
            mpa_move!(Xn, X, best_pos, _stage(Iter, cfg.t1, cfg.t2),
                cfg.P, CF, _wfac(cfg.nmpa, Iter, Max_iter), rng, cfg.sigma, cfg.inv_beta)
            X, Xn = Xn, X
            cfg.Fixbox && clamp_cols!(X, lb, ub)
        end

        # ---------------- sweep B : evaluate moved population --------------
        evaluate_so!(fitness, X, ctx)

        if root && cfg.saveHDF
            transpose_into!(view(hist_X, :, :, 2Iter + 2), X)
            @inbounds copyto!(view(hist_f, :, 2Iter + 2), fitness)
        end

        mn, ind = findmin(fitness)
        if mn < best_fit
            best_fit = mn
            @inbounds for i in 1:dim
                best_pos[i] = X[i, ind]
            end
        end
        _memory_step!(X, fitness, X_old, fit_old, dim, n)

        # ---------------- FADs / eddy formation ----------------------------
        if root
            fads!(Xn, X, lb, ub, cfg.FADs, CF, rng, perm1, perm2)
            X, Xn = Xn, X
            cfg.Fixbox && clamp_cols!(X, lb, ub)
        end

        # ---------------- bookkeeping --------------------------------------
        Iter += 1
        curve[Iter] = best_fit
        iter_duration = (time_ns() - iter_t0) / 1e9
        total_time += iter_duration

        if root
            if cfg.disp && (Iter % cfg.disp_every == 0 || Iter == Max_iter)
                print_buffer = "SO-MPA Iter: $Iter/$Max_iter | BestFit: $(round(best_fit, digits=4)) | IterTime: $(round(iter_duration, digits=3))s | TotalTime: $(round(total_time, digits=2))s"
                if cfg.disp_param
                    print_buffer *= "\n  BestParams SO: $(round.(best_pos, digits=3))"
                end
                println(print_buffer)
            end
            if cfg.write_csv_log
                push!(csvbuf, (text=string("Iter: $Iter, BestFit: $(round(best_fit, digits=4)), BestParams: $(round.(best_pos, digits=3))"),))
                (Iter % cfg.csv_flush_every == 0) && _flush_csv!(csvbuf, cfg.csv_log_filepath)
            end
            if cfg.saveHDF && (Iter % cfg.history_save_interval == 0 || Iter == Max_iter)
                SaveMPAHistory(cfg.hdf_filepath,
                    view(hist_X, :, :, 1:2Iter), view(hist_f, :, 1:2Iter), false;
                    convergence_curve=view(curve, 1:Iter),
                    compress=cfg.hdf_compress, is_root_process=true)
            end

            if prev_best - best_fit > cfg.ftol
                stagnation = 0
                prev_best = best_fit
            else
                stagnation += 1
            end
            stop = stagnation >= cfg.patience ||
                   (time_ns() - t_run0) / 1e9 >= cfg.max_time
            if cb !== nothing
                r = cb(Iter, best_fit, best_pos, view(curve, 1:Iter))
                (r === false || r === :stop) && (stop = true)
            end
            stop && cfg.disp && println("SO-MPA: early stop at iteration $Iter.")
        end
        stop = bcast_scalar(Int(stop), ctx) != 0
    end

    # keep the historical curve length; pad with the final best
    @inbounds for k in (Iter+1):Max_iter
        curve[k] = best_fit
    end

    if root
        if cfg.write_csv_log
            push!(csvbuf, (text=string(Dates.now(),
                    " SingleObjectiveMPA ($variant, $parallelism) finish!\n",
                    "______________________________"),))
            _flush_csv!(csvbuf, cfg.csv_log_filepath)
        end
        if cfg.saveHDF && Iter > 0 && Iter % cfg.history_save_interval != 0 && Iter != Max_iter
            SaveMPAHistory(cfg.hdf_filepath,
                view(hist_X, :, :, 1:2Iter), view(hist_f, :, 1:2Iter), false;
                convergence_curve=view(curve, 1:Iter),
                compress=cfg.hdf_compress, is_root_process=true)
        end
        cfg.disp && println("SingleObjectiveMPA ($variant, $parallelism) finished. Best Fitness: $best_fit")
    end

    if round_digits >= 0
        return round(best_fit, digits=round_digits), round.(best_pos, digits=round_digits), curve
    end
    return best_fit, best_pos, curve
end

"Greedy memory: any agent that got worse is restored to its previous state."
@inline function _memory_step!(X::Matrix{Float64}, fitness::Vector{Float64},
    X_old::Matrix{Float64}, fit_old::Vector{Float64}, dim::Int, n::Int)
    @inbounds for j in 1:n
        if fit_old[j] < fitness[j]
            fitness[j] = fit_old[j]
            @simd for i in 1:dim
                X[i, j] = X_old[i, j]
            end
        else
            fit_old[j] = fitness[j]
            @simd for i in 1:dim
                X_old[i, j] = X[i, j]
            end
        end
    end
    return
end

# =============================================================================
#  Multi objective
# =============================================================================

"""
    MOMPA(; fobj, lb, ub, SearchAgents_no, Max_iter, num_objectives, kwargs...)
        -> (archive_positions, archive_objectives, convergence_curve)

Multi-objective Marine Predators Algorithm with an external archive.
`archive_positions` is `narchive × dim`, `archive_objectives` is
`narchive × num_objectives`, `convergence_curve[k]` is the archive size after
iteration `k`.  Under MPI, non-root ranks return `(nothing, nothing, nothing)`
(unchanged).

Shares every keyword of [`SOMPA`](@ref) plus:

* `archive_size_factor`  -- archive capacity = `factor · SearchAgents_no`
* `archive_mode`         -- `:pareto` (default, archive == non-dominated set) or
  `:fronts` (MPAOP ≤ 0.2 behaviour, keeps dominated solutions once the first
  front is smaller than the archive)
* `elite_selection`      -- `:crowding` (default, binary tournament on crowding
  distance) or `:random` (MPAOP ≤ 0.2 behaviour)
* `hv_ref`               -- reference point; when given, the hypervolume is
  reported alongside the archive size

!!! note "Behaviour change"
    With the defaults the returned archive is guaranteed to contain **no
    dominated solution**, and leaders are drawn from sparse regions of the
    front.  Pass `archive_mode = :fronts, elite_selection = :random` to
    reproduce results from MPAOP ≤ 0.2 exactly.
"""
function MOMPA(;
    fobj::Function,
    lb::Vector{Float64},
    ub::Vector{Float64},
    SearchAgents_no::Int64,
    Max_iter::Int64,
    num_objectives::Int,
    p0_optional::Vector=[],
    variant::Symbol=:standard_mpa,
    first_stage_ratio::Float64=1 / 3,
    second_stage_ratio::Float64=2 / 3,
    FADs0::Float64=0.2,
    P0::Float64=0.5,
    archive_size_factor=1.0,
    parallelism::Symbol=:serial,
    disp::Bool=true,
    disp_param::Bool=false,
    Fixbox::Bool=true,
    write_csv_log::Bool=false,
    csv_log_filepath::String="mompa_log.csv",
    saveHDF::Bool=false,
    hdf_filepath::String="mompa_history.h5",
    history_save_interval::Int=typemax(Int),
    # ---------------- new ----------------
    archive_mode::Symbol=:pareto,
    elite_selection::Symbol=:crowding,
    hv_ref::Union{Nothing,Vector{Float64}}=nothing,
    fobj_batch=nothing,
    seed::Union{Nothing,Integer}=nothing,
    rng::Union{Nothing,AbstractRNG}=nothing,
    init::Symbol=:uniform,
    beta::Float64=1.5,
    max_time::Float64=Inf,
    callback=nothing,
    disp_every::Int=1,
    csv_flush_every::Int=1,
    hdf_compress::Int=0,
    nthreads::Int=Threads.nthreads(),
    chunks_per_thread::Int=2,
    reuse_buffer::Bool=true,
)
    dim = length(lb)
    length(ub) == dim || throw(DimensionMismatch("lb and ub must have the same length"))
    SearchAgents_no > 0 || throw(ArgumentError("SearchAgents_no must be positive"))
    Max_iter > 0 || throw(ArgumentError("Max_iter must be positive"))
    num_objectives >= 2 ||
        throw(ArgumentError("num_objectives must be ≥ 2 (use SOMPA for one objective)"))
    variant in (:standard_mpa, :nmpa) ||
        throw(ArgumentError("variant must be :standard_mpa or :nmpa"))
    archive_mode in (:pareto, :fronts) ||
        throw(ArgumentError("archive_mode must be :pareto or :fronts"))
    elite_selection in (:crowding, :random) ||
        throw(ArgumentError("elite_selection must be :crowding or :random"))
    parallelism in (:serial, :threads, :mpi, :mpi_threads) ||
        throw(ArgumentError("parallelism must be :serial, :threads, :mpi or :mpi_threads"))
    @inbounds for i in 1:dim
        lb[i] <= ub[i] || throw(ArgumentError("lb[$i] > ub[$i]"))
    end

    comm, rank, nprocs = _setup_mpi(parallelism)
    n = SearchAgents_no
    ctx = build_evalctx(fobj, fobj_batch, parallelism, comm, rank, nprocs,
        n, dim, num_objectives, nthreads, reuse_buffer, chunks_per_thread)

    cfg = MPARun(dim, n, Max_iter,
        Int(round(Max_iter * first_stage_ratio)),
        Int(round(Max_iter * second_stage_ratio)),
        variant === :nmpa, FADs0, P0, Fixbox,
        disp, disp_param, max(1, disp_every),
        write_csv_log, csv_log_filepath, max(1, csv_flush_every),
        saveHDF, hdf_filepath, history_save_interval, hdf_compress,
        0.0, typemax(Int), max_time, levy_sigma(beta), 1 / beta)

    max_archive = max(1, Int(round(archive_size_factor * SearchAgents_no)))

    return _mompa_core(ctx, cfg, callback, make_rng(rng, seed), lb, ub, p0_optional,
        init, num_objectives, max_archive, archive_mode, elite_selection,
        hv_ref, parallelism)
end

function _mompa_core(ctx::EvalCtx, cfg::MPARun, cb::CB, rng::AbstractRNG,
    lb::Vector{Float64}, ub::Vector{Float64}, p0, init_method::Symbol,
    M::Int, max_archive::Int, archive_mode::Symbol, elite_selection::Symbol,
    hv_ref::Union{Nothing,Vector{Float64}}, parallelism::Symbol) where {CB}

    dim, n, Max_iter = cfg.dim, cfg.n, cfg.Max_iter
    root = ctx.root

    X = Matrix{Float64}(undef, dim, n)
    Xn = Matrix{Float64}(undef, dim, n)
    E = Matrix{Float64}(undef, dim, n)
    Fo = Matrix{Float64}(undef, M, n)
    curve = zeros(Float64, Max_iter)
    perm1 = collect(1:n)
    perm2 = collect(1:n)
    A = MOArchive(dim, M, max_archive, n; mode=archive_mode)

    nsnap = cfg.saveHDF && root ? 2 * Max_iter : 0
    hist_X = Array{Float64,3}(undef, n, dim, nsnap)
    hist_O = Array{Float64,3}(undef, n, M, nsnap)

    csvbuf = NamedTuple{(:text,),Tuple{String}}[]

    if root
        init_population!(X, lb, ub, p0, rng; method=init_method)
        if cfg.write_csv_log
            push!(csvbuf, (text=string(Dates.now(),
                    " MultiObjectiveMPA ($parallelism) start running!\n",
                    "______________________________"),))
            _flush_csv!(csvbuf, cfg.csv_log_filepath)
        end
    end
    sync_positions!(X, ctx)

    Iter = 0
    total_time = 0.0
    stop = false
    t_run0 = time_ns()

    while Iter < Max_iter && !stop
        iter_t0 = time_ns()
        CF = _cf(cfg.nmpa, Iter, Max_iter)

        # ---------------- sweep A ------------------------------------------
        cfg.Fixbox && clamp_cols!(X, lb, ub)
        evaluate_mo!(Fo, X, ctx)

        if root && cfg.saveHDF
            transpose_into!(view(hist_X, :, :, 2Iter + 1), X)
            transpose_into!(view(hist_O, :, :, 2Iter + 1), Fo)
        end

        if root
            update!(A, X, Fo)
            fill_elites!(E, A, rng, elite_selection)
            mpa_move!(Xn, X, E, _stage(Iter, cfg.t1, cfg.t2),
                cfg.P, CF, _wfac(cfg.nmpa, Iter, Max_iter), rng, cfg.sigma, cfg.inv_beta)
            X, Xn = Xn, X
            cfg.Fixbox && clamp_cols!(X, lb, ub)
        end

        # ---------------- sweep B ------------------------------------------
        evaluate_mo!(Fo, X, ctx)

        if root && cfg.saveHDF
            transpose_into!(view(hist_X, :, :, 2Iter + 2), X)
            transpose_into!(view(hist_O, :, :, 2Iter + 2), Fo)
        end

        if root
            update!(A, X, Fo)
            fads!(Xn, X, lb, ub, cfg.FADs, CF, rng, perm1, perm2)
            X, Xn = Xn, X
            cfg.Fixbox && clamp_cols!(X, lb, ub)
        end

        # ---------------- bookkeeping --------------------------------------
        Iter += 1
        iter_duration = (time_ns() - iter_t0) / 1e9
        total_time += iter_duration

        if root
            curve[Iter] = A.n
            hv = NaN
            if hv_ref !== nothing && A.n > 0
                _, O = archive_matrices(A)
                hv = hypervolume(O, hv_ref)
            end

            if cfg.disp && (Iter % cfg.disp_every == 0 || Iter == Max_iter)
                print_buffer = "MO-MPA Iter: $Iter/$Max_iter | ArchiveSize: $(A.n) | IterTime: $(round(iter_duration, digits=3))s | TotalTime: $(round(total_time, digits=2))s"
                isnan(hv) || (print_buffer *= " | HV: $(round(hv, digits=6))")
                if cfg.disp_param && A.n > 0
                    print_buffer *= "\n  Sample Archive Params MO: $(round.(view(A.prey, :, 1), digits=3))"
                end
                println(print_buffer)
            end

            if cfg.write_csv_log
                if A.n > 0
                    for k in 1:M
                        best = 1
                        @inbounds for i in 2:A.n
                            A.obj[k, i] < A.obj[k, best] && (best = i)
                        end
                        push!(csvbuf, (text=string(
                                "*As of the ", Iter, "-th iteration, for Objective ", k,
                                "\n the minimum value is ", round.(view(A.obj, :, best), digits=4),
                                "\n achieved with parameters: [",
                                join(round.(view(A.prey, :, best), digits=3), ", "), "]\n"),))
                    end
                else
                    push!(csvbuf, (text=string("Iter:$Iter,ArchiveStatus:Empty"),))
                end
                (Iter % cfg.csv_flush_every == 0) && _flush_csv!(csvbuf, cfg.csv_log_filepath)
            end

            if cfg.saveHDF && (Iter % cfg.history_save_interval == 0 || Iter == Max_iter)
                AP, AO = archive_matrices(A)
                SaveMPAHistory(cfg.hdf_filepath,
                    view(hist_X, :, :, 1:2Iter), view(hist_O, :, :, 1:2Iter), true;
                    archive_prey_at_save=AP, archive_objectives_at_save=AO,
                    convergence_curve=view(curve, 1:Iter),
                    compress=cfg.hdf_compress, is_root_process=true)
            end

            stop = (time_ns() - t_run0) / 1e9 >= cfg.max_time
            if cb !== nothing
                cbP, cbO = archive_matrices(A)
                r = cb(Iter, cbP, cbO, view(curve, 1:Iter))
                (r === false || r === :stop) && (stop = true)
            end
            stop && cfg.disp && println("MO-MPA: early stop at iteration $Iter.")
        end
        stop = bcast_scalar(Int(stop), ctx) != 0
    end

    if !root
        return nothing, nothing, nothing
    end

    if cfg.write_csv_log
        push!(csvbuf, (text=string(Dates.now(),
                " MultiObjectiveMPA ($parallelism) finish!\n",
                "______________________________"),))
        _flush_csv!(csvbuf, cfg.csv_log_filepath)
    end

    AP, AO = archive_matrices(A)

    if cfg.saveHDF && Iter > 0 && Iter % cfg.history_save_interval != 0 && Iter != Max_iter
        SaveMPAHistory(cfg.hdf_filepath,
            view(hist_X, :, :, 1:2Iter), view(hist_O, :, :, 1:2Iter), true;
            archive_prey_at_save=AP, archive_objectives_at_save=AO,
            convergence_curve=view(curve, 1:Iter),
            compress=cfg.hdf_compress, is_root_process=true)
    end
    cfg.disp && println("MultiObjectiveMPA ($parallelism) finished. Final Archive Size: $(size(AP, 1))")

    return AP, AO, curve[1:Iter]
end
