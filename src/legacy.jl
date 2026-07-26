# =============================================================================
#  legacy.jl -- compatibility layer for scripts written against MPAOP 0.1
#
#  The published 0.1 API was positional (`MPA` / `MPA_MPI`).  These thin
#  wrappers forward to the current engine, so old scripts keep running while
#  getting the new performance for free.  `Convergence_curve` keeps its old
#  `1 × Max_iter` shape so that `plot(CV')` still works.
# =============================================================================

"""
    MPA(SearchAgents_no, Max_iter, p0, lb, ub, dim, fobj;
        disp = true, Fixbox = true, Threads_parallel = false,
        Write = false, FADs0 = 0.2, P0 = 0.5)
        -> (Top_predator_fit, Top_predator_pos, Convergence_curve)

Legacy positional interface of MPAOP 0.1.  Equivalent to

```julia
MOMPA(fobj = fobj, lb = lb, ub = ub, SearchAgents_no = SearchAgents_no,
      Max_iter = Max_iter, num_objectives = 1, p0_optional = p0,
      disp = disp, Fixbox = Fixbox,
      parallelism = Threads_parallel ? :threads : :serial)
```

New code should call [`MOMPA`](@ref) directly -- it exposes seeding, early
stopping, batch evaluation, hybrid MPI and multi-objective optimisation.
"""
function MPA(SearchAgents_no::Int64, Max_iter::Int64, p0, lb::Vector{Float64},
    ub::Vector{Float64}, dim::Int64, fobj::Function;
    disp::Bool=true, Fixbox::Bool=true, Threads_parallel::Bool=false,
    Write::Bool=false, FADs0=0.2, P0=0.5, kwargs...)

    dim == length(lb) ||
        @warn "MPA: `dim` ($dim) does not match `length(lb)` ($(length(lb))); using length(lb)."
    fit, pos, curve = MOMPA(fobj=fobj, lb=lb, ub=ub, num_objectives=1,
        SearchAgents_no=SearchAgents_no, Max_iter=Max_iter,
        p0_optional=collect(p0), disp=disp, Fixbox=Fixbox,
        FADs0=Float64(FADs0), P0=Float64(P0),
        parallelism=Threads_parallel ? :threads : :serial,
        write_csv_log=Write, csv_log_filepath="fitting_process",
        kwargs...)
    return fit, pos, reshape(curve, 1, :)
end

"""
    MPA_MPI(SearchAgents_no, Max_iter, p0, lb, ub, dim, fobj;
            disp = true, Fixbox = true, Write = false, FADs0 = 0.2, P0 = 0.5)

Legacy MPI entry point; identical to `MPA` with `parallelism = :mpi`.
Run it with `mpiexec -n <ranks> julia script.jl`.
"""
function MPA_MPI(SearchAgents_no::Int64, Max_iter::Int64, p0, lb::Vector{Float64},
    ub::Vector{Float64}, dim::Int64, fobj::Function;
    disp::Bool=true, Fixbox::Bool=true, Write::Bool=false,
    FADs0=0.2, P0=0.5, kwargs...)

    fit, pos, curve = MOMPA(fobj=fobj, lb=lb, ub=ub, num_objectives=1,
        SearchAgents_no=SearchAgents_no, Max_iter=Max_iter,
        p0_optional=collect(p0), disp=disp, Fixbox=Fixbox,
        FADs0=Float64(FADs0), P0=Float64(P0), parallelism=:mpi,
        write_csv_log=Write, csv_log_filepath="fitting_process",
        kwargs...)
    return fit, pos, reshape(curve, 1, :)
end

"""
    initialization(SearchAgents_no, dim, ub, lb) -> Matrix

Four-argument form of [`initialization`](@ref).

!!! warning
    MPAOP 0.1 returned this matrix transposed (`dim × SearchAgents_no`).  Since
    0.2 both methods return `SearchAgents_no × dim`, i.e. one agent per row.
"""
initialization(SearchAgents_no::Int64, dim::Int64, ub::Vector{Float64}, lb::Vector{Float64}) =
    initialization(SearchAgents_no, dim, ub, lb, [])

# --- normal quantile without pulling in Distributions ----------------------
@inline _norm_quantile(p::Float64) = sqrt(2.0) * erfinv(2p - 1)

"""
    confidence_interval(param, func, level) -> Matrix

Asymptotic confidence interval of a least-squares / maximum-likelihood fit,
obtained from the finite-difference Hessian of `func` at `param`:

    C = H⁻¹,   σᵢ = √Cᵢᵢ,   CIᵢ = paramᵢ ± z_{(1+level)/2} · σᵢ

`H⁻¹` is projected onto the positive-semidefinite cone (eigenvalues clamped at
zero) so that a slightly indefinite numerical Hessian cannot produce a negative
variance.  Returns an `n × 2` matrix of lower/upper bounds.

Since 0.3 this needs neither `FiniteDiff`, `PositiveFactorizations` nor
`Distributions` -- it is implemented with `LinearAlgebra` + `SpecialFunctions`
only, which is why loading MPAOP is now several times faster.
"""
function confidence_interval(param::Vector{Float64}, func::Function, level::Float64)
    0 < level < 1 || throw(ArgumentError("level must be in (0, 1)"))
    n = length(param)
    H = _fd_hessian(func, param)
    C = inv(Symmetric(H))
    F = eigen(Symmetric(Matrix(C)))
    vals = max.(F.values, 0.0)
    var = [sum(F.vectors[i, k]^2 * vals[k] for k in 1:n) for i in 1:n]
    sd = sqrt.(var)
    z = _norm_quantile((1 + level) / 2)
    ci = Matrix{Float64}(undef, n, 2)
    @inbounds for i in 1:n
        ci[i, 1] = param[i] - z * sd[i]
        ci[i, 2] = param[i] + z * sd[i]
    end
    return ci
end

"Central-difference Hessian with a per-coordinate step ∝ cbrt(eps)."
function _fd_hessian(f::F, x::Vector{Float64}) where {F}
    n = length(x)
    h = [cbrt(eps(Float64)) * max(abs(xi), 1.0) for xi in x]
    H = Matrix{Float64}(undef, n, n)
    xp = copy(x)
    f0 = f(x)
    @inbounds for i in 1:n
        xp .= x
        xp[i] = x[i] + h[i]
        fp = f(xp)
        xp[i] = x[i] - h[i]
        fm = f(xp)
        H[i, i] = (fp - 2 * f0 + fm) / (h[i]^2)
        for j in (i+1):n
            xp .= x
            xp[i] = x[i] + h[i]
            xp[j] = x[j] + h[j]
            fpp = f(xp)
            xp[j] = x[j] - h[j]
            fpm = f(xp)
            xp[i] = x[i] - h[i]
            fmm = f(xp)
            xp[j] = x[j] + h[j]
            fmp = f(xp)
            H[i, j] = H[j, i] = (fpp - fpm - fmp + fmm) / (4h[i] * h[j])
        end
    end
    return H
end
