# =============================================================================
#  Tour of the features added in MPAOP 0.3.
#  Run with:   julia -t auto --project=. example/example_v03_features.jl
# =============================================================================
using MPAOP
using Random
using Printf

# ---------------------------------------------------------------- objectives
sphere(x::Vector{Float64}) = sum(abs2, x)

# shifted sphere -- the optimum is NOT at the centre of the box, so opposition
# based learning and the different initialisations actually have something to do
shifted(x::Vector{Float64}) = sum(abs2, x .- 3.7)

function zdt1(x::Vector{Float64})
    f1 = x[1]
    g = 1.0 + 9.0 / (length(x) - 1) * sum(@view x[2:end])
    return [f1, g * (1 - sqrt(f1 / g))]
end

const D = 20
lb_so, ub_so = fill(-10.0, D), fill(10.0, D)
lb_mo, ub_mo = zeros(D), ones(D)

# ---------------------------------------------------------------------------
println("\n=== 1. Reproducibility: the same seed gives the same answer ===")
a = SOMPA(fobj=sphere, lb=lb_so, ub=ub_so, SearchAgents_no=40, Max_iter=100,
    disp=false, seed=2025)
b = SOMPA(fobj=sphere, lb=lb_so, ub=ub_so, SearchAgents_no=40, Max_iter=100,
    disp=false, seed=2025)
@printf("serial  : %.6e\nrepeat  : %.6e   identical = %s\n", a[1], b[1], a[1] == b[1])

# ...and threads reproduce it exactly, because every random decision is taken
# on rank 0 / the main task and only the evaluations are distributed.
t = SOMPA(fobj=sphere, lb=lb_so, ub=ub_so, SearchAgents_no=40, Max_iter=100,
    disp=false, seed=2025, parallelism=:threads)
@printf("threads : %.6e   identical = %s   (%d threads)\n",
    t[1], t[1] == a[1], Threads.nthreads())

# ---------------------------------------------------------------------------
println("\n=== 2. Better initial populations (mean over 20 seeds, shifted sphere) ===")
for ini in (:uniform, :lhs, :center), opp in (false, true)
    vals = [SOMPA(fobj=shifted, lb=lb_so, ub=ub_so, SearchAgents_no=40, Max_iter=100,
        disp=false, seed=s, init=ini, opposition=opp)[1] for s in 1:20]
    @printf("  init=%-8s opposition=%-5s -> mean %.4e   best %.4e\n",
        ini, opp, sum(vals) / length(vals), minimum(vals))
end
# note: on a symmetric objective centred in the box (plain `sphere` with
# lb = -ub) opposition is a no-op, because f(lb+ub-x) == f(x) there.

# ---------------------------------------------------------------------------
println("\n=== 3. Early stopping ===")
r = SOMPA(fobj=sphere, lb=lb_so, ub=ub_so, SearchAgents_no=40, Max_iter=5000,
    disp=false, seed=3, ftol=1e-14, patience=60)
@printf("  stopped early, best = %.4e (curve is padded to %d entries)\n",
    r[1], length(r[3]))

# a callback can stop the run and record whatever you like
trace = Float64[]
SOMPA(fobj=sphere, lb=lb_so, ub=ub_so, SearchAgents_no=40, Max_iter=1000,
    disp=false, seed=3,
    callback=(iter, best, pos, curve) -> begin
        push!(trace, best)
        return best > 1e-12          # false ends the run
    end)
@printf("  callback trace length = %d, final = %.4e\n", length(trace), trace[end])

# ---------------------------------------------------------------------------
println("\n=== 4. Batch (vectorised) objective ===")
# One call per generation instead of one per agent -- the place to plug in a
# GPU kernel, a BLAS-heavy model, or a solver with batch mode.
batch_sphere(X::Matrix{Float64}) = vec(sum(abs2, X, dims=2))
r = SOMPA(fobj=sphere, fobj_batch=batch_sphere, lb=lb_so, ub=ub_so,
    SearchAgents_no=200, Max_iter=100, disp=false, seed=11)
@printf("  batch result = %.4e\n", r[1])

# ---------------------------------------------------------------------------
println("\n=== 5. NMPA variant on a 20-D sphere ===")
for v in (:standard_mpa, :nmpa)
    res = SOMPA(fobj=sphere, lb=lb_so, ub=ub_so, SearchAgents_no=40,
        Max_iter=300, disp=false, seed=5, variant=v)
    @printf("  %-14s -> %.4e\n", v, res[1])
end

# ---------------------------------------------------------------------------
println("\n=== 6. Multi objective: archive modes and quality indicators ===")
println("    (mean over 15 seeds -- a single run is far too noisy to compare)")
ref = [1.1, 1.1]
for (am, es) in ((:pareto, :crowding), (:pareto, :random),
    (:fronts, :crowding), (:fronts, :random))
    hvs = Float64[]
    sps = Float64[]
    nd = true
    for s in 1:15
        AP, AO, _ = MOMPA(fobj=zdt1, lb=lb_mo, ub=ub_mo, SearchAgents_no=100,
            Max_iter=200, num_objectives=2, disp=false, seed=s,
            archive_mode=am, elite_selection=es)
        push!(hvs, hypervolume(AO, ref))
        push!(sps, spacing_metric(AO))
        nd &= all(pareto_filter(AO))
    end
    @printf("  archive=%-7s elite=%-9s HV=%.5f  SP=%.4f  all non-dominated=%s\n",
        am, es, sum(hvs) / length(hvs), sum(sps) / length(sps), nd)
end

# true front of ZDT1, for IGD
f1 = collect(range(0, 1, length=500))
true_front = hcat(f1, 1 .- sqrt.(f1))
AP, AO, _ = MOMPA(fobj=zdt1, lb=lb_mo, ub=ub_mo, SearchAgents_no=100,
    Max_iter=200, num_objectives=2, disp=false, seed=17)
@printf("  IGD = %.5f   GD = %.5f\n", igd(AO, true_front), gd(AO, true_front))

# ---------------------------------------------------------------------------
println("\n=== 7. Merging independent runs into one front ===")
allO = reduce(vcat, (MOMPA(fobj=zdt1, lb=lb_mo, ub=ub_mo, SearchAgents_no=60,
    Max_iter=100, num_objectives=2, disp=false, seed=s)[2] for s in 1:5))
merged = allO[pareto_filter(allO), :]
@printf("  %d solutions from 5 runs -> %d non-dominated, HV = %.5f\n",
    size(allO, 1), size(merged, 1), hypervolume(merged, ref))

# ---------------------------------------------------------------------------
println("\n=== 8. History with compression ===")
mktempdir() do dir
    path = joinpath(dir, "history.h5")
    _, _, curve = SOMPA(fobj=sphere, lb=lb_so, ub=ub_so, SearchAgents_no=40,
        Max_iter=100, disp=false, seed=1,
        saveHDF=true, hdf_filepath=path, hdf_compress=6,
        history_save_interval=50)
    pos, fit, ismo, _, _ = ReadMPAHistory(path)
    @printf("  positions %s   fitness %s   multi-objective = %s   file = %.1f KiB\n",
        size(pos), size(fit), ismo, filesize(path) / 1024)
    @printf("  stored convergence curve matches: %s\n",
        ReadMPAConvergence(path) ≈ curve)
end

# ---------------------------------------------------------------------------
println("\n=== 9. Legacy 0.1 API still works ===")
fit, pos, CV = MPA(40, 100, [], lb_so, ub_so, D, sphere; disp=false)
@printf("  MPA(...) -> %.4e, curve size %s\n", fit, size(CV))
ci = confidence_interval(pos[1:3], p -> sum(abs2, p .- 1.0), 0.95)
println("  confidence_interval sample:\n", ci)
