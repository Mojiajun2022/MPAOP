using MPAOP
using Test
using Random
using LinearAlgebra

# ---------------------------------------------------------------- objectives
sphere(x::Vector{Float64}) = sum(abs2, x)
rosenbrock(x::Vector{Float64}) = sum(100.0 * (x[i+1] - x[i]^2)^2 + (1 - x[i])^2 for i in 1:length(x)-1)
rastrigin(x::Vector{Float64}) = 10 * length(x) + sum(xi^2 - 10 * cos(2π * xi) for xi in x)

function zdt1(x::Vector{Float64})
    f1 = x[1]
    g = 1.0 + (9.0 / (length(x) - 1)) * sum(@view x[2:end])
    return [f1, g * (1.0 - sqrt(f1 / g))]
end

function zdt3(x::Vector{Float64})
    f1 = x[1]
    g = 1.0 + (9.0 / (length(x) - 1)) * sum(@view x[2:end])
    h = 1.0 - sqrt(f1 / g) - (f1 / g) * sin(10π * f1)
    return [f1, g * h]
end

function dtlz2(x::Vector{Float64})
    M = 3
    k = length(x) - M + 1
    g = sum((x[i] - 0.5)^2 for i in (length(x)-k+1):length(x))
    f = ones(M)
    for i in 1:M
        f[i] *= (1 + g)
        for j in 1:(M-i)
            f[i] *= cos(x[j] * π / 2)
        end
        i > 1 && (f[i] *= sin(x[M-i+1] * π / 2))
    end
    return f
end

# --------------------------------------------------------------------- tests
@testset "MPAOP.jl" begin

    @testset "kernels" begin
        X = [3.0 -4.0; 0.5 9.0]
        lb = [0.0, 1.0]
        ub = [1.0, 5.0]
        MPAOP.clamp_cols!(X, lb, ub)
        @test X == [1.0 0.0; 1.0 5.0]

        rng = Xoshiro(1)
        Z = zeros(50, 4)
        MPAOP.levy!(Z, rng)
        @test all(isfinite, Z)
        @test size(levy(7, 3, 1.5)) == (7, 3)

        @test MPAOP.chunk_ranges(10, 3) == [1:4, 5:7, 8:10]
        @test MPAOP.chunk_ranges(2, 5) == [1:1, 2:2]
        c, d = MPAOP.block_counts(10, 4)
        @test sum(c) == 10 && d[1] == 0 && all(diff(Int.(d)) .== Int.(c[1:end-1]))

        A = reshape(collect(1.0:6.0), 2, 3)
        B = Matrix{Float64}(undef, 3, 2)
        MPAOP.transpose_into!(B, A)
        @test B == permutedims(A)
    end

    @testset "initialisation" begin
        rng = Xoshiro(42)
        lb = [-1.0, 2.0, 0.0]
        ub = [1.0, 3.0, 10.0]
        X = Matrix{Float64}(undef, 3, 25)
        for m in (:uniform, :lhs, :center)
            MPAOP.init_population!(X, lb, ub, [], rng; method=m)
            @test all(lb[i] <= X[i, j] <= ub[i] for i in 1:3, j in 1:25)
        end
        # LHS must use every stratum exactly once
        MPAOP.init_population!(X, lb, ub, [], rng; method=:lhs)
        bins = [floor(Int, (X[1, j] - lb[1]) / (ub[1] - lb[1]) * 25) for j in 1:25]
        @test sort(bins) == collect(0:24)

        MPAOP.init_population!(X, lb, ub, [0.5, 2.5, 5.0], rng)
        @test X[:, 1] == [0.5, 2.5, 5.0]

        P = initialization(12, 3, ub, lb, [])
        @test size(P) == (12, 3)

        Xo = similar(X)
        MPAOP.opposite_population!(Xo, X, lb, ub)
        @test Xo ≈ (lb .+ ub) .- X
    end

    @testset "dominance / sorting" begin
        @test is_dominated([2.0, 2.0], [1.0, 1.0])        # b dominates a
        @test is_dominated([2.0, 2.0], [1.0, 2.0])
        @test !is_dominated([1.0, 1.0], [2.0, 2.0])
        @test !is_dominated([1.0, 2.0], [2.0, 1.0])
        @test !is_dominated([1.0, 1.0], [1.0, 1.0])       # equal, no dominance

        rng = Xoshiro(7)
        for trial in 1:20
            n = rand(rng, 3:40)
            M = rand(rng, 2:4)
            O = round.(rand(rng, n, M), digits=2)
            ranks, fronts = non_dominated_sort(O)

            # brute-force reference ranks
            ref = zeros(Int, n)
            remaining = Set(1:n)
            r = 0
            while !isempty(remaining)
                r += 1
                cur = [i for i in remaining
                       if !any(is_dominated(O[i, :], O[j, :]) for j in remaining if j != i)]
                for i in cur
                    ref[i] = r
                    delete!(remaining, i)
                end
            end
            @test ranks == ref
            @test sort(vcat(fronts...)) == collect(1:n)
            @test all(all(ranks[i] == k for i in fronts[k]) for k in eachindex(fronts))
        end
    end

    @testset "crowding distance" begin
        O = [0.0 1.0; 0.25 0.75; 0.5 0.5; 0.75 0.25; 1.0 0.0]
        cd = zeros(5)
        calculate_crowding_distance!(cd, O, collect(1:5), 2)
        @test cd[1] == Inf && cd[5] == Inf
        @test all(isfinite, cd[2:4])
        @test cd[3] ≈ 1.0                      # evenly spaced front
        cd2 = zeros(2)
        calculate_crowding_distance!(cd2, O[1:2, :], [1, 2], 2)
        @test all(cd2 .== Inf)
    end

    @testset "archive" begin
        rng = Xoshiro(3)
        dim, M, n = 4, 2, 30
        X = rand(rng, dim, n)
        F = rand(rng, M, n)

        A = MPAOP.MOArchive(dim, M, 10, n; mode=:pareto)
        MPAOP.update!(A, X, F)
        @test A.n <= 10
        P, O = MPAOP.archive_matrices(A)
        @test size(P) == (A.n, dim) && size(O) == (A.n, M)
        # every archived point must be non-dominated within the archive
        @test all(pareto_filter(O))

        # and non-dominated w.r.t. the whole population it saw
        Fall = permutedims(F)
        for i in axes(O, 1)
            @test !any(is_dominated(O[i, :], Fall[j, :]) for j in axes(Fall, 1))
        end

        # legacy fronts mode fills up to capacity
        A2 = MPAOP.MOArchive(dim, M, 10, n; mode=:fronts)
        MPAOP.update!(A2, X, F)
        @test A2.n == 10

        # public wrapper round trip
        ap, ao = update_archive(zeros(0, dim), zeros(0, M),
            permutedims(X), permutedims(F), 10, dim, M)
        @test size(ap, 1) == size(ao, 1) <= 10
        @test all(pareto_filter(ao))

        E = select_elite_from_archive(ap, ao, 7, dim)
        @test size(E) == (7, dim)
        @test all(any(all(E[i, :] .== ap[j, :]) for j in axes(ap, 1)) for i in 1:7)
        E2 = select_elite_from_archive(ap, ao, 7, dim; strategy=:crowding, rng=Xoshiro(1))
        @test size(E2) == (7, dim)
        @test size(select_elite_from_archive(zeros(0, dim), zeros(0, M), 5, dim)) == (5, dim)
    end

    @testset "metrics" begin
        # unit square, one point at origin -> hv = 1
        @test hypervolume([0.0 0.0], [1.0, 1.0]) ≈ 1.0
        @test hypervolume([0.0 1.0; 1.0 0.0], [2.0, 2.0]) ≈ 3.0
        @test hypervolume([0.0 0.0 0.0], [1.0, 1.0, 1.0]) ≈ 1.0
        @test hypervolume([2.0 2.0], [1.0, 1.0]) == 0.0     # nothing dominates ref
        # MC estimate close to exact
        pts = [0.1 0.9; 0.4 0.5; 0.8 0.2]
        hv_exact = hypervolume(pts, [1.0, 1.0]; method=:exact)
        hv_mc = hypervolume(pts, [1.0, 1.0]; method=:mc, samples=400_000, rng=Xoshiro(9))
        @test isapprox(hv_exact, hv_mc; rtol=0.05)

        front = [0.0 1.0; 0.5 0.5; 1.0 0.0]
        @test igd(front, front) ≈ 0.0
        @test gd(front, front) ≈ 0.0
        @test igd([0.0 1.0], front) > 0
        @test spacing_metric(front) ≈ 0.0 atol = 1e-12
        @test max_spread(front) ≈ sqrt(2)
        @test pareto_filter([0.0 0.0; 1.0 1.0]) == [true, false]
    end

    @testset "SOMPA convergence" begin
        for (f, lo, hi, d, tol) in ((sphere, -10.0, 10.0, 10, 1e-4),
            (rastrigin, -5.12, 5.12, 5, 5.0),
            (rosenbrock, -5.0, 10.0, 5, 5.0))
            fit, pos, curve = SOMPA(fobj=f, lb=fill(lo, d), ub=fill(hi, d),
                SearchAgents_no=40, Max_iter=150, disp=false, seed=2024)
            @test length(pos) == d
            @test length(curve) == 150
            @test fit <= tol
            @test fit ≈ f(pos)
            @test issorted(curve, rev=true)          # best-so-far is monotone
            @test all(fill(lo, d) .<= pos .<= fill(hi, d))
        end
    end

    @testset "SOMPA options" begin
        d = 6
        lb = fill(-5.0, d)
        ub = fill(5.0, d)

        # reproducibility
        a = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=30, disp=false, seed=99)
        b = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=30, disp=false, seed=99)
        @test a[1] == b[1] && a[2] == b[2] && a[3] == b[3]
        c = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=30, disp=false, seed=100)
        @test a[1] != c[1]

        # threads must give exactly the same answer as serial for a given seed
        t = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=30,
            disp=false, seed=99, parallelism=:threads)
        @test t[1] == a[1] && t[2] == a[2]

        # batch objective == scalar objective
        batch(X) = [sphere(collect(view(X, i, :))) for i in axes(X, 1)]
        bb = SOMPA(fobj=sphere, fobj_batch=batch, lb=lb, ub=ub,
            SearchAgents_no=20, Max_iter=30, disp=false, seed=99)
        @test bb[1] ≈ a[1] && bb[2] ≈ a[2]

        # p0 seeding
        p0 = fill(0.25, d)
        f0, _, _ = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=8, Max_iter=1,
            p0_optional=p0, disp=false, seed=5)
        @test f0 <= sphere(p0)

        # variants, init schemes, opposition
        for v in (:standard_mpa, :nmpa), ini in (:uniform, :lhs, :center)
            r = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=40,
                disp=false, seed=7, variant=v, init=ini)
            @test isfinite(r[1])
        end
        ro = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=40,
            disp=false, seed=7, opposition=true)
        @test isfinite(ro[1])

        # Fixbox = false lets the search leave the box
        rf = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=20,
            disp=false, seed=7, Fixbox=false)
        @test isfinite(rf[1])

        # early stopping
        stopped = Ref(0)
        cbr = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=500,
            disp=false, seed=7,
            callback=(it, bf, bp, cv) -> (stopped[] = it; it < 12))
        @test stopped[] == 12
        @test length(cbr[3]) == 500                    # curve is padded, not truncated
        pat = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=2000,
            disp=false, seed=7, ftol=1e-12, patience=5)
        @test isfinite(pat[1])

        # legacy rounding on request
        rr = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=20,
            disp=false, seed=7, round_digits=4)
        @test rr[1] == round(rr[1], digits=4)

        # NaN from the objective must not break the search
        nanobj(x) = x[1] > 0 ? NaN : sphere(x)
        rn = SOMPA(fobj=nanobj, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=20,
            disp=false, seed=7)
        @test isfinite(rn[1])

        # errors
        @test_throws DimensionMismatch SOMPA(fobj=sphere, lb=lb, ub=ub[1:end-1],
            SearchAgents_no=5, Max_iter=5, disp=false)
        @test_throws ArgumentError SOMPA(fobj=sphere, lb=lb, ub=ub,
            SearchAgents_no=5, Max_iter=5, disp=false, variant=:nope)
    end

    @testset "MOMPA" begin
        d = 10
        lb = zeros(d)
        ub = ones(d)

        AP, AO, curve = MOMPA(fobj=zdt1, lb=lb, ub=ub, SearchAgents_no=60,
            Max_iter=120, num_objectives=2, disp=false, seed=11)
        @test size(AP, 2) == d
        @test size(AO, 2) == 2
        @test size(AP, 1) == size(AO, 1) == length(pareto_filter(AO))
        @test all(pareto_filter(AO))                       # no dominated member
        @test length(curve) == 120
        @test all(0.0 .<= AP .<= 1.0)
        # objectives are consistent with the positions
        for i in axes(AP, 1)
            @test AO[i, :] ≈ zdt1(collect(AP[i, :]))
        end
        # converged reasonably close to the analytic front f2 = 1 - sqrt(f1)
        err = sum(abs(AO[i, 2] - (1 - sqrt(AO[i, 1]))) for i in axes(AO, 1)) / size(AO, 1)
        @test err < 0.05

        # reproducibility
        r1 = MOMPA(fobj=zdt1, lb=lb, ub=ub, SearchAgents_no=30, Max_iter=25,
            num_objectives=2, disp=false, seed=4)
        r2 = MOMPA(fobj=zdt1, lb=lb, ub=ub, SearchAgents_no=30, Max_iter=25,
            num_objectives=2, disp=false, seed=4)
        @test r1[1] == r2[1] && r1[2] == r2[2]

        # legacy archive/elite behaviour still available
        l1, l2, _ = MOMPA(fobj=zdt1, lb=lb, ub=ub, SearchAgents_no=30, Max_iter=40,
            num_objectives=2, disp=false, seed=4,
            archive_mode=:fronts, elite_selection=:random)
        @test size(l1, 1) == 30

        # three objectives
        t1, t2, _ = MOMPA(fobj=dtlz2, lb=zeros(8), ub=ones(8), SearchAgents_no=40,
            Max_iter=50, num_objectives=3, disp=false, seed=6)
        @test size(t2, 2) == 3
        @test all(pareto_filter(t2))

        # disconnected front
        z3p, z3o, _ = MOMPA(fobj=zdt3, lb=lb, ub=ub, SearchAgents_no=50, Max_iter=80,
            num_objectives=2, disp=false, seed=8)
        @test all(pareto_filter(z3o))

        # threads == serial for the same seed
        s = MOMPA(fobj=zdt1, lb=lb, ub=ub, SearchAgents_no=30, Max_iter=25,
            num_objectives=2, disp=false, seed=4, parallelism=:threads)
        @test s[2] == r1[2]

        # batch objective
        mbatch(X) = permutedims(reduce(hcat, [zdt1(collect(view(X, i, :))) for i in axes(X, 1)]))
        bt = MOMPA(fobj=zdt1, fobj_batch=mbatch, lb=lb, ub=ub, SearchAgents_no=30,
            Max_iter=25, num_objectives=2, disp=false, seed=4)
        @test bt[2] ≈ r1[2]

        # hypervolume reporting path + callback
        seen = Ref(0)
        MOMPA(fobj=zdt1, lb=lb, ub=ub, SearchAgents_no=20, Max_iter=10,
            num_objectives=2, disp=false, seed=4, hv_ref=[1.5, 1.5],
            callback=(it, P, O, cv) -> (seen[] = it; true))
        @test seen[] == 10

        @test_throws ArgumentError MOMPA(fobj=zdt1, lb=lb, ub=ub, SearchAgents_no=10,
            Max_iter=5, num_objectives=1, disp=false)
    end

    @testset "HDF5 round trip" begin
        mktempdir() do dir
            so = joinpath(dir, "so.h5")
            fit, pos, curve = SOMPA(fobj=sphere, lb=fill(-5.0, 4), ub=fill(5.0, 4),
                SearchAgents_no=12, Max_iter=10, disp=false, seed=1,
                saveHDF=true, hdf_filepath=so, history_save_interval=5)
            @test isfile(so)
            P, F, ismo, ap, ao = ReadMPAHistory(so)
            @test size(P) == (12, 4, 20)
            @test size(F) == (12, 20)
            @test ismo == false
            @test ReadMPAConvergence(so) ≈ curve

            mo = joinpath(dir, "mo.h5")
            AP, AO, _ = MOMPA(fobj=zdt1, lb=zeros(6), ub=ones(6), SearchAgents_no=12,
                Max_iter=10, num_objectives=2, disp=false, seed=1,
                saveHDF=true, hdf_filepath=mo, history_save_interval=10, hdf_compress=3)
            P2, O2, ismo2, ap2, ao2 = ReadMPAHistory(mo)
            @test ismo2 == true
            @test size(P2) == (12, 6, 20)
            @test size(O2) == (12, 2, 20)
            @test ap2 ≈ AP && ao2 ≈ AO

            # early-stopped run must still produce a consistent file
            es = joinpath(dir, "es.h5")
            SOMPA(fobj=sphere, lb=fill(-5.0, 3), ub=fill(5.0, 3), SearchAgents_no=8,
                Max_iter=50, disp=false, seed=1, saveHDF=true, hdf_filepath=es,
                callback=(it, a, b, c) -> it < 7)
            P3, F3, _, _, _ = ReadMPAHistory(es)
            @test size(P3, 3) == 14

            @test ReadMPAHistory(joinpath(dir, "missing.h5"))[1] === nothing
        end
    end

    @testset "CSV log" begin
        mktempdir() do dir
            p = joinpath(dir, "log.csv")
            SOMPA(fobj=sphere, lb=fill(-1.0, 3), ub=fill(1.0, 3), SearchAgents_no=8,
                Max_iter=6, disp=false, seed=1, write_csv_log=true, csv_log_filepath=p)
            @test isfile(p)
            txt = read(p, String)
            @test occursin("start running", txt)
            @test occursin("Iter: 6", txt)
            @test occursin("finish", txt)

            q = joinpath(dir, "mo.csv")
            MOMPA(fobj=zdt1, lb=zeros(4), ub=ones(4), SearchAgents_no=8, Max_iter=4,
                num_objectives=2, disp=false, seed=1, write_csv_log=true,
                csv_log_filepath=q, csv_flush_every=2)
            @test occursin("Objective 2", read(q, String))
        end
    end

    @testset "legacy 0.1 API" begin
        d = 5
        lb = fill(-10.0, d)
        ub = fill(10.0, d)
        fit, pos, CV = MPA(30, 80, [], lb, ub, d, sphere; disp=false)
        @test size(CV) == (1, 80)          # 0.1 returned a row matrix
        @test length(pos) == d
        @test fit < 1e-3
        fit2, _, _ = MPA(20, 40, zeros(d), lb, ub, d, sphere; disp=false,
            Threads_parallel=true, FADs0=0.3, P0=0.6)
        @test isfinite(fit2)

        @test size(initialization(9, 3, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0])) == (9, 3)

        # confidence interval of a quadratic:  f = Σ (x-μ)²/σ² has H = 2/σ² I
        μ = [1.0, -2.0]
        σ = [0.5, 2.0]
        quad(x) = sum(((x .- μ) ./ σ) .^ 2)
        ci = confidence_interval(copy(μ), quad, 0.95)
        @test size(ci) == (2, 2)
        @test all(ci[:, 1] .< μ .< ci[:, 2])
        # analytic: half width = 1.959964 * σ/√2
        @test (ci[1, 2] - ci[1, 1]) / 2 ≈ 1.959963985 * σ[1] / sqrt(2) rtol = 1e-3
        @test (ci[2, 2] - ci[2, 1]) / 2 ≈ 1.959963985 * σ[2] / sqrt(2) rtol = 1e-3
        @test_throws ArgumentError confidence_interval(copy(μ), quad, 1.5)
    end

    @testset "steady-state allocations" begin
        lb = fill(-5.0, 8)
        ub = fill(5.0, 8)
        SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=16, Max_iter=2, disp=false, seed=1)
        a = @allocated SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=16,
            Max_iter=200, disp=false, seed=1)
        # the whole 200-iteration run must stay far below the ~200 MB the
        # pre-0.3 implementation needed for the same work
        @test a < 2_000_000

        MOMPA(fobj=zdt1, lb=zeros(8), ub=ones(8), SearchAgents_no=16, Max_iter=2,
            num_objectives=2, disp=false, seed=1)
        m = @allocated MOMPA(fobj=zdt1, lb=zeros(8), ub=ones(8), SearchAgents_no=16,
            Max_iter=200, num_objectives=2, disp=false, seed=1)
        @test m < 20_000_000
    end
end
