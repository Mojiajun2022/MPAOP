# Run with:   mpiexec -n 4 julia --project=. test/mpi_test.jl
#
# Checks that the MPI path returns exactly the same answer as a serial run with
# the same seed (all stochastic decisions are taken on rank 0 and broadcast, so
# the result must not depend on the number of ranks).

using MPAOP, MPI, Test

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NP = MPI.Comm_size(COMM)

sphere(x::Vector{Float64}) = sum(abs2, x)
function zdt1(x::Vector{Float64})
    f1 = x[1]
    g = 1.0 + (9.0 / (length(x) - 1)) * sum(@view x[2:end])
    return [f1, g * (1.0 - sqrt(f1 / g))]
end

d = 8
lb = fill(-5.0, d)
ub = fill(5.0, d)
# 17 agents on purpose: not divisible by the rank count -- the new Allgatherv
# split handles it without padding the population.
N = 17

ref_fit, ref_pos, ref_curve = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=N,
    Max_iter=40, disp=false, seed=1234, parallelism=:serial)

mpi_fit, mpi_pos, mpi_curve = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=N,
    Max_iter=40, disp=false, seed=1234, parallelism=:mpi)

hyb_fit, hyb_pos, _ = SOMPA(fobj=sphere, lb=lb, ub=ub, SearchAgents_no=N,
    Max_iter=40, disp=false, seed=1234, parallelism=:mpi_threads)

mo_ref = MOMPA(fobj=zdt1, lb=zeros(d), ub=ones(d), SearchAgents_no=N, Max_iter=30,
    num_objectives=2, disp=false, seed=99, parallelism=:serial)
mo_mpi = MOMPA(fobj=zdt1, lb=zeros(d), ub=ones(d), SearchAgents_no=N, Max_iter=30,
    num_objectives=2, disp=false, seed=99, parallelism=:mpi)

@testset "MPI rank $RANK / $NP" begin
    @test mpi_fit == ref_fit
    @test mpi_pos == ref_pos
    @test mpi_curve == ref_curve
    @test hyb_fit == ref_fit
    @test hyb_pos == ref_pos
    if RANK == 0
        @test mo_mpi[1] == mo_ref[1]
        @test mo_mpi[2] == mo_ref[2]
        @test mo_mpi[3] == mo_ref[3]
    else
        @test mo_mpi[1] === nothing
    end
end

MPI.Barrier(COMM)
RANK == 0 && println("MPI test passed on $NP ranks (population $N, not divisible by $NP).")
MPI.Finalize()
