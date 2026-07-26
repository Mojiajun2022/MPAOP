# test_mpi_ranks.jl
using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

println("Hello from rank $(rank) out of $(size) processes.")

# 每个进程基于自己的rank做一点不同的事情
for i = 1:3
    if MPI.Comm_rank(comm) == i-1 # 模拟不同进程的特定任务
        println("Greetings specifically from rank $(MPI.Comm_rank(comm))!")
    end
    MPI.Barrier(comm) # 等待所有进程到达这一点
end

if rank == 0
    println("Rank 0: MPI test seems to be working if you see different ranks above.")
end

MPI.Finalize()