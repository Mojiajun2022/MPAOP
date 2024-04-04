using MPAOP
using Optim
using Plots
using Statistics
using MPI
MPI.Initialized() || MPI.Init()
root = 0
comm = MPI.COMM_WORLD
Size = MPI.Comm_size(MPI.COMM_WORLD)
rank = MPI.Comm_rank(comm)
function fobj(x)
    f1 = abs(x[1] + x[2]) - abs(x[3])
    f2 = x[1] * x[2] * x[3] + 18
    f3 = x[1]^2 * x[2] + 3 * x[3]
    f = abs(f1) + abs(f2) + abs(f3)
    return f
end
# p0 = [0.01; 0.01; 0.01]
# p0 = rand(3,1)
ntest = 100
F = zeros(Float64, 1, ntest)
narvs = 3
lb = -10.0 .* [1.0, 1.0, 1.0]
ub = 10.0 .* [1.0, 1.0, 1.0]

# @time begin
    for i = 1:ntest
    fval, x, CV = MPA_MPI(200, 20,[],lb, ub,  length(ub), fobj;FADs0 = 0.2,P0=0.5,disp= false, Write = false)
    if rank == root
    f1 = abs(x[1] + x[2]) - abs(x[3])
    f2 = x[1] * x[2] * x[3] + 18
    f3 = x[1]^2 * x[2] + 3 * x[3]
    F[i] = abs(f1) + abs(f2) + abs(f3)
    else
        break

    # f1 = Plots.plot([1:length(CV);],vec(CV))
    # println(   minimum(CV))
    # savefig(f1,"1.png")
    end
    end
# end
if rank == root
println(string("平均结果为", mean(F)))
end 
MPI.Finalize()