using MPAOP
using Optim
using Plots
using Statistics

#Compare with optimize

function fobj(x) # The function to be fitted
    f1 = abs(x[1] + x[2]) - abs(x[3])
    f2 = x[1] * x[2] * x[3] + 18
    f3 = x[1]^2 * x[2] + 3 * x[3]
    f = abs(f1) + abs(f2) + abs(f3)
    return f
end
# p0 = [0.01; 0.01; 0.01]
# p0 = rand(3,1)
ntest = 500
F = zeros(Float64, 1, ntest)

lb = -10.0 .* [1.0, 1.0, 1.0] #fitting lower bound
ub = 10.0 .* [1.0, 1.0, 1.0] #fitting upper bound
narvs = length(lb)# Number of fitting parameters


# Fitting results by using MAPOP 
@time begin
    for i in 1:ntest
        p0 = [0, 0, 0] #initial value
        # p0 = [] #If you do not specify an initial value, enter []
        N1 = 36 #Number of search agents（Similar to the number of birds in the particle swarm algorithm, the general range is set to between 20 and 80.
        N2 = 200 # The total number of iterations, that is, the number of times all agents have finished running.
        #If using mpi, use the MPAF_MPI function, and use the mpirun -np N julia test.jl command directly in the terminal.
        fval, x, CV = MPA(N1, N2, p0, lb, ub, narvs, fobj, disp=false, Write=false) #if disp = true, the best result and best parameters obtained until the current iteration will be displayed during each iteration, and when Write = true, if Wirite=true, the corresponding fitting process will be saved in the "fitting_process" file. 
        # fval = optimize(fobj, lb, ub, p0, (NelderMead()), Optim.Options(g_tol=1e-8, iterations=10000)) #
        #  x = Optim.minimizer(fval)
        f1 = abs(x[1] + x[2]) - abs(x[3])
        f2 = x[1] * x[2] * x[3] + 18
        f3 = x[1]^2 * x[2] + 3 * x[3]
        F[i] = abs(f1) + abs(f2) + abs(f3)
    end
    Plots.plot(F')
    #plot(CV')
end
println(string("The mean of results by using MPAOP is", mean(F)))


# Fitting results by using NelderMead
@time begin
    for i in 1:ntest
        p0 = [0, 0, 0] #initial value
        #If using mpi, use the MPAF_MPI function, and use the mpirun -np N julia test.jl command directly in the terminal.
        fval = optimize(fobj, lb, ub, p0, (NelderMead()), Optim.Options(g_tol=1e-8, iterations=10000)) #
        x = Optim.minimizer(fval)
        f1 = abs(x[1] + x[2]) - abs(x[3])
        f2 = x[1] * x[2] * x[3] + 18
        f3 = x[1]^2 * x[2] + 3 * x[3]
        F[i] = abs(f1) + abs(f2) + abs(f3)
    end
    Plots.plot(F')
end
println(string("The mean of results by using NelderMead is", mean(F)))