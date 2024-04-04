

# MPAOP.jl

The package MPAOP.jl provides a global optimization algorithm——*Marine Predators Algorithm*，*MPA* for Julia. 

The Marine Predator Algorithm(MPA) is first proposed by Afshin Faramarzi et al.(Faramarzi, Afshin, et al. Expert systems with applications 152 (2020): 113377).

The original version is the MATLAB version. This project converts the Matlab version to the Julia version.

There are some corresponding improvements have been made, such as it can start from initial values and it is ability to run MPI parallel calculations, etc.

## Description

We can use MPAOP.jl to search the global optimal solution. The function form is: 

```julia
Top_predator_fit, Top_predator_pos, CV = MPAOP(SearchAgents_no, Max_iterations, p0, lb, ub, narvs, fobj;disp=true, Fixbox=true, Write=false, FADs0=0.2, P0=0.5)
```

###### Output results：

`Top_predator_fit`:  The minimum value searched by MPA.

`Top_predator_fit`: The parameter corresponding to the minimum value.

`CV`: Iterative evolution change curve.

###### Input parameters：

`SearchAgents_no`: Number of search agents. After trying it, I think it is most suitable to set it between 20 and 80.

` Max_iterations: `Maximum number of iterations. After trying, I think setting SearchAgents_no * Max_iterations=10000 or more can get better results.

` p0`: Starting value when searching, If there is no definite starting value, you can enter ` p0=[]` as the starting value, so that the starting value will be randomly generated.

`lb`: Lower bound on fitting parameters.

`ub`: Upper bound on fitting parameters.

`narvs`: Number of fitting parameters, if nothing else happens, it should be equal to ` length(lb)`.

`fobj`: Fitted objective function. The output value must be Float, not a Vector or Matrix.

`disp`: This item is an optional value, indicating that the best parameters and best error values obtained at each iteration are output in the run box. The default is `disp = true`, if you do not want to display the output, you can modify it to `disp = false`.

`Fixbox` \: If you do not want to constrain the boundaries of this fitting, that is, do not constrain the fitting parameters between lb and ub, then set `Fixbox=false`, and the default is `Fixbox=true`.

`Write`: If Write=true, the output of each iteration will be stored in a file named "fitting_process", which will record the time of each run and the results of each iteration, and the default is `Write = false`.

`FADs0`: Fish Aggregating Devices, it will affect the algorithm optimization process is usually taken as 0.2, If the fitting needs to be more inclined to local search, increase it, otherwise decrease it, and the range of FADs0 is 0 to 1. The default is `FADs0 = 0.2`.

`P0`: A constant used to regulate predator behavior, usually between 0 and 1.  If the fitting needs to be more inclined to local search, decrease it, otherwise increase it. The default is `P0 = 0.5`.

## Example

###### Run in serial

If we construct an objective function for searching, such as

```julia
function fobj(x) # The function need to be fitted
    f1 = abs(x[1] + x[2]) - abs(x[3])
    f2 = x[1] * x[2] * x[3] + 18
    f3 = x[1]^2 * x[2] + 3 * x[3]
    f = abs(f1) + abs(f2) + abs(f3)
    return f
end
```

The above function has a lot of local optimal solutions. Then we can use the following code for searching for the minimum value. 

```julia
using MPAOP
p0 = [] #Random initial value
SearchAgents_no = 48 #Number of search agents（Similar to the number of birds in the particle swarm algorithm, the general range is set to between 20 and 80.
Max_iterations = 200 # The total number of iterations, that is, the number of times all agents have finished running.
#If using mpi, use the MPA_MPI function, and use the mpirun -np N julia test.jl command directly in the terminal.

lb = -10.0 .* [1.0, 1.0, 1.0] #fitting lower bound
ub = 10.0 .* [1.0, 1.0, 1.0] #fitting upper bound
narvs = length(lb)# Number of fitting parameters

Top_predator_fit, Top_predator_pos, CV = MPA(SearchAgents_no, Max_iterations, p0, lb, ub, narvs, fobj; disp=false)
print("The best value from the search is ", Top_predator_fit)
print("The best parameters from the search are ", Top_predator_pos)

using Plots
plot(CV', xlabel="Iterations", ylabel="Best value", grid=nothing, leg=false, framestyle=:box)
```

Output:

 ```julia
 The best value from the search is 0.0
 The best parameters from the search are [1.421, -4.34, 2.92]
 ```

<img width="956" alt="image" src="https://github.com/Mojiajun2022/Optimization-algorithm/assets/130334983/0c2bfaaf-c983-48e0-887f-140e5a3a05c7">



###### Run in parallel

When you need to use parallel methods to perform operations on each "Search agent", you can call MPI to run, the The function to run MPI is MPA_MPI. There is no need to modify any part of the above code during execution. You only need to modify MPA to MPA_MPI, and add MPI initialization conditions. If you run the program fit.jl file is as follows:

```julia
using MPI, MPAOP, Plots
MPI.Initialized() || MPI.Init()
root = 0
comm = MPI.COMM_WORLD
Size = MPI.Comm_size(MPI.COMM_WORLD)
rank = MPI.Comm_rank(comm)
function fobj(x) # The function need to be fitted
    f1 = abs(x[1] + x[2]) - abs(x[3])
    f2 = x[1] * x[2] * x[3] + 18
    f3 = x[1]^2 * x[2] + 3 * x[3]
    f = abs(f1) + abs(f2) + abs(f3)
    return f
end

using MPAOP
p0 = [] #Random initial value
SearchAgents_no = 48 #Number of search agents（Similar to the number of birds in the particle swarm algorithm, the general range is set to between 20 and 80.
Max_iterations = 200 # The total number of iterations, that is, the number of times all agents have finished running.
#If using mpi, use the MPA_MPI function, and use the mpirun -np N julia test.jl command directly in the terminal.

lb = -10.0 .* [1.0, 1.0, 1.0] #fitting lower bound
ub = 10.0 .* [1.0, 1.0, 1.0] #fitting upper bound
narvs = length(lb)# Number of fitting parameters

Top_predator_fit, Top_predator_pos, CV = MPA_MPI(SearchAgents_no, Max_iterations, p0, lb, ub, narvs, fobj; disp=true, Write=true)
if rank == root
  print("The best value from the search is ", Top_predator_fit)
print("The best parameters from the search are ", Top_predator_pos)
p = plot(CV', xlabel="Iterations", ylabel="Best value", grid=nothing, leg=false, framestyle=:box)
  savefig(p, "CV.png")
end
```

Then you can run the above code in terminal, i.e.,

```
mpirun -np N julia fit.jl
```

where N is the number of cores for parallel processing. The output of each iteration will be displayed on the terminal if `disp=true`, and the process of each iteration will be saved in the "fitting_process" file if `write=true`.



###### Pre-test

Due to randomness, the above results are different every time it is run, so sometimes it is necessary to test the stability of the parameters to determine the SearchAgents_no and Max_iterations that need to be set. A test method is provided below：

```julia
using MPAOP,Statistics

function fobj(x) # The function to be fitted
    f1 = abs(x[1] + x[2]) - abs(x[3])
    f2 = x[1] * x[2] * x[3] + 18
    f3 = x[1]^2 * x[2] + 3 * x[3]
    f = abs(f1) + abs(f2) + abs(f3)
    return f
end

ntest = 500 #Set the number of repeated runs. The more times, the more accurate the statistical results will be.

F = zeros(Float64, 1, ntest)
lb = -10.0 .* [1.0, 1.0, 1.0] #fitting lower bound
ub = 10.0 .* [1.0, 1.0, 1.0] #fitting upper bound
narvs = length(lb)# Number of fitting parameters
for i in 1:ntest
    p0 = [0, 0, 0] #initial value
    # p0 = [] #If you do not specify an initial value, enter []
    SearchAgents_no = 36 #Number of search agents（Similar to the number of birds in the particle swarm algorithm, the general range is set to between 20 and 80.
    Max_iterations = 200 # The total number of iterations, that is, the number of times all agents have finished running.
    #If using mpi, use the MPAF_MPI function, and use the mpirun -np N julia test.jl command directly in the terminal.
    fval, x, CV = MPA(SearchAgents_no, Max_iterations, p0, lb, ub, narvs, fobj; disp=false, Write=false, FADs0=0.3, P0=0.6) #if disp = true, the best result and best parameters obtained until the current iteration will be displayed during each iteration, and when Write = true, if Wirite=true, the corresponding fitting process will be saved in the "fitting_process" file. 
    F[i] = fval
end

println(string("The mean of results by using MPAOP is", mean(F)))

```

By setting different `SearchAgents_no`,` Max_iterations`, `FADs0`,` P0` and other parameters, compare different `mean(F)` to determine the most appropriate parameters.



###### 





